# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2025 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
Options and literals for specifying and representing Multistrand simulations.
"""

import copy
from enum import IntEnum
from collections import namedtuple
from itertools import chain
from typing import List, Optional
import os.path
import importlib.resources

import numpy as np

from .interface import Interface
from ..utils.thermo import C2K
from ..objects import Strand, Complex, StopCondition
from ..__init__ import __version__


class Literals:
    """
    Preset tags that are used in the MergeResult objects (FirstStepRate,
    FirstStepLeakRate, FirstPassageRate) and multistrand.experiment
    """

    failure = "FAILURE"
    success = "SUCCESS"
    alt_success = "ALT_SUCCESS"

    """
     protected results may occupy the [x.tag in a x in options.interface.results]
    """
    time_out = "timeout"
    no_initial_moves = "noinitial"
    sim_error = "error"

    """ rate_method """
    metropolis = 1
    kawasaki = 2
    arrhenius = 3

    """ Nupack dangle options """
    dangles_none = 0
    dangles_some = 1
    dangles_all = 2

    """ Substrate type.    """
    substrateRNA = 1
    substrateDNA = 2

    """ Simulation modes """
    first_passage_time = 16  # 0x0010
    first_step = 48  # 0x0030
    trajectory = 128  # 0x0080
    transition = 256  # 0x0100

    """
    FD, May 8th, 2018:
    StopCondition and Macrostate definitions:
        Compuational expense, low to high: dissoc -- exact -- count -- loose
    """
    exact_macrostate = 0
    """match a secondary structure exactly (i.e. any system state that has a
    complex with this exact structure)"""
    bound_macrostate = 1
    """match any system state in which the given strand is bound to another
    strand"""
    dissoc_macrostate = 2
    ordered_macrostate = 2
    """match any system state in which there exists a complex with exactly the
    given strands, in that order (previous name: dissoc macrostate)"""
    loose_macrostate = 3
    """match a secondary structure, allowing a certain number of disagreements,
    allowing domain state to be unspecified (* notation)"""
    count_macrostate = 4
    """match a secondary structure, allowing a certain number of
    disagreements"""


class EnergyType(IntEnum):
    loop = 0
    """[default]: no volume or association terms included"""
    volume = 1
    """include dG_volume (no clear interpretation)"""
    complex = 2
    """include dG_assoc (NUPACK complex microstate energy, sans symmetry terms)"""
    tube = 3
    """include dG_volume + dG_assoc (system state energy when summed over complexes)"""


class TransitionType(IntEnum):
    unimol = 0
    """[default]: unimolecular transition"""
    bimol_join = 1
    """bimolecular join (energies are irrelevant)"""
    bimol_break = 2
    """bimolecular break (energies are relevant)"""


TrajLogComplex = namedtuple(
    "TrajLogComplex", [
        "seed", "id", "strand_names",
        "sequence", "structure", "energy", "enthalpy"])
TrajLogTransition = namedtuple(
    "TrajLogTransition", ["time", "stop_conditions"])
TrajLogResult = namedtuple(
    "TrajLogResult", [
        "seed", "stop_flag", "time",
        "collision_rate", "stop_tag"])
TrajLogResultNoCollision = namedtuple(
    "TrajLogResultNoCollision",
    TrajLogResult._fields[:-2] + (TrajLogResult._fields[-1],))


class Options:
    """
    The main wrapper for controlling a Multistrand simulation. Has an interface
    for returning results.
    """
    dangleToString = ["None", "Some", "All"]
    RateMethodToString = ["None", "Metropolis", "Kawasaki", "Arrhenius"]

    # Parameter type. Vienna is depreciated.
    viennaModel = 0
    nupackModel = 1
    parameterTypeToString = ["Vienna", "Nupack" ]
    substrateToString = ["Invalid", "RNA", "DNA"]
    arrheniusParams = [
        f"{feature}{eqcls}" for feature in ["lnA", "E"] for eqcls in
        "Stack Loop End StackLoop StackEnd LoopEnd StackStack".split()]

    # translation
    simulationMode = {"Normal"    :         Literals.first_passage_time,
                      "First Step":         Literals.first_step,
                      "Transition":         Literals.transition,
                      "Trajectory":         Literals.trajectory,
                      "First Passage Time": Literals.first_passage_time}

    cotranscriptional_rate_default = 0.001  # 1 nt added every 1 ms

    activestatespace = False

    def __init__(self, *args, **kargs):
        """
        Keyword Arguments:
        ------------------
        start_state: List[Complex]    -- Initial state of the system.
        substrate_type: str (DNA|RNA) -- Material used in the energy model.
        dangles: str (None|Some|All)  -- Dangle terms used in the energy model.
        uniscale: float               -- Scaling of the unimolecular rate coefficient.
        biscale: float                -- Scaling of the bimolecular rate coefficient.
        simulation_time: float        -- Cap on the simulation time per trajectory.
        step_count: int               -- Cap on the Markov jumps per trajectory.
        num_simulations: int          -- Number of trajectories to run.
        rate_method: str (Kawasaki|Metropolis|Arrhenius)

        Optional Arguments:
        -------------------
        initial_seed: int             -- PRNG seed at the start of the first trajectory.
        state_seed: (int,int,int)     -- full PRNG buffer to restart a trajectory mid-way.
        simulation_start_time: float  -- Initial time point (in combination with `state_seed`).

        If rate_method == 'Arrhenius', please also set:
            lnAEnd, lnALoop, lnAStack, lnAStackStack, lnALoopEnd, lnAStackEnd,
            lnAStackLoop, EEnd, ELoop, EStack, EStackStack, ELoopEnd, EStackEnd, EStackLoop,
            bimolecular_rate.
        """

        ##################################################
        #                                                #
        # Data Members                                   #
        # ->Members new to the python implementation     #
        #                                                #
        ##################################################

        """ Pipe to let Multistrand know the version from ../__init__.py """
        self.ms_version = float(__version__)

        self.errorlog = []
        """ Keeps lines relating to possible errors or warnings that
        should be reported to the user. Usually issues relating to the
        input file or parameters with odd values.

        TODO: implement some functions to report the errors found here.
        """
        self.full_trajectory = []
        self.full_trajectory_times = []
        self.full_trajectory_arrType = []
        self.trajectory_complexes = []
        self.trajectory_state_count = 0
        self._current_end_state = []
        self._current_transition_list = []
        self.trajectory_current_time = 0.0
        self.current_graph = None

        self._reusable_sim_system = None
        """
        Reference to the first `SimSystem` instance created from `self`,
        providing a reusable `EnergyModel` (C++) for subsequent instances.
        """

        self.verbosity = 1
        """ Indicates how much output will be generated for each trajectory run.
        Value = 0:  No end state reported, no warnings for timeout and nonitial steps
        Value = 1:  No end states reports, warnings active
        Value = 2:  Warnings and end states reports to stdout
        Value = 3:  Print debugging information from SimulationSystem to stdout
        """

        self.print_initial_first_step = False
        """
        If True, this value will print the initial state in First Step Mode to
        the trajectory with a timestamp of -1.0
        """

        self.cotranscriptional = False
        """
        If True, enables the cotranscriptional simulation mode. The mode works
        only when a single strand is supplied. Starting with the initial 8
        nucleotides, the simulation adds a nucleotide on the 3' end every 1
        millisecond.
        """

        self.cotranscriptional_rate = self.cotranscriptional_rate_default
        """
        By default, the cotranscriptional mode adds one nucleotide every 1
        millisecond.
        """

        #############################################
        #                                           #
        # Data Members: Energy Model                #
        #                                           #
        #############################################
        # See accessors below
        self._start_state: List[Complex] = []

        self.gt_enable = True
        """
        If `False`, GT base pairs are penalized by 10000 kcal/mol
        in the energy model for open loops and multi-loops.
        """

        self.log_ml = False
        """ Use logarithm to compute multiloop energy contributions?
        int          False (0): Do not use the logarithm.

        If True, uses log to compute one component of the multiloop energy, for
        loops of length > 6. Otherwise, uses the usual linear approximation.
        Using the log formula is slightly more expensive as it makes computation
        time for a multiloop scale with the number of adjoining helices.
        """

        self._join_concentration: float = 1.0
        """ concentration for V calcs
        Units are in M (molar), and represent the concentration of a single
        unique strand in the system. The volume simulated is then chosen using
        this parameter.
        """

        self._temperature_celsius = 37.0
        self._temperature_kelvin = 310.15

        self.rate_scaling = None
        """FD: This is a legacy option that sets unimolecular and bimolecular
        scaling automatically if set"""

        self._unimolecular_scaling: float = -1.0
        """ Rate scaling factor for unimolecular reactions."""

        self._bimolecular_scaling: float = -1.0
        """ Rate scaling factor for bimolecular reactions."""

        self._rate_method: int = Literals.kawasaki
        """ Choice of methods for determining forward/reverse rates. """

        self._dangles: int = Literals.dangles_some
        """ Dangles options for the energy model.

        None [0]: Do not include any dangles terms in the energy model.
        Some [1]: Some dangles terms.  (Nupack Default)
        All  [2]: Include all dangles terms, including odd overlapping ones.
        """

        self._parameter_type: int = self.nupackModel
        """ Which type of energy model parameter file to use.

        Vienna [0]: No longer well tested. Recommend not using.
        Nupack [1]: Includes some multi-complex parameters, otherwise
                    nearly the same as mfold style files.
        """

        self.substrate_type: int = Literals.substrateDNA
        """ What substrate's parameter files to use.

        Invalid [0]: Indicates we should not auto-search for a param file.
        RNA     [1]: RNA parameters are available from Nupack.
        DNA     [2]: DNA parameters are available from Nupack.

        The Vienna model support is deprecated.
        """

        ####################
        #
        # BEGIN simmode
        #
        ####################

        self._simulation_mode = Literals.first_passage_time
        """ The simulation mode: how we want the simulation system to
        perform the main loop.
        """

        self._simulation_start_time: float = 0.0
        self._simulation_time: float = 600.0
        self._step_count: int = 0
        self._num_simulations: int = 1

        self._initial_seed: Optional[int] = None
        self._state_seed: Optional[Tuple[int, int, int]] = None

        self.name_dict = {}
        """ Dictionary from strand name to a list of unique strand objects
        having that name.

        Type         Default
        dict         {}

        Modified when start state is added. Used as a lookup when stop states
        are added.
        """

        # allocate Arrhenius model parameters
        for par in self.arrheniusParams:
            setattr(self, par, -0.1)

        """ These are undocumented adjustments to the energy model """
        self.dSA = -0.0
        self.dHA = -0.0

        """
        Buffer conditions
        """
        self.sodium = 1.0
        self.magnesium = 0.0

        ####################
        #
        # BEGIN startstop
        #
        ####################

        # See accessors below
        self._stop_conditions = []
        self._use_stop_conditions = False

        self.stop_count = 0
        """ The number of stop states. Equivalent to 'len(self.stop_conditions)'.

        Type         Default
        int          0

        Incremented automatically when a stop state is added. Should not be
        modified externally.
        """

        self.output_time = -1.0
        """ The amount of time (in seconds) to wait between outputs of
        trajectory information.

        Type         Default
        float        -1.0

        A value of -1.0 corresponds to not basing outputs on output_time
        (but perhaps outputting based on some other condition). A value of 0
        means output as often as possible.
        """

        self._output_interval: int = -1
        """ The number of states between outputs of trajectory information.

        Type         Default
        int          -1

        A value of -1 corresponds to not basing outputs on output_interval
        (but perhaps outputting based on some other condition). A value of 1
        means output every state, 2 means every other state, and so on.
        """

        self.current_interval = 0
        """ Current value of output state counter.

        Type         Default
        int          0

        When current_interval is equal to output_interval, the output state is
        True, and otherwise the output state is False. This is modified by
        increment_output_state, and probably shouldn't be used externally."""

        self.output_state = False
        """ Indicates whether output should be reported.

        Type         Default
        boolean      False

        Value should be True if self.current_interval == self.output_interval
        and False otherwise.
        """

        self.interface = Interface()

        ##############################
        #
        # End of __init__: call the keyword hook fn.
        #
        ##############################

        self.__init_keyword_args(self, *args, **kargs)

    def clear(self) -> None:
        """
        Erase all dynamic simulation state, *excluding* the `SimSystem` that is
        stored in order to reuse its `EnergyModel`.

        Note that `Options.initial_seed` and `Options.state_seed` are considered
        static configuration to be kept.
        """
        self.full_trajectory = []
        self.full_trajectory_times = []
        self.full_trajectory_arrType = []
        self.trajectory_complexes = []
        self.trajectory_state_count = 0
        self.trajectory_current_time = 0.0
        self._current_end_state = []
        self._current_transition_list = []
        self.current_graph = None
        self.current_interval = 0
        self.output_state = False
        self.errorlog = []

        self.interface = Interface()

    def free_sim_system(self) -> None:
        """
        Release the `SimSystem` object (C++) that is stored in order to reuse
        its `EnergyModel`. Even without calling this method, the `SimSystem`
        will be freed when `self` is garbage-collected. Hence, this method is
        only required if an `Options` object needs to be pickled, e.g., when
        aggregating results from parallel simulations.
        """
        self._reusable_sim_system = None

    def restart_from_checkpoint(self, idx: int) -> "Options":
        """
        Create a new `Options` instance to restart a trajectory mid-way, using
        the checkpoint `self.full_trajectory[idx]`.
        """
        # retrieve necessary data from checkpoint
        cp_log = self.full_trajectory[idx]
        cp_state = self.checkpoint_state(cp_log)
        cp_time = self.full_trajectory_times[idx]
        cp_state_seed = cp_log[0].seed
        assert all(cmplx.seed == cp_state_seed for cmplx in cp_log)

        # construct a new `Options` instance with modified initial conditions
        self.free_sim_system()
        opt = copy.deepcopy(self)
        opt.clear()
        assert opt == self
        opt._start_state = []
        opt.start_state = cp_state
        opt.simulation_start_time = cp_time
        opt.state_seed = cp_state_seed
        return opt

    def strand_names(self):
        names = {f"{s.id}:{s.name}": s
                 for s in chain(*(c.strand_list for c in self.start_state))}
        assert len(names) == len(np.unique([s.id for s in names.values()]))
        return names

    def checkpoint_state(self, checkpoint: List[TrajLogComplex]):
        """
        Reconstruct a new `Options.start_state` representation from a
        checkpoint in the previous trajectory log. Used by
        `Options.restart_from_checkpoint()`.
        """
        assert all(isinstance(cmplx, TrajLogComplex) for cmplx in checkpoint)
        strand_names = self.strand_names()
        return [
            Complex(strands=[strand_names[n]
                             for n in cmplx.strand_names.split(',')],
                    structure=cmplx.structure)
            # the input complex list is reversed on the way to `SComplexList`
            for cmplx in checkpoint[::-1]]

    def __eq__(self, other: "Options") -> bool:
        """
        Compare configurations syntactically, ignoring random seeds and
        simulator state.
        """
        return (
            self.ms_version,
            self.verbosity, self.print_initial_first_step, self.activestatespace,
            self.substrate_type, self.parameter_type, self.parameter_file,
            self.gt_enable, self.log_ml, self.dangles,
            self.cotranscriptional, self.cotranscriptional_rate,
            self.join_concentration, self.temperature,
            self.rate_scaling,
            self.rate_method, self.unimolecular_scaling, self.bimolecular_scaling,
            self.simulation_mode, self.num_simulations,
            self.simulation_start_time, self.simulation_time, self.step_count,
            self.dSA, self.dHA, self.sodium, self.magnesium,
            self.start_state, self.stop_conditions,
            self.output_time, self.output_interval, self.output_state,
        ) == (
            other.ms_version,
            other.verbosity, other.print_initial_first_step, other.activestatespace,
            other.substrate_type, other.parameter_type, other.parameter_file,
            other.gt_enable, other.log_ml, other.dangles,
            other.cotranscriptional, other.cotranscriptional_rate,
            other.join_concentration, other.temperature,
            other.rate_scaling,
            other.rate_method, other.unimolecular_scaling, other.bimolecular_scaling,
            other.simulation_mode, other.num_simulations,
            other.simulation_start_time, other.simulation_time, other.step_count,
            other.dSA, other.dHA, other.sodium, other.magnesium,
            other.start_state, other.stop_conditions,
            other.output_time, other.output_interval, other.output_state,
        ) and (
            self.rate_method != Literals.arrhenius or
            all(getattr(self, par) == getattr(other, par)
                for par in self.arrheniusParams)
        )

    def legacyRates(self):
        warningmsg = "Warning! rate_scaling is set, enabling support for legacy code. Now setting rate defaults for "

        if self.temperature == 298.15 and self.rate_method == Literals.kawasaki:
            warningmsg += "Kawasaki 25 C"
            self.JSKawasaki25()
        elif self.temperature == 310.15 and self.rate_method == Literals.kawasaki:
            warningmsg += "Kawasaki 37 C"
            self.JSKawasaki37()
        elif self.temperature == 298.15 and self.rate_method == Literals.metropolis:
            warningmsg += "Metropolis 25 C"
            self.JSMetropolis25()
        elif self.temperature == 310.15 and self.rate_method == Literals.metropolis:
            warningmsg += "Metropolis 37 C"
            self.JSMetropolis37()
        else:
            warningmsg += "JS-Default"
            self.JSDefault()

        print(warningmsg)
        self.rate_scaling = None

    def JSDefault(self):
        """
        Default rates (Kawasaki at 37 degree Celcius) from Joseph Schaeffer's
        thesis.
        """
        self.JSKawasaki37()

    def JSMetropolis25(self):
        """
        Default rates for Metropolis at 25 degree Celcius, from Joseph
        Schaeffer's thesis.
        """
        self.substrate_type = Literals.substrateDNA
        self.rate_method = Literals.metropolis

        self.unimolecular_scaling = 4.4e8
        self.bimolecular_scaling = 1.26e6

    def JSKawasaki25(self):
        """
        Default rates for Kawasaki at 25 degree Celcius, from Joseph Schaeffer's
        thesis.
        """
        self.substrate_type = Literals.substrateDNA
        self.rate_method = Literals.kawasaki

        self.unimolecular_scaling = 6.1e7
        self.bimolecular_scaling = 1.29e6

    def JSKawasaki37(self):
        """
        Default rates for Kawasaki at 37 degree Celcius, from Joseph Schaeffer's
        thesis.
        """
        self.substrate_type = Literals.substrateDNA
        self.rate_method = Literals.kawasaki

        self.unimolecular_scaling = 1.5e8
        self.bimolecular_scaling = 1.38e6

    def JSMetropolis37(self):
        """
        Default rates for Metropolis at 37 degree Celcius, from Joseph
        Schaeffer's thesis.
        """
        self.substrate_type = Literals.substrateDNA
        self.rate_method = Literals.metropolis

        self.unimolecular_scaling = 7.3e8
        self.bimolecular_scaling = 1.40e6

    def DNA23Metropolis(self):
        """
        Parameters for Metropolis at 25 degree Celcius, from the DNA23
        conference (55th walker).

        Reference:
        ----------
        Sedigheh Zolaktaf, Frits Dannenberg, Xander Rudelis, Anne Condon,
        Joseph M. Schaeffer, Mark Schmidt, Chris Thachuk, and Erik Winfree.
        2017. ‘Inferring Parameters for an Elementary Step Model of DNA
        Structure Kinetics with Locally Context-Dependent Arrhenius Rates’. In
        DNA Computing and Molecular Programming, edited by Robert Brijder and
        Lulu Qian, 172–87. Lecture Notes in Computer Science. Cham: Springer
        International Publishing. https://doi.org/10.1007/978-3-319-66799-7_12.
        """
        self.substrate_type = Literals.substrateDNA
        self.rate_method = Literals.metropolis

        self.unimolecular_scaling = 2.41686715e+6
        self.bimolecular_scaling = 8.01171383e+5

    def DNA23Arrhenius(self):
        """
        Arrhenius model parameters presented at the DNA23 conference.

        Reference:
        ----------
        Sedigheh Zolaktaf, Frits Dannenberg, Xander Rudelis, Anne Condon,
        Joseph M. Schaeffer, Mark Schmidt, Chris Thachuk, and Erik Winfree.
        2017. ‘Inferring Parameters for an Elementary Step Model of DNA
        Structure Kinetics with Locally Context-Dependent Arrhenius Rates’. In
        DNA Computing and Molecular Programming, edited by Robert Brijder and
        Lulu Qian, 172–87. Lecture Notes in Computer Science. Cham: Springer
        International Publishing. https://doi.org/10.1007/978-3-319-66799-7_12.
        """
        self.substrate_type = Literals.substrateDNA
        self.rate_method = Literals.arrhenius

        self.unimolecular_scaling = 1.0
        self.bimolecular_scaling = 1.60062641e-2

        self.lnAStack = 1.41839430e+1
        self.EStack = 5.28692038e+0

        self.lnALoop = 1.64236969e+1
        self.ELoop = 4.46143369e+0

        self.lnAEnd = 1.29648159e+1
        self.EEnd = 3.49798154e+0

        self.lnAStackLoop = 5.81061725e+0
        self.EStackLoop = -1.12763854e+0

        self.lnAStackEnd = 1.75235569e+1
        self.EStackEnd = 2.65589869e+0

        self.lnALoopEnd = 2.42237267e+0
        self.ELoopEnd = 8.49339120e-2

        self.lnAStackStack = 8.04573830e+0
        self.EStackStack = -6.27121400e-1

    def DNA29Arrhenius(self, sample_idx: Optional[int] = None):
        """
        Arrhenius model parameters presented at the DNA29 conference.

        Arguments:
        ----------
        sample_idx : Optional[int]
          The default parameter vector is the one used in the case study
          (Sec. 5). If an index is provided, then it will be retrieved from a
          file storing all MCMC samples of the approximate Bayesian posterior
          over parameter vectors (Sec. 4.3, target: PE, inference method: RWM).

        Reference:
        ----------
        Jordan Lovrod, Boyan Beronov, Chenwei Zhang, Erik Winfree, and Anne
        Condon. 2023. “Revisiting Hybridization Kinetics with Improved
        Elementary Step Simulation.” In 29th International Conference on DNA
        Computing and Molecular Programming (DNA 29), edited by Ho-Lin Chen and
        Constantine G. Evans, 276:5:1-5:24. Leibniz International Proceedings in
        Informatics (LIPIcs). Dagstuhl, Germany: Schloss Dagstuhl –
        Leibniz-Zentrum für Informatik. https://doi.org/10.4230/LIPIcs.DNA.29.5.
        """
        self.substrate_type = Literals.substrateDNA
        self.rate_method = Literals.arrhenius
        self.unimolecular_scaling = 1.0

        if sample_idx is None:
            self.bimolecular_scaling = 5.40306772408701e-2

            self.lnAStack = 7.317929742353791e+0
            self.EStack = 1.371171987160233e+0

            self.lnALoop = 9.398865994892086e+0
            self.ELoop = 6.675295990888666e-1

            self.lnAEnd = 1.4182060613367684e+1
            self.EEnd = 3.428091116849789e+0

            self.lnAStackLoop = 1.124903532063121e+1
            self.EStackLoop = 1.1475124044307237e+0

            self.lnAStackEnd = 1.2975095181504653e+1
            self.EStackEnd = -1.2455835738995455e+0

            self.lnALoopEnd = 4.453062436943478e-1
            self.ELoopEnd = -2.0459543870895307e+0

            self.lnAStackStack = 1.15977270694611e+1
            self.EStackStack = -2.4564902855673165e+0
        else:
            assert isinstance(sample_idx, int)
            from pandas import read_csv
            df = read_csv(importlib.resources.files(
                "multistrand._options.parameters") / "dna29_arrhenius.csv")
            params = df.iloc[sample_idx]
            for p in (["bimolecular_scaling"] + self.arrheniusParams):
                setattr(self, p, params[p])

    # FD: After temperature, substrate (RNA/DNA) or danlges is updated, we
    # attempt to update boltzmann samples.
    def updateBoltzmannSamples(self):
        for c in self._start_state:
            c.set_boltzmann_parameters(
                self.dangleToString[self.dangles],
                self.substrateToString[self.substrate_type],
                self._temperature_celsius, self.sodium, self.magnesium)
            self.warn_Boltzmann_sample_wo_GT(c)

    def warn_Boltzmann_sample_wo_GT(self, c: Complex):
        if c.boltzmann_sample and not self.gt_enable:
            raise Warning(
                "Attempting to use Boltzmann sampling, but GT pairing is "
                "disabled. Energy model of Multistrand will not match that of "
                "the NUPACK sampling method.")

    @property
    def simulation_start_time(self):
        """
        Initial time point for each trajectory (default: 0.0). Used by
        `Options.restart_from_checkpoint()`.
        """
        return self._simulation_start_time

    @simulation_start_time.setter
    def simulation_start_time(self, value):
        self._simulation_start_time = float(value)

    @property
    def simulation_time(self):
        """
        Maximum time (in seconds) allowed for each trajectory.

        Type         Default
        double       600.0
        """
        return self._simulation_time

    @simulation_time.setter
    def simulation_time(self, value):
        self._simulation_time = float(value)

    @property
    def step_count(self):
        """
        Maximum number of Markov jumps allowed for each trajectory.

        Type         Default
        int          0 (unlimited)
        """
        return self._step_count

    @step_count.setter
    def step_count(self, value):
        self._step_count = int(value)

    @property
    def num_simulations(self):
        """ Total number of trajectories to run. """
        return self._num_simulations

    @num_simulations.setter
    def num_simulations(self, value):
        self._num_simulations = int(value)

    @property
    def output_interval(self):
        return self._output_interval

    @output_interval.setter
    def output_interval(self, value):
        self._output_interval = int(value)

    @property
    def bimolecular_scaling(self):
        if self.rate_scaling != None :
            self.legacyRates()
        return self._bimolecular_scaling

    @bimolecular_scaling.setter
    def bimolecular_scaling(self, value):
        self._bimolecular_scaling = float(value)

    @property
    def unimolecular_scaling(self):
        if self.rate_scaling != None:
            self.legacyRates()
        return self._unimolecular_scaling

    @unimolecular_scaling.setter
    def unimolecular_scaling(self, value):
        self._unimolecular_scaling = float(value)

    @property
    def join_concentration(self):
        return self._join_concentration

    @join_concentration.setter
    def join_concentration(self, value):
        self._join_concentration = float(value)

    # FD: Shadow variables for danlges because we need to observe changes
    # (and update boltzmann samples accordingly)
    @property
    def dangles(self):
        return self._dangles

    @dangles.setter
    def dangles(self, value):
        if isinstance(value, str):
            value = self.dangleToString.index(value)
        self._dangles = int(value)
        assert self.dangles in range(3)
        self.updateBoltzmannSamples()

    # FD: Shadow parameter so that boltzmann samples can be updated when this
    # parameter is set. In a better control flow, complexes themselves might fetch the right
    # constants just before evaluating their boltzmann samples.
    @property
    def substrate_type(self):
        return self._substrate_type

    @substrate_type.setter
    def substrate_type(self, value):
        if isinstance(value, str):
            value = self.substrateToString.index(value)
        self._substrate_type = int(value)
        assert self._substrate_type in range(1, 3)
        self.updateBoltzmannSamples()

    @property
    def parameter_type(self):
        return self._parameter_type

    @parameter_type.setter
    def parameter_type(self, value):
        if isinstance(value, str):
            value = self.parameterTypeToString.index(value)
        self._parameter_type = int(value)
        assert self.parameter_type in range(2)

    @property
    def parameter_file(self):
        pfile = ""
        if self.parameter_type == Options.nupackModel:
            substrate = ("rna06" if self.substrate_type == Literals.substrateRNA
                         else "dna04")
            pfile = str(importlib.resources.files("nupack") / "parameters" /
                        f"{substrate}-nupack3.json")
        else:
            raise NotImplementedError(
                "Only NUPACK energy models are currently supported.")
        if not os.path.exists(pfile):
            raise EnvironmentError("Cannot import NUPACK dependency!")
        return pfile

    @property
    def simulation_mode(self):
        return self._simulation_mode

    @simulation_mode.setter
    def simulation_mode(self, value):
        if isinstance(value, str):
            value = self.simulationMode[value]
        self._simulation_mode = int(value)
        assert self.simulation_mode in [16, 48, 256, 128]

    @property
    def rate_method(self):
        return self._rate_method

    @rate_method.setter
    def rate_method(self, value):
        if isinstance(value, str):
            value = self.RateMethodToString.index(value)
        self._rate_method = int(value)
        assert self.rate_method in range(1, 4)

    # FD: Following same listener pattern for sodium, magnesium, so that changes
    # are propagated to complexes.
    @property
    def sodium(self):
        return self._sodium

    @sodium.setter
    def sodium(self, value):
        self._sodium = float(value)
        self.updateBoltzmannSamples()

    @property
    def magnesium(self):
        return self._magnesium

    @magnesium.setter
    def magnesium(self, value):
        self._magnesium = float(value)
        self.updateBoltzmannSamples()

    @property
    def boltzmann_sample(self):
        raise ValueError('Options.boltzmann_sample is now depreciated. Use Complex.boltzmann_sample instead.')

    @boltzmann_sample.setter
    def boltzmann_sample(self, val):
        raise ValueError('Options.boltzmann_sample is now depreciated. Use Complex.boltzmann_sample instead.')

    @property
    def start_state(self):
        """ Get the start state, i.e. a list of Complex objects.

        Type         Default
        list         []

        This should be used by ssystem.cc to get the (potentially sampled)
        start state.
        """
        return self._start_state

    @start_state.setter
    def start_state(self, *args):
        """ Set the start state, i.e. a list of Complex objects.

        Type         Default
        list         []

        The start state should be set (e.g. by the parser) so trajectories know
        how to start.
        """
        # Error checking first
        if self._start_state != []:
            raise Exception("Start state should only be set once.")
        if len(args) == 0 or len(args[0]) == 0:
            raise ValueError("No start state given.")

        # deduce our input from the type of args[0].
        # Copy the input list because it's easy to do and it's safer

        if isinstance(args[0], Complex):
            # args is a list of complexes
            vals = copy.deepcopy(args)
        elif len(args) == 1 and hasattr(args[0], "__iter__"):
            vals = copy.deepcopy(args[0])
        else:
            raise ValueError("Could not comprehend the start state you gave me.")

        # vals is now an iterable over our starting configuration
        for i in vals:
            if isinstance(i, Complex):
                self._add_start_complex(i)
            else:
                raise TypeError(f"Start states must be Complexes. "
                                f"Received something of type {type(i)}.")

    def _add_start_complex(self, c: Complex):
        self._start_state.append(c)
        c.set_boltzmann_parameters(
            self.dangleToString[self.dangles],
            self.substrateToString[self.substrate_type],
            self._temperature_celsius, self._sodium, self._magnesium)
        self.warn_Boltzmann_sample_wo_GT(c)

    @property
    def initial_seed(self):
        """
        Configuration of the initial PRNG seed (32 bits) to use at the start
        of the current trajectory. If `None` (default), a random seed will be
        chosen. See also: `Options.state_seed`.

        C type: long
        """
        return self._initial_seed if self.initial_seed_flag else None

    @initial_seed.setter
    def initial_seed(self, seed):
        self._initial_seed = int(seed)
        assert 0 <= abs(self._initial_seed) <= 1 << 31

    @property
    def initial_seed_flag(self):
        return self._initial_seed != None

    @property
    def state_seed(self):
        """
        Full PRNG buffer (48 bits) to set at the start of the current
        trajectory. Used by `Options.restart_from_checkpoint()`.
        If `None` (default), then fall back to `Options.initial_seed`.

        C type: unsigend short[3]
        """
        return self._state_seed if self.state_seed_flag else None

    @state_seed.setter
    def state_seed(self, seed):
        self._state_seed = tuple(map(int, seed))
        assert len(self._state_seed) == 3
        assert all(0 <= s < 1 << 16 for s in self._state_seed)

    @property
    def state_seed_flag(self):
        return self._state_seed != None

    @property
    def stop_conditions(self):
        """ The stop states, i.e. a list of StopCondition objects.

        Type         Default
        list         []

        Stop states should be added to this list (e.g. by the parser) so
        trajectories know when to end.
        """
        return self._stop_conditions

    @stop_conditions.setter
    def stop_conditions(self, stop_list):
        """ The stop states, i.e. a list of StopCondition objects.

        Type         Default
        list         []

        Stop states should be added to this list (e.g. by the parser) so
        trajectories know when to end.
        """
        # Error checking
        if self._stop_conditions != []:
            raise Exception("Stop conditions should be set only once.")

        # Type checking
        for item in stop_list:
            if not isinstance(item, StopCondition):
                raise TypeError(f"All items must be 'StopCondition', not '{type(item)}'.")

        # Copy the input list because it's easy to do and it's safer
        stop_list = copy.deepcopy(stop_list)

        for scond in stop_list:
            scond.verify_no_boltzmann()

        # Set the internal data members
        self.stop_count = len(stop_list)
        self._stop_conditions = stop_list
        self._use_stop_conditions = True

    @property
    def use_stop_conditions(self):
        """ Indicates whether trajectories should end when stop states
        are reached.

        Type            Default
        boolean         False: End trajectory upon reaching max time only.

        Defaults to ending trajectories only on the max simulation
        time, but setting any stop conditions automatically changes
        this to True, and it will stop whenever it reaches a stop
        condition, or at max time [whichever comes first].
        """
        return self._use_stop_conditions

    @use_stop_conditions.setter
    def use_stop_conditions(self, val):
        if val == True and len(self._stop_conditions) == 0:
             raise Warning("Options.use_stop_conditions was set to True, but no stop conditions have been defined!")
        self._use_stop_conditions = val

    @property
    def increment_output_state(self):
        """ Modifies self.current_interval and self.output_state as
        necessary based on self.output_interval.
        """
        if self.output_interval == None or self.output_interval < 0:
            raise ValueError("output_interval has invalid value: %s" % self.output_interval)

        elif self.current_interval > self.output_interval:
            raise ValueError("current_interval has invalid value: %s" % self.current_interval)

        elif self.current_interval == self.output_interval:
            self.current_interval == 0

        else:
            self.current_interval += 1

        self.output_state = (self.current_interval == self.output_interval)
        return None

    @property
    def temperature(self):
        """
        Temperature, in degrees Kelvin.

        Arguments:
        temperature [type=float,default=310.15] -- Standard units of Kelvin.
               Default value corresponds to 37(C).

        This is used mostly in scaling energy model terms, and those
        scaling factors always use Kelvin. Note that when set, the
        Options object will try to make sure it's actually a 'sane'
        value, as follows:

        Temperatures in the range [0,100] are assumed to be Celsius,
        and are converted to Kelvin.

        Temperatures in the range [273,373] are assumed to be in
        Kelvin.

        Any temperature outside these ranges is set as the temperature
        in Kelvin, and a warning is raised. If any conversion takes
        place, a message is added to the Options object's errorlog.
        """
        return self._temperature_kelvin

    @temperature.setter
    def temperature(self, val):
        """ performs some sanity checking and provides a log error if
        perhaps the user was confused.

        Assumptions: Input should be in Kelvin, if it's not in a
        'reasonable' range for Kelvin, convert to Celsius if it's in a
        reasonable range for C [sending an output to the error log
        warning that it did so], otherwise error.

        Reasonable ranges:
            [0,100]  : Celsius
            [273,373]: Kelvin
            Others:    If you want a Fahrenheit reasonable range, I think you
                       might be unreasonable. Also, it overlaps with Celsius a bit much.

        Yes, these ranges are quite generous.
        """
        if C2K < val < C2K + 100:
            self._temperature_kelvin = val
            self._temperature_celsius = val - C2K
            self.updateBoltzmannSamples()

        elif 0.0 < val < 100.0:
            self._temperature_celsius = val
            self._temperature_kelvin = val + C2K
            self.updateBoltzmannSamples()
            self.errorlog.append("Warning: Temperature was set at the value [{0}]. We expected a value in Kelvin, or with appropriate units.\n         Temperature was automatically converted to [{1}] degrees Kelvin.\n".format(val, self._temperature_kelvin))

        else:
            self._temperature_kelvin = val
            self._temperature_celsius = val - C2K
            self.updateBoltzmannSamples()
            self.errorlog.append("Warning: Temperature was set at the value [{0}]. This is outside the normal range of temperatures we expect, so it was assumed to be in Kelvin.\n".format(val))
            raise Warning("Temperature did not fall in the usual expected ranges. Temperatures should be in units Kelvin, though the range [0,100] is assumed to mean units of Celsius.")

    def make_unique(self, strand):
        """Returns a new Strand object with a unique identifier replacing the
        old id. Also adds the new strand to self.name_dict[strand.name].
        """
        new_strand = Strand(
            self.unique_id, strand.name, strand.sequence, strand.domain_list)
        self.unique_id += 1

        try:
            self.name_dict[strand.name].append(new_strand)
        except KeyError:
            self.name_dict[strand.name] = [new_strand]
        return new_strand

    @property
    def add_result_status_line(self):
        return None

    @add_result_status_line.setter
    def add_result_status_line(self, val):
        log = TrajLogResultNoCollision._make(val)
        self.interface.add_result(log, res_type='status_line')
        if len(self._current_end_state) > 0:
            self.interface.end_states.append(self._current_end_state)
            self._current_end_state = []

    @property
    def add_result_status_line_firststep(self, val):
        return None

    @add_result_status_line_firststep.setter
    def add_result_status_line_firststep(self, val):
        log = TrajLogResult._make(val)
        self.interface.add_result(log, res_type='firststep')
        if len(self._current_end_state) > 0:
            self.interface.end_states.append(self._current_end_state)
            self._current_end_state = []

    @property
    def add_complex_state_line(self):
        return None

    @add_complex_state_line.setter
    def add_complex_state_line(self, val):
        log = TrajLogComplex._make(val)
        self._current_end_state.append(log)
        if self.verbosity > 1:
            print(f"{log.seed}: [{log.id}] "
                  f"'{log.strand_names}': {log.energy}"
                  f"\n{log.sequence}\n{log.structure}\n")

    @property
    def add_transition_info(self):
        return None

    @add_transition_info.setter
    def add_transition_info(self, val):
        log = TrajLogTransition._make(val)
        self._current_transition_list.append(log)

    @property
    def add_trajectory_complex(self):
        return None

    @add_trajectory_complex.setter
    def add_trajectory_complex(self, val):
        log = TrajLogComplex._make(val)
        self.trajectory_complexes.append(log)

    @property
    def add_trajectory_current_time(self):
        return None

    @add_trajectory_current_time.setter
    def add_trajectory_current_time(self, val):
        self.trajectory_current_time = val
        self.trajectory_state_count += 1
        self.full_trajectory.append(self.trajectory_complexes)
        self.full_trajectory_times.append(self.trajectory_current_time)
        self.trajectory_complexes = []

    @property
    def add_trajectory_arrType(self):
        return None

    @add_trajectory_arrType.setter
    def add_trajectory_arrType(self, val):
        self.full_trajectory_arrType.append(val)

    @property
    def interface_trajectory_seed(self) -> Optional[int]:
        """
        This is the C++ runtime value of the initial PRNG seed for the current
        trajectory, as reported by the simulator after initialisation. For
        configuring this value, use `Options.initial_seed`.
        """
        return self.interface.trajectory_seed

    @interface_trajectory_seed.setter
    def interface_trajectory_seed(self, val):
        self.interface.trajectory_seed = int(val)
        get_structure = lambda s: (s._last_boltzmann_structure
                                   if s.boltzmann_sample else s._fixed_structure)
        self.interface.start_structures[val] = list(
            map(get_structure, self._start_state))

    @property
    def increment_trajectory_count(self):
        self.interface.increment_trajectory_count()
        if len(self._current_transition_list) > 0:
            self.interface.transition_lists.append(self._current_transition_list)
            self._current_transition_list = []

    def __init_keyword_args(self, *args, **kargs):
        """ Create an options object [with default (presumably useful) values]

        Now with new and improved argument lists!

        Any default Options attribute can be set just by passing a
        keyword argument with the desired value, e.g.
        Options(simulation_time=100.0)

        Listed below are some shortcuts, for attributes that have a
        range of options, or shortened names for long attributes.


        Keyword Argument  |  Options
        dangles           |  'None', 'Some', 'All'
        parameter_type    |  'Nupack', 'Vienna'
        substrate_type    |  'DNA','RNA'

        sim_time          | [simulation_time] Max time to simulate
        num_sims          | [num_simulations] Number of trajectories to run
        biscale           | [bimolecular_scaling] Bimolecular scaling constant
        uniscale          | [unimolecular_scaling] Unimolecular scaling constant

        start_state       |  List of Complexes

        ...
        More to come!"""
        arg_lookup_table = {
            'biscale': lambda x: self.__setattr__('bimolecular_scaling', x),
            'uniscale': lambda x: self.__setattr__('unimolecular_scaling', x),
            'num_sims': lambda x: self.__setattr__('num_simulations', x),
            'sim_time': lambda x: self.__setattr__('simulation_time', x),
            'concentration': lambda x: self.__setattr__('join_concentration', x)
            }

        # FD: Start throwing errors if not in the right format
        for key, value in kargs.items():
            if key == "sim_time":
                if not isinstance(value, (float)):
                    raise Warning("Please provide sim_time as float")
            if key == "num_sims":
                if not isinstance(value, (int)):
                    raise Warning("Please provide num_sims as int")
            if key == "biscale":
                if not isinstance(value, (float)):
                    raise Warning("Please provide biscale as float")
            if key == "uniscale":
                if not isinstance(value, (float)):
                    raise Warning("Please provide uniscale as float")
            if key == "concentration":
                if not isinstance(value, (float)):
                    raise Warning("Please provide concentration as float")

        for key, value in kargs.items():
            if key in arg_lookup_table:
                arg_lookup_table[key](value)
            # FD: Do some additional parsing for legacy support
            else:
                self.__setattr__(key, value)
