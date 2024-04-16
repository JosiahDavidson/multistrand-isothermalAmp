# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2024 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from contextlib import nullcontext
from typing import List, Tuple, Optional

import numpy as np
import pytest

from multistrand.objects import Complex, Domain, Strand
from multistrand.options import Options, Literals
from multistrand.utils.utility import generate_sequence, printTrajectory
from multistrand.system import SimSystem
from multistrand.concurrent import MergeSim


@pytest.mark.parametrize(
    "rate_model",
    ["JSMetropolis25", "JSMetropolis37", "JSKawasaki25", "JSKawasaki37",
     "DNA23Metropolis", "DNA23Arrhenius"])
class Test_Determinism:

    @pytest.mark.parametrize(
        "toehold_seq, bm_design_B, toehold_extra, structure",
        [("GTGGGT", "ACCGCACGTCACTCACCTCG", "TTT",
          "..(.((.....)).).....((((((+))))))((((((((((((((((((((+))))))))))))))))))))")
        ])
    @pytest.mark.parametrize("seed", list(np.random.randint(1 << 31, size=5)))
    def test_trajectories(cls, rate_model: str,
                          toehold_seq: str, bm_design_B: str, toehold_extra: str,
                          structure: str, seed: int) -> None:
        """
        Compare multiple trajectories from identically configured Multistrand
        simulations, first replaying from the start (with sparse output),
        and then replaying from a checkpoint near the end (with full output).
        """
        # time stepping
        intvl_full, intvl_tail = int(1e2), 1
        intvl_idx, num_intvl = -5, 3

        # build initial state
        toehold = Domain(
            name="toehold", sequence=toehold_seq, length=len(toehold_seq))
        branch_migration_B = Domain(
            name="bm_B", sequence=bm_design_B, seq_length=len(bm_design_B))
        substrate_B = toehold + branch_migration_B
        incumbent_B = Strand(name="incumbent", domains=[branch_migration_B.C])
        incoming_B = substrate_B.C
        start_complex = Complex(
            strands=[incoming_B, substrate_B, incumbent_B], structure=structure)

        # simulate full trajectory twice
        o_full = cls.create_config(rate_model, [start_complex], seed)
        o_full.output_interval = intvl_full
        seeds1, structs1, energies1, times1, end1 = cls.simulate(o_full, False)
        seeds2, structs2, energies2, times2, end2 = cls.simulate(o_full, False)

        # compare trajectories (full / full)
        assert (seeds1 == seeds2).all()
        assert (structs1 == structs2).all()
        assert (energies1 == energies2).all()
        assert (times1 == times2).all()
        assert end1 == end2
        assert len(end1) > 0

        # replay trajectory tail twice
        o_tail = o_full.restart_from_checkpoint(intvl_idx)
        del o_full
        o_tail.output_interval = intvl_tail
        seeds3, structs3, energies3, times3, end3 = cls.simulate(o_tail, False)
        seeds4, structs4, energies4, times4, end4 = cls.simulate(o_tail, False)
        del o_tail

        # compare trajectories (tail / tail)
        assert (seeds3 == seeds4).all()
        assert (structs3 == structs4).all()
        assert (energies3 == energies4).all()
        assert (times3 == times4).all()
        assert end3 == end4
        assert len(end3) > 0

        # compare trajectories (full / tail)
        I_full = np.s_[intvl_idx : intvl_idx + num_intvl : 1]
        I_tail = np.s_[0 : num_intvl * intvl_full : intvl_full]
        assert (seeds1[I_full] == seeds3[I_tail]).all()
        assert (structs1[I_full] == structs3[I_tail]).all()
        assert (energies1[I_full] == energies3[I_tail]).all()
        assert np.allclose(times1[I_full], times3[I_tail], atol=1e-6)

        # FIXME:
        #   It appears that there are still undefined behaviours left
        #   in the C++ simulator, and this test should be helpful for finding
        #   them eventually. For now, we only have approximate consistency of
        #   trajectory replays.
        # assert end1 == end3

    @staticmethod
    def create_config(
        rate_model: str, start_state: List[Complex], seed: int) -> Options:
        opt = Options()
        opt.simulation_mode = Literals.trajectory
        opt.temperature = 37.0
        opt.dangles = 1
        opt.start_state = start_state
        opt.initial_seed = seed
        opt.num_simulations = 1
        opt.simulation_time = 5e-4
        getattr(opt, rate_model)()
        return opt

    @classmethod
    def simulate(cls, opt: Options, verbose: bool) -> Tuple[np.ndarray,...]:
        sys = SimSystem(opt)
        sys.start()
        if verbose:
            printTrajectory(opt, show_seed=True)
        return cls.pack_trajectory(opt)

    @classmethod
    def pack_trajectory(cls, opt: Options) -> Tuple:
        seeds = np.array(
            [s[0].seed for s in opt.full_trajectory], dtype=np.uint16)
        structs = np.array(
            [' '.join([c.structure for c in s]) for s in opt.full_trajectory],
            dtype=str)
        energies = np.array(
            [sum(cmplx.energy for cmplx in s) for s in opt.full_trajectory],
            dtype=np.float64)
        times = np.array(opt.full_trajectory_times, dtype=np.float64)
        end = opt.interface.end_states
        return (seeds, structs, energies, times, end)

    @pytest.mark.parametrize(
        "random_seq, raising_random",
        [(True, pytest.raises(ValueError)), (False, nullcontext())])
    @pytest.mark.parametrize(
        "wrong_mode, raising_mode",
        [(True, pytest.raises(TypeError)), (False, nullcontext())])
    @pytest.mark.parametrize("seqlen", [7, 23])
    @pytest.mark.parametrize("seed", list(np.random.randint(1 << 31, size=2)))
    def test_options_factory(
            self, rate_model: str, seqlen: int, seed: int,
            random_seq: bool, raising_random, wrong_mode: bool, raising_mode):
        """
        Exercise precondition checks for the user-provided `OptionsFactory()`,
        whose output should be deterministic.
        """
        def build_start_state(seq: str) -> List[Complex]:
            d_top = Domain(name="orig",
                           sequence=seq or generate_sequence(seqlen))
            s_top = Strand(name="top", domains=[d_top])
            s_bot = s_top.C
            startTop = Complex(strands=[s_top], structure=".")
            startBot = Complex(strands=[s_bot], structure=".")
            return [startTop, startBot]

        def options_factory(num_sims: int, seq: str) -> Options:
            o = self.create_config(rate_model, build_start_state(seq), seed)
            o.simulation_mode = Literals.first_step
            o.num_sims = num_sims
            o.simulation_time = 1e-6
            return o

        num_sims = 4
        seq = generate_sequence(seqlen)

        sim = MergeSim()
        sim.setNumOfThreads(num_sims)
        with raising_random:
            sim.setOptionsFactory2(options_factory, num_sims,
                                   "" if random_seq else seq)
        sim.setOptionsFactory2(options_factory, num_sims, seq)
        if wrong_mode:
            sim.setPassageMode()
        else:
            sim.setFirstStepMode()
        with raising_mode:
            sim.run()


# ==============================================================================


if __name__ == "__main__":
    toehold_seq = "GTGGGT"
    bm_design_B = "ACCGCACGTCACTCACCTCG"
    toehold_extra = "TTT"
    structure = "..(.((.....)).).....((((((+))))))((((((((((((((((((((+))))))))))))))))))))"
    seed = np.random.randint(1 << 31)
    Test_Determinism().test_trajectories(
        "DNA23Arrhenius", toehold_seq, bm_design_B, toehold_extra, structure,
        seed)
