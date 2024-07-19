# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2024 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

from typing import Iterator, Tuple
from functools import partial

import pytest
import numpy as np

from multistrand.options import Options, TransitionType
from multistrand.system import calculate_rate


@pytest.mark.parametrize("energies", [np.random.rand(100, 2)])
@pytest.mark.parametrize(
    "rate_model",
    ["JSMetropolis25", "JSMetropolis37", "DNA23Metropolis",
     "DNA23Arrhenius", "DNA29Arrhenius"])
class Test_Kinetic_Model:

    @staticmethod
    def iterate_energies(energies: np.ndarray) -> Iterator[Tuple[float, float]]:
        for a, b in energies:
            yield min(a, b), max(a, b)

    @staticmethod
    def bidir_rates(opt: Options, transition_type: TransitionType,
                    energies: Tuple[float, float]):
        lo, hi = energies
        return [calculate_rate(a, b, opt, transition_type)
                for a, b in [(lo, hi), (hi, lo)]]

    @classmethod
    def test_base_rate(cls, energies: np.ndarray, rate_model: str):
        """
        Basic check of kinetic model interface. Currently ignores transition
        context features.
        """
        opt = Options()
        getattr(opt, rate_model)()

        for up, down in map(partial(cls.bidir_rates, opt,
                                    TransitionType.unimol),
                            cls.iterate_energies(energies)):
            print(np.array([up, down, opt.unimolecular_scaling]))
            assert down == opt.unimolecular_scaling
            assert 0 < up and up < down
        for up, down in map(partial(cls.bidir_rates, opt,
                                    TransitionType.bimol_break),
                            cls.iterate_energies(energies)):
            print(np.array([up, down, opt.bimolecular_scaling]))
            assert 0 < up and up < down
        for up, down in map(partial(cls.bidir_rates, opt,
                                    TransitionType.bimol_join),
                            cls.iterate_energies(energies)):
            print(np.array([up, down, opt.bimolecular_scaling]))
            assert up == down == opt.bimolecular_scaling


class Test_Parameter_Loading:
    """
    For parameter families with multiple samples (i.e., obtained from Bayesian
    posterior inference), check that the hard-coded default sample matches the
    intended one when loading from the parameter file with all samples.
    """
    @staticmethod
    def test_DNA29():
        opt_default, opt_index = Options(), Options()
        opt_default.DNA29Arrhenius()
        opt_index.DNA29Arrhenius(sample_idx=184)
        assert opt_default == opt_index
