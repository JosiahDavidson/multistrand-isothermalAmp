# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2025 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This demonstrates the pair-type functionality
For a given complex, Pairtype returns a unique representation.
This is important for hashing functions
"""

from multistrand.system import SimSystem
from multistrand.utils.utility import printTrajectory
from multistrand.experiment import standardOptions, hybridization

    
def doSims(strandSeq, numTraj=2):
    o = standardOptions()
    o.num_simulations = numTraj
    o.output_interval = 1
    o.simulation_time = 0.01
    hybridization(o, strandSeq)
    o.DNA23Metropolis()

    s = SimSystem(o)
    s.start()
    printTrajectory(o, show_seed=True)
    return o, s


def main():
    return doSims("GTGAAACGC")


if __name__ == '__main__':
    main()
