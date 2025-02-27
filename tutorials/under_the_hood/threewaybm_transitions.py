# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2025 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)

"""
This example is similar to hairpin_transistion.py, except that multistranded
complexes are handled. Mainly, we think about what the exact and loose
macrostate differences are -- what do they tell us?

Try it like this, e.g.:
  python -i threewaybm_transitions.py
"""

import numpy as np

from multistrand.objects import *
from multistrand.options import Options, Literals
from multistrand.system import SimSystem


def setup_options_threewaybm(toehold_seq="GTGGGT", bm_design="ACCGCACGTCACTCACCTCG"):

    # the structures below are hard-coded for these lengths
    assert len(toehold_seq)==6
    assert len(bm_design)==20

    toehold = Domain(name="toehold", sequence=toehold_seq)
    branch_migration = Domain(name="bm", sequence=bm_design)
    
    substrate = toehold + branch_migration
    incumbent = Strand(name="incumbent", domains=[branch_migration.C])

    incoming = substrate.C

    # start with 6-base toehold fully bound
    start_complex = Complex(strands=[incoming, substrate, incumbent], structure=".(+)(+)")

    initial_structure_dp   = "....................((((((+))))))((((((((((((((((((((+))))))))))))))))))))"
    six_bases_structure_dp = "..............((((((((((((+))))))))))))((((((((((((((+))))))))))))))......"
    six_bases_loose_dp     = "**************((**********+**********))((************+************))******"
    twelve_bases_struc_dp  = "........((((((((((((((((((+))))))))))))))))))((((((((+))))))))............"
    twelve_bases_loose_dp  = "********((*****************+***************))((******+******))************"
    eighteen_structure_dp  = "..((((((((((((((((((((((((+))))))))))))))))))))))))((+)).................."
    eighteen_loose_dp      = "**((**********************+**********************))((+))******************"

    six_bases_complex           = Complex(strands=[incoming,substrate,incumbent], structure=six_bases_structure_dp)
    twelve_bases_complex        = Complex(strands=[incoming,substrate,incumbent], structure=twelve_bases_struc_dp)
    eighteen_bases_complex      = Complex(strands=[incoming,substrate,incumbent], structure=eighteen_structure_dp)
    six_basesloose_complex      = Complex(strands=[incoming,substrate,incumbent], structure=six_bases_loose_dp)
    twelve_basesloose_complex   = Complex(strands=[incoming,substrate,incumbent], structure=twelve_bases_loose_dp)
    eighteen_basesloose_complex = Complex(strands=[incoming,substrate,incumbent], structure=eighteen_loose_dp)

    disassoc_complex            = Complex(strands=[incumbent], structure=".")   # succesful strand displacement
    failed_complex              = Complex(strands=[incoming], structure="..")   # failed strand displacement attempt

    start_sc          = Macrostate("INITIAL", [(start_complex,Literals.count_macrostate,2)])                 # Within distance 2 of the start_complex state.
    six_sc_exact      = Macrostate("SIX_EXACT", [(six_bases_complex,Literals.exact_macrostate,0)])           # the third parameter is ignored; not needed for exact macrostates
    six_sc_loose      = Macrostate("SIX_LOOSE", [(six_basesloose_complex,Literals.loose_macrostate,2)])      # 8 base pairs defined; must have at least 6 to match.
    twelve_sc_exact   = Macrostate("TWELVE_EXACT", [(twelve_bases_complex,Literals.exact_macrostate,0)])
    twelve_sc_loose   = Macrostate("TWELVE_LOOSE", [(twelve_basesloose_complex,Literals.loose_macrostate,2)])
    eighteen_sc_exact = Macrostate("EIGHTEEN_EXACT", [(eighteen_bases_complex,Literals.exact_macrostate,0)])
    eighteen_sc_loose = Macrostate("EIGHTEEN_LOOSE", [(eighteen_basesloose_complex,Literals.loose_macrostate,2)])

    # why bother giving a list of just one macrostate-def tuple? A Macrostate
    # with a list of multiple tuples give the AND (i.e. intersection) of
    # microstates.

    completed_sc      = StopCondition("stop:COMPLETE", [(disassoc_complex,Literals.dissoc_macrostate,0)])  # incumbent strand fell off
    rejected_sc       = StopCondition("stop:REJECTED", [(failed_complex,Literals.dissoc_macrostate,0)])    # incoming strand fell off

    # join_concentration is not defined, because in this simulation we stop before there's any chance for association steps
    o_exact = Options(simulation_mode="Transition", dangles="Some",
                      num_simulations=5, simulation_time=.01, temperature=310.15,
                      start_state=[start_complex], verbosity=0)
    o_exact.stop_conditions = [start_sc, six_sc_exact, twelve_sc_exact,
                               eighteen_sc_exact, completed_sc, rejected_sc]
    o_exact.JSMetropolis37()

    o_loose = Options(simulation_mode="Transition", dangles="Some",
                      num_simulations=5, simulation_time=.01, temperature=310.15,
                      start_state=[start_complex], verbosity=0)
    o_loose.stop_conditions = [start_sc, six_sc_loose, twelve_sc_loose,
                               eighteen_sc_loose, completed_sc, rejected_sc]
    o_loose.JSMetropolis37()

    return o_exact, o_loose


# mol will be a list of True/False for which transition macrostates the system
# has entered so in_state(mol) returns True if the system is in at least one of
# the listed macrostates.
def in_state(mol): return sum(mol) > 0


# mol is a Boolean descriptor of macrostate occupancy, like mol above.
# a short-hand name for this macrostate (based on the order given in stop_conditions) is provided.
def mol_name(mol):
    charindex = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0'
    names = [charindex[j] for j,i in enumerate(mol) if i]
    if names == []:
        names = charindex[26]
    else:
        names = ",".join(names)
    return names


# t0 and t1 are Boolean descriptors of macrostate occupancy, like mol above.
# here, we provide a printable name for the transition between two macrostate occupancy lists.
def trans_name(t0,t1):
    return mol_name(t0) + ' -> ' + mol_name(t1)


def print_transitions(transition_traj):
    for t in transition_traj:
        print(f"{t[0]:12g} : {mol_name(t[1])}")


# for each simulation, the transition trajectory reports the tuple (time_entered, which_macrostates_the_system_is_now_in)
def parse_transition_lists(transition_traj_list):
    transition_dict = {}

    # the mol1 --> mol2 transition times represent (time of entry into mol1) to (time of entry into mol2)
    for transition_traj in transition_traj_list:
        truncated = [i for i in transition_traj if in_state(i[1])]
        tt = truncated # only keep the entry times and mol states for non-trivial mols

        for i in range(len(tt)-1):
            nm = trans_name(tt[i][1],tt[i+1][1])
            if nm in transition_dict:
                transition_dict[nm].append(tt[i+1][0] - tt[i][0])
            else:
                transition_dict[nm] = [tt[i+1][0] - tt[i][0]]

    return transition_dict


def parse_transition_list(transition_traj_list):
    return parse_transition_lists([transition_traj_list])

    
def print_transition_dict(transition_dict, options=None):
    k = list(transition_dict.keys())
    k.sort() 
    for i in k:
        transition_times = np.array( transition_dict[i] )
        print(f"{i}: {np.mean(transition_times):.2e} ({len(transition_dict[i])})")

    # also print the true names of the macrostates, if an Options object is provided
    charindex = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0'
    if options:
        for idx, i in enumerate(options.stop_conditions):
            print(f"{i.tag}: {charindex[idx]}")


#### Stuff to try... automatically or by hand, line-by-line 

# takes a little time to run, but not too long.
print("--- Running Simulations ---")
o_exact, o_loose = setup_options_threewaybm()
# try this too.  toehold dissociates?... much slower!
# o_exact,o_loose = setup_options_threewaybm(bm_design="ACCGCACGTCCACGGTGTCG")

s=SimSystem(o_exact)
s.start()
s=SimSystem(o_loose)
s.start()
print("--- Finished simulations ---")


def print_results(o):
    print()
    print("--- Analysis of simulations by transitional states ---")
    # print("  Coarse-grained trajectory of simulation #1:")
    # print_transitions(o1.interface.transition_lists[0])
    print("  Transitions from simulation #1:")
    parsedlist = parse_transition_list(o.interface.transition_lists[0])
    print_transition_dict(parsedlist)
    print("  Transitions averaged over all %d simulations:" % o.num_simulations)
    parsedlist = parse_transition_lists(o.interface.transition_lists)
    print_transition_dict(parsedlist,o) # adds names for macrostates

print_results(o_exact)
print_results(o_loose)

# The lesson here is similar to hairpin_transitions.py:
# exact macrostates don't track system behavior well, because too much of the
# time the system is not in any macrostate -- thus it's easy to avoid the
# "checkpoints" and go straight from start to finish.
