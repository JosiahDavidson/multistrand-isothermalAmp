# Multistrand nucleic acid kinetic simulator
# Copyright (c) 2008-2024 California Institute of Technology. All rights reserved.
# The Multistrand Team (help@multistrand.org)


import os, random
import errno
from functools import reduce
from typing import Tuple, List

import numpy as np


def concentration_string(concentration):
    """
    An easy print function to format concentration in M
    """
    if concentration < 1e-12:
        return "{} fM".format(1e15 * concentration)
    if concentration < 1e-9:
        return "{} pM".format(1e12 * concentration)
    if concentration < 1e-6:
        return "{} nM".format(1e9 * concentration)
    if concentration < 1e-3:
        return "{} uM".format(1e6 * concentration)
    if concentration < 1:
        return "{} mM".format(1e3 * concentration)
    return "{} M".format(concentration)


def seqComplement(sequence):
    complement = {'G': 'C',
                  'C': 'G',
                  'A': 'T',
                  'T': 'A'}
    return "".join([complement[i] for i in reversed(sequence)])


def standardFileName(SCRIPT_DIR, mySeq=None, extraTitle=None, runs=None):
    fileName = str(SCRIPT_DIR)

    if not mySeq == None:
        fileName += str("/" + mySeq + '/' + mySeq)
    else:
        fileName += "/"

    if not runs == None:
        fileName += "-" + str(runs)
    if not extraTitle == None:
        fileName += "-" + extraTitle
    if not os.path.exists(os.path.dirname(fileName)):
        try:
            os.makedirs(os.path.dirname(fileName))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    return fileName


def generate_sequence(n, allowed_bases=['G', 'C', 'T', 'A'], base_probability=None):
    """ Generate a sequence of N base pairs.

    Bases are chosen from the allowed_bases [any sequence type], and
    according to the probability distribution base_probability - if
    none is specified, uses uniform distribution."""

    result = ""
    if base_probability is None:
        return result.join([random.choice(allowed_bases) for _ in range(n)])
    else:
        def uniform_seq(r):
            """ This function returns a lambda to be used on a sequence of tuples of
            (p,item) e.g. via reduce(uniform_seq(.75), seq, (0.0,'none')).

            In this situation, the lambda computes (p',item'), where
            p' is the sum of consecutive probabilities in the sequence
            <= r, and item' is the corresponding item.

            It achieves this by updating the result with the new sum
            of probabilities and last item checked until the the
            summed probability has exceeded the target.

            Note: r should be in [0.0,1.0) as produced by random(),
            any input r >= 1.0 (assuming the sequence given sums to
            1.0) will give the final element as the result which is technically incorrect.
            """
            return lambda x, y: (x[0] + y[0], y[1]) if r >= x[0] else (x[0], x[1])

        return result.join([
            reduce(uniform_seq(random.random()),
                   zip(base_probability, allowed_bases),
                   # note this subscript [1] pulls out the item selected by
                   # the reduce since the result was a tuple.
                   (0.0, 'Invalid Probabilities'))[1]
            for _ in range(n)])


def uniqueStateID(idsList, structsList):
    """
    Takes a list of ids and structures, and computes the pairtype for each of
    the complexes. Then returns a list of pairtypes that is alphabetically
    ordered.
    """
    pairTypes = []

    for ids, struct in zip(idsList, structsList):
        myPairType = pairType(ids, struct)
        pairTypes.append(myPairType)

    # now sort the list of lists by the first element of each list (which is
    mySortedList = sorted(pairTypes, key=lambda x: x[0])

    # to make this hashable, we make it into a tuple.
    return tuple(mySortedList)


def generatePairing(dotParen, stack, offset, output) -> None:
    """
    Utility function.
    """
    index = 1
    for c in dotParen:
        if c == '(':
            # pushing the first end of the basepair
            stack.append(offset + index)
        elif c == ')':
            # popping the stack, setting two locations
            currIndex = offset + index
            otherIndex = stack.pop()

            output[currIndex - 1] = otherIndex
            output[otherIndex - 1] = currIndex
        elif not c == '.':
            raise Warning('generatePairing: There is an error in the dot paren structure.')
        index += 1


def pairType(ids, structs):
    """
    Given identifiers and dot-parens for a complex,
    pairType returns a unique identifier for that secondary structure.
    """
    idList = ids.split(',')
    if all(id_.count(":") == 0 for id_ in idList):
        pass
    elif all(id_.count(":") == 1 for id_ in idList):
        nList = []
        for id_ in idList:
            id_ = id_.split(":")[1]
            nList.append(id_)
        idList = nList
    else:
        raise ValueError("Unsupported strand labelling.")

    dotParens = structs.split('+')
    N = len(dotParens)

    # the new ordering, for example: 3 0 1 2, so that idList[3] < idList[0] < idList[1] < idList[2]
    ordering = sorted(range(len(idList)), key=idList.__getitem__)

    idString = ''

    newLengths = [len(dotParens[ordering[i]]) for i in range(N)]
    newOffsets = [sum(newLengths[0:i]) for i in range(N)]  # the offsets under the new ordering
    offsets = [0, ] * N  # the new offsets under the old ordering
    for i in range(N):
        newPosition = ordering[i]
        offsets[newPosition] = newOffsets[i]
        idString += idList[newPosition]

    myStack = []
    output = [0, ] * sum([len(dp) for dp in dotParens])
    for index in range(len(idList)):
        generatePairing(dotParens[index], myStack, offsets[index], output)
    return (tuple(idString), tuple(output))


def normalizeCyclicPermutation(
    strandkey: List[int], seq: str, dpp: str) -> Tuple[str, str]:
    """
    Given the primary and secondary structure of a complex, both encoded as
    flat strings, return their minimal strand-level cyclic permutation in the
    same format. The order over cyclic permutations is induced by an order
    over strands, as defined by `strandkey`.

    Note:
      This is a purely syntactic transformation, and should leave the physical
      state unchanged.

      In general, the DPP syntax is equivariant w.r.t. the strand permutation
      *only* at intra-strand base pairs. Inter-strand base pairs are
      additionally flipped in orientation if they cross the index of the
      cyclic permutation.

    Comments:
      This function is intended for post-processing purposes, such as in
      `printTrajectory()`. It is semantically equivalent to
      `StrandOrdering::reorder()` in the C++ simulator, except that the latter
      is specialised to the strand order of a matching stop complex.
    """
    S = len(strandkey)
    assert S == len(np.unique(strandkey))
    assert S - 1 == seq.count('+') == dpp.count('+')

    # find minimal cyclic permutation
    prm_parts = lambda v: lambda p: (v[p:len(v)], v[0:p])
    prm = lambda v: lambda p: sum(prm_parts(v)(p), start=[])
    argprm = min(np.atleast_1d(np.argmin(strandkey)), key=prm(strandkey))

    # apply cyclic permutation
    if argprm == 0:
        # identity
        n_seq = seq
        n_dpp = dpp
    else:
        # transform primary structure
        n_seq = '+'.join(prm(seq.split('+'))(argprm))
        # transform intra-strand secondary structure
        n_dpp = [list('+'.join(p)) for p in prm_parts(dpp.split('+'))(argprm)]
        # transform inter-strand secondary structure
        n_dpp = '+'.join(map(''.join, map(flipInterPB, n_dpp, [True, False])))
        assert all(n_dpp.count(c) == dpp.count(c) for c in "().+")
    return n_seq, n_dpp


def flipInterPB(dpp: List[str], fw: bool) -> List[str]:
    """
    Detect inter-strand paired bases and flip their orientation in the DPP
    syntax, operating *in-place* and only in the designated direction.

    Used by `normalizeCyclicPermutation()`.
    """
    assert all(len(b) == 1 and b in "().+" for b in dpp)
    intra, n = 0, len(dpp)
    # recurse right/left on the new strand list's head/tail
    if fw:
        (p1, p2), bs = '()', enumerate(dpp)
    else:
        (p1, p2), bs = ')(', zip(range(n)[::-1], dpp[::-1])
    # find inter-strand paired bases with deprecated orientation
    for i, b in bs:
        if b == p1:
            intra += 1
        elif b == p2:
            if intra > 0:
                intra -= 1
            else:
                dpp[i] = p1
    return dpp


def printTrajectory(
    opts, timescale=(1e6, "us"),
    show_seed: bool=False, show_permutation: bool=False,
    show_uid: bool=False, feature=None):
    """
    Pretty-print the simulation trajectories stored in an `Options` object.
    The output representation of primary and secondary structures is normalized
    w.r.t. permutations among complexes and w.r.t. cyclic permutations of
    strands within each complex, using an ordering over strand lists based on
    the unique strand IDs.

    Custom output per simulation step, beyond the state seeds and
    complex/strand permutations, can be added via an argument of type
    `feature: Callable[[Options, int], str]`.
    """
    assert type(opts).__name__ == "Options"

    I, prev_tubeseq, prev_num_cmplx, prev_time = 0, "", 0, 0.0
    width_num_cmplx = len(str(max(map(len, opts.full_trajectory))))
    strands = opts.strand_names()

    for i in range(len(opts.full_trajectory)):
        curr_time = opts.full_trajectory_times[i]
        time = timescale[0] * curr_time
        state = opts.full_trajectory[i]

        # start of a new trajectory
        if i == 0 or curr_time < prev_time:
            if i > 0:
                # increment trajectory index
                I += 1
                # clean up overhang
                for prev_cmplx in opts.full_trajectory[i-1]:
                    if prev_cmplx in state:
                        state.remove(prev_cmplx)
            # print trajectory header
            if I > 0:
                print(); print()
            width_struct = sum(len(c.sequence) for c in state) + len(state) - 1
            t_col = f"  t[{timescale[1]}]  "
            dG_col = "dG[kcal/mol]"
            width_dG = len(dG_col)
            indent_hdr = width_num_cmplx + 3 + width_struct
            hdr = " | " + t_col + " | " + dG_col
            if show_seed:
                hdr += " | " + 4 * " " + "state_seed" + 5 * " "
            width_hdr = indent_hdr + len(hdr)
            print(width_hdr * "=")
            print(f"trajectory_seed = {opts.interface.results[I].seed}")
            print(width_hdr * "-")
            print(indent_hdr * " " + hdr)
            print(width_hdr * "-")

        # assemble normalized representation for each complex
        cids, seqs, structs, dG, pairTypes = [], [], [], 0.0, []
        for cmplx in state:
            sids = [strands[n].id for n in cmplx.strand_names.split(',')]
            cids.append(sids)
            sq, st = normalizeCyclicPermutation(
                sids, cmplx.sequence, cmplx.structure)
            seqs.append(sq)
            structs.append(st)
            dG += cmplx.energy
            if show_uid:
                uniqueID = pairType(cmplx.strand_names, cmplx.structure)
                pairTypes.append(
                    ''.join(uniqueID[0]) + '_' + ','.join(map(str, uniqueID[1])))

        # normalize the full state representation
        cids_ = sorted(cids)
        cixs = [cids_.index(cid) for cid in cids]
        tubeseq = ' '.join(np.array(seqs)[cixs])
        tubestruct = ' '.join(np.array(structs)[cixs])

        # print the full state on the next line(s)
        if show_permutation:
            # compare original and normalized complex states
            print(indent_hdr * "-")
            print(' # ' + width_num_cmplx * ' ' +
                  ' '.join((cmplx.sequence for cmplx in state)))
            print(f" {'#':<{width_num_cmplx}}  " +
                  f"{' '.join((cmplx.structure for cmplx in state))}")
            print(indent_hdr * "-")
        if i == 0 or tubeseq != prev_tubeseq:
            # display strand order when it changes
            print((width_num_cmplx + 3) * ' ' + tubeseq)
        if show_seed:
            seed = ",".join(f"{s:>5d}" for s in state[0].seed)
        if show_uid:
            identities = '+'.join(pairTypes[cixs])
        print(
            f"[{len(state):>{width_num_cmplx}}] " +
            f"{tubestruct} | {time:0<#9.4g} |   {dG:>+{width_dG-5}.3f}   " +
            ("" if not show_seed else f" | ({seed})") +
            ("" if feature is None else f" | {feature(opts, i)}") +
            ("" if not show_uid else f" | uID='{identities}'"))

        if opts.output_interval == 1 and tubeseq != prev_tubeseq:
            # if `normalizeCyclicPermutation()` is implemented correctly,
            # then this should happen *only* upon association/dissociation
            assert prev_num_cmplx != len(state)
        prev_time = curr_time
        prev_tubeseq, prev_num_cmplx = tubeseq, len(state)