# Reference: RNA SECONDARY STRUCTURES AND THEIR PREDICTION (pg 603, table 1)
# (Energies are in tenths of a kcal/mole)

"""
A dictionary where stacking_energies[x][y] is the energy contributed
by the exterior closing pair x and interior closing pair y.
"""
stacking_energies = {"GU": {"GU": -3, "AU": -3, "UA": -3, "CG": -13, "GC": -13},
                     "AU": {"GU": -3, "AU": -12, "UA": -18, "CG": -21, "GC": -21},
                     "UA": {"GU": -3, "AU": -18, "UA": -12, "CG": -21, "GC": -21},
                     "CG": {"GU": -13, "AU": -21, "UA": -21, "CG": -48, "GC": -43},
                     "GC": {"GU": -13, "AU": -21, "UA": -21, "CG": -30, "GC": -48}}

"""
A list where bulge_loop_energies[i] is the energy contributed by a loop of size i+1.
(Some energies are interpolated)
"""
bulge_loop_energies = [28, 39, 45, 50, 52, 53, 55, 56, 57, 58, 58, 59,
                       60, 61, 61, 62, 62, 63, 63, 64, 64, 64, 64, 64, 65, 65, 65, 66, 66, 67]

"""
A dictionary where hairpin_loop_energies[x][i] is the energy contributed
by a loop of size i+1 closed by x.
"""
hairpin_loop_energies = {"CG": [float('inf'), float('inf'), 84, 59, 41, 43, 45, 46, 48, 49, 49, 50, 51, 52, 52, 53, 53, 54, 54, 55, 55, 55, 56, 56, 57, 57, 57, 58, 58, 59],
                         "AU": [float('inf'), float('inf'), 80, 75, 69, 64, 66, 68, 69, 70, 70, 71, 72, 73, 73, 74, 74, 75, 75, 76, 76, 76, 76, 76, 77, 77, 77, 78, 78, 79]}

"""
A dictionary where interior_loop_energies[x][i] is the energy contributed
by a loop of size i+1 closed by x.
"""
interior_loop_energies = {"CG-CG": [float('inf'), 1, 9, 16, 21, 25, 26, 27, 28, 29, 30, 31, 31, 32, 32, 33, 33, 34, 34, 35, 35, 35, 36, 36, 37, 37, 37, 38, 38, 39],
                          "CG-AU": [float('inf'), 10, 18, 25, 30, 34, 35, 36, 37, 38, 38, 39, 39, 40, 40, 41, 41, 42, 42, 43, 43, 43, 44, 44, 45, 45, 45, 46, 46, 47],
                          "AU-CG": [float('inf'), 10, 18, 25, 30, 34, 35, 36, 37, 38, 38, 39, 39, 40, 40, 41, 41, 42, 42, 43, 43, 43, 44, 44, 45, 45, 45, 46, 46, 47],
                          "AU-AU": [float('inf'), 18, 26, 33, 38, 42, 43, 44, 45, 46, 47, 48, 48, 49, 49, 50, 50, 51, 51, 52, 52, 52, 53, 53, 54, 54, 54, 55, 55, 56]}


def eH(i, j, S):
    """Calculates the energy of a hairpin loop from S[i..j]."""
    size = j - i - 2
    if size < 0 or size > 29:
        return float('inf')
    closing = str(S[i]) + str(S[j])
    closing = ''.join(sorted(closing))

    if closing not in hairpin_loop_energies:
        return float('inf')
    energy = hairpin_loop_energies[closing][size]
    return energy/10.0


def eL(i, j, ip, jp, S):
    """Calculates the energy of a bulge or interior loop formed by S[i..j] and S[i'...j'],
    where (i,j) is the exterior closing pair and (i',j') is the interior closing pair."""
    # i paired with j, ip paired with jp inside i and j
    # size does not include the ending pairs
    size = ip - i + j - jp - 3
    if size < 0 or size > 29:
        return float('inf')
    if i+1 == ip or j-1 == jp:
        # BULGE
        closing_ext = str(S[i]) + str(S[j])
        if closing_ext == "UG":
            closing_ext = "GU"
        closing_int = str(S[ip]) + str(S[jp])
        if closing_int == "UG":
            closing_int = "GU"

        if (closing_int not in stacking_energies or closing_ext not in stacking_energies):
            return float('inf')
        return (bulge_loop_energies[size] + stacking_energies[closing_ext][closing_int])/10.0
    else:
        # INTERIOR
        closing = ''.join(sorted(str(S[i])+str(S[j]))) + \
            "-" + ''.join(sorted(str(S[ip])+str(S[jp])))

        if closing not in interior_loop_energies:
            return float('inf')
        return interior_loop_energies[closing][size]/10.0


def eS(i, j, S):
    """Calculates the energy of a stacking loop from S[i..j]."""
    ext_closing = str(S[i]) + str(S[j])
    if ext_closing == "UG":
        ext_closing = "GU"

    int_clos_i = i + 1
    int_clos_j = j - 1

    int_closing = str(S[int_clos_i]) + str(S[int_clos_j])
    if int_closing == "UG":
        int_closing = "GU"
    if (ext_closing not in stacking_energies) or (int_closing not in stacking_energies[ext_closing]):
        return float('inf')
    energy = stacking_energies[ext_closing][int_closing]
    return energy/10.0


if __name__ == "__main__":
    # Test hairpin loop
    assert eH(0, 3, "AAUU") == hairpin_loop_energies["AU"][1]/10.0

    # Test bulge loop
    assert eL(0, 4, 2, 3, "GAAUC") == 0.7  # example from Salser

    # Test interior loop
    assert eL(0, 5, 2, 3, "GAAUAC") == interior_loop_energies["AU-CG"][1]/10.0

    # Test stacking loop
    assert eS(0, 3,"GAUC") == stacking_energies["GC"]["AU"]/10.0
