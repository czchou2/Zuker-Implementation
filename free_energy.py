stacking_energies = {"GU": {"GU": -3, "AU": -3, "UA": -3, "CG": -13, "GC": -13},
                     "AU": {"GU": -3, "AU": -12, "UA": -18, "CG": -21, "GC": -21},
                     "UA": {"GU": -3, "AU": -18, "UA": -12, "CG": -21, "GC": -21},
                     "CG": {"GU": -13, "AU": -21, "UA": -21, "CG": -48, "GC": -43},
                     "GC": {"GU": -13, "AU": -21, "UA": -21, "CG": -30, "GC": -48}}

# Some energies are interpolated
bulge_loop_energies = [28, 39, 45, 50, 52, 53, 55, 56, 57, 58,58, 59,60, 61, 61, 62,62, 63, 63, 64,64,64,64,64, 65,65,65,66,66, 67]

hairpin_loop_energies = {"CG": [float('inf'), float('inf'), 84, 59, 41, 43, 45, 46, 48, 49,49, 50,51, 52,52, 53,53, 54,54, 55,55,55,56,56, 57,57,57,58,58, 59],
                         "AU" : [float('inf'), float('inf'), 80, 75, 69, 64, 66, 68, 69, 70,70, 71,72, 73,73, 74,74,75,75,76,76,76,76,76,77,77,77,78,78,79]}
interior_loop_energies = {"CG-CG": [float('inf'), 1, 9, 16, 21, 25, 26, 27, 28, 29,30, 31,31, 32,32, 33,33, 34,34, 35,35,35,36,36, 37,37,37,38,38, 39], 
                          "CG-AU":[float('inf'), 10, 18, 25, 30, 34, 35, 36, 37, 38,38, 39,39, 40,40, 41,41, 42,42, 43,43,43,44,44, 45,45,45,46,46, 47],
                          "AU-CG":[float('inf'), 10, 18, 25, 30, 34, 35, 36, 37, 38,38, 39,39, 40,40, 41,41, 42,42, 43,43,43,44,44, 45,45,45,46,46, 47],
                          "AU-AU":[float('inf'), 18, 26, 33, 38, 42, 43, 44, 45, 46,47, 48,48, 49,49, 50,50, 51,51, 52,52,52,53,53, 54,54,54,55,55, 56]}
def eH(i,j,S):
    size = j - i - 2
    if size < 0 or size > 29:
        return float('inf')
    closing = str(S[i]) + str(S[j])
    closing = ''.join(sorted(closing))
    energy = hairpin_loop_energies[closing][size]
    return energy/10.0 # energy is in tenths of kcal/mole
def eL(i,j,ip, jp,S):
    # i paired with j, ip paired with jp inside i and j
    # size does not include the ending pairs
    size = j - i - 4
    if size < 0 or size > 29:
        return float('inf')
    if(i+1 == ip or j-1 == jp):
        #BULGE
        closing_ext = str(S[i]) + str(S[j])
        if(closing_ext == "UG"):
            closing_ext = "GU"
        closing_int = str(S[ip]) + str(S[jp])
        if(closing_int == "UG"):
            closing_int = "GU"
        return (bulge_loop_energies[size] + stacking_energies[closing_ext][closing_int])/10.0
    else:
        #INTERNAL
        closing = ''.join(sorted(str(S[i])+str(S[j]))) + "-" + ''.join(sorted(str(S[ip])+str(S[jp])))
        return interior_loop_energies[closing][size]/10.0


if __name__ == "__main__":
    print(eH(0,3,"AAUU"))
    print(eL(0,4,2,3,"GAAUC")) # example from Salser
    print(eL(0,5,2,3,"GAAUAC")) 
    print(len(hairpin_loop_energies["CG"]))
    print(len(hairpin_loop_energies["AU"]))
    print(len(interior_loop_energies["AU-AU"]))