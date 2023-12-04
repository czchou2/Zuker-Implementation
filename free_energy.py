stacking_energies = {"GU": {"GU": -3, "AU": -3, "UA": -3, "CG": -13, "GC": -13},
                     "AU": {"GU": -3, "AU": -12, "UA": -18, "CG": -21, "GC": -21},
                     "UA": {"GU": -3, "AU": -18, "UA": -12, "CG": -21, "GC": -21},
                     "CG": {"GU": -13, "AU": -21, "UA": -21, "CG": -48, "GC": -43},
                     "GC": {"GU": -13, "AU": -21, "UA": -21, "CG": -30, "GC": -48}}

bulge_loop_energies = [28, 39, 45, 50, 52, 53, 55, 56, 57, 58, 59, 61, 62, 63, 64, 65, 67]

hairpin_loop_energies = {"CG": [float('inf'), float('inf'), 84, 59, 41, 43, 45, 46, 48, 49, 50, 52, 53, 54, 55, 57, 59],
                         "AU" : [float('inf'), float('inf'), 80, 75, 69, 64, 66, 68, 69, 70, 71, 73, 74,75,76,77]}
# interior_loop_energies = 
def eH(i,j,S):
    size = j - i - 2
    if size < 0:
        return float('inf')
    closing = str(S[i]) + str(S[j])
    closing = ''.join(sorted(closing))
    energy = hairpin_loop_energies[closing][size]
    return energy/10.0 # energy is in tenths of kcal/mole
def eL(i,j,ip, jp,S):
    # i paired with j, ip paired with jp inside i and j
    # size does not include the ending pairs
    size = j - i - 4
    if size < 0:
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
        pass

if __name__ == "__main__":
    print(eH(0,3,"AAUU"))
    print(eL(0,4,2,3,"GAAUC")) # example from Salser