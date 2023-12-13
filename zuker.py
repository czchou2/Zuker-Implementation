from free_energy import eH, eS, eL
# import sys
import timeit


def parse(filename):
    """Parses input FASTA file for RNA sequence."""
    seq = ""
    with open(filename, "r", encoding="UTF-8") as file:
        lines = [line.rstrip() for line in file]
        assert (len(lines) == 2)
        seq = lines[1]
        file.close()
    return seq


def dbn_output(name, dbn, rna, dbn_file):
    '''
    input: name string, dot-bracket string, sequence string, file path
    '''
    f = open(dbn_file, "w")
    f.write("#Name: "+name + "\n")
    f.write("#Length: "+str(len(rna)) + "\n")
    f.write("#PageNumber: NA\n")
    f.write(rna + "\n")
    f.write(dbn)
    f.close()


# def a(i, j, S):
#     """Returns -1 or 0 depending on whether or not a base pairing between the ith and jth nucleotides in S is allowed or not."""
#     if S[i] == 'A' and S[j] == 'U' \
#             or S[i] == 'U' and S[j] == 'A' \
#             or S[i] == 'C' and S[j] == 'G' \
#             or S[i] == 'G' and S[j] == 'C':
#         return -1
#     return 0


def zuker(S):
    n = len(S)

    # TODO:
    a = 0  # energy contribution for closing of loop
    b = 0
    c = 0

    # minimal loop size (slide 7)
    m = 3

    # minimal energy of general substructure i...j
    W = [[0 for _ in range(n)] for _ in range(n)]
    # minimal energy of closed substructure i...j
    V = [[float("inf") for _ in range(n)] for _ in range(n)]
    # minimal energy of true part of a multi-loop i...j
    WM = [[float("inf") for _ in range(n)] for _ in range(n)]

    # initialization:
    # for j - i <= m, W(i,j)=0, V(i,j)=inf, WM(i,j)=inf

    # recursion:
    # for i < j - m
    for i in range(n):
        for j in range(n):
            if i >= j - m:
                continue

            # Recursion equation for W
            eq0 = min([W[i][k-1] + V[k][j]
                       for k in range(i, j-m)])
            # corresponds to j unpaired, j paired
            W[i][j] = min(W[i][j-1], eq0)

            # Recursion equation for V
            # corresponds to:
            # hairpin loop, stacking loop
            # interior loop/bulge
            # multi-loop
            eq1 = min(eH(i, j, S), V[i+1][j-1] + eS(i, j, S))
            eq2 = min([V[ip][jp] + eL(i, j, ip, jp, S)
                       for ip in range(i+1, j)
                       for jp in range(ip+1, j)])
            eq3 = min([WM[i+1][k] + WM[k+1][j-1] + a
                      for k in range(i+1, j)])
            V[i][j] = min(eq1, eq2, eq3)

            # Recursion equation for WM
            # corresponds to:
            # j unpaired, i unpaired, closed
            # non-closed
            eq4 = min(WM[i][j-1] + c, WM[i+1][j] + c, V[i][j] + b)
            eq5 = min([WM[i][k] + WM[k+1][j] for k in range(i+1, j)])
            WM[i][j] = min(eq4, eq5)

    return W, V, WM


# x = 0 (W), 1 (V), WM (2)
def backtrace(i, j, W, V, WM, x):
    # temp
    m = 3
    a = 0
    b = 0
    c = 0

    i = 1
    j = len(W)
    P = []
    while j > 1:

        # backtrack matrix W
        if x == 0:
            if W[i][j] == W[i][j-1]:
                # j unpaired
                j -= 1
                continue

            score, k = min([(W[i][k-1] + V[k][j], k)
                            for k in range(i, j-m)])
            if W[i][j] == score:
                # j paired
                backtrace(i, k-1, W, V, WM, 0)

                backtrace(k, j, W, V, WM, 1)
                # P.append((j, k))

        # backrack matrix V
        elif x == 1:
            if V[i][j] == eH(i, j):
                # hairpin loop
                # done?
                return
            elif V[i][j] == V[i+1][j-1] + eS(i, j):
                # stacking loop
                i += 1
                j -= 1
            else:
                score, ip, jp = min([(V[ip][jp] + eL(i, j, ip, jp, S), ip, jp)
                       for ip in range(i+1, j)
                       for jp in range(ip+1, j)])
                if V[i][j] == score:
                    # interior loop or bulge
                    backtrace(ip, jp, W, V, WM, 1)
                    # append P?

                else:
                    # multiloop
                    score, k = min([(WM[i+1][k] + WM[k+1][j-1] + a, k)
                      for k in range(i+1, j)])
                    backtrace(i+1, k, W, V, WM, 2)
                    backtrace(k+1, j-1, W, V, WM, 2)
                    # append P?
                    

        # backtrack matrix WM
        else:
            if WM[i][j] == WM[i][j-1] + c:
                # j unpaired
            elif WM[i][j] == WM[i+1][j] + c:
                # i unpaired
            elif WM[i][j] == V[i][j] + b:
                # closed
            else:
                # unclosed
                # backtrace

    return P


def main():
    file = "bpRNA_CRW_296.fasta"
    seq = parse(file)

    print(len(seq))

    S = seq[:150]
    starttime = timeit.default_timer()
    W, V, WM = zuker(S)
    print(W[0][len(S)-1], timeit.default_timer() - starttime)

    structure = backtrace(0, len(S)-1, W, V, WM, 0)
    # print(structure)

    # dbn_file = sys.argv[1]
    # dbn = "((..))"
    # rna = "AAGATT"
    # dbn_output("test",dbn,rna,dbn_file)


if __name__ == "__main__":
    main()
