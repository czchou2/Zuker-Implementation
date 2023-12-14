from free_energy import eH, eS, eL
import sys
import argparse
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
    b = 0 # energy contribution for paired base (secondary structure) in multiloop
    c = 0 # energy contribution for unpaired base in multiloop

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
    for j in range(n):
        for i in reversed(range(j-m)):

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

    for j in range(n):
        for i in reversed(range(j-m)):
            # Recursion equation for W
            eq0 = min([W[i][k-1] + V[k][j] if k>0 else V[k][j] for k in range(i,j-m)])
            # corresponds to j unpaired, j paired
            W[i][j] = min([W[i][j-1], eq0])

    return W, V, WM



# x = 0 (W), 1 (V), WM (2)
def backtrace(i, j, W, V, WM, x,dbn,S):
    '''
    input: i index, j index, W matrix, V matrix, WM matrix, x (matrix recursing on), dbn (dot bracket to be modified), S rna sequence
    output: dbn with pairs added
    '''
    # temp
    m = 3
    a = 0
    b = 0
    c = 0

    if j - m > i: 
        # backtrack matrix W
        if x == 0:
            if W[i][j] == W[i][j-1]:
                # j unpaired
                dbn = backtrace(i, j-1, W, V, WM, 0,dbn,S)
            else:
                score, k = min([(W[i][k-1] + V[k][j],k) if k>0 else (V[k][j],k) for k in range(i,j-m)])

                if W[i][j] == score:
                    # j paired
                    dbn = backtrace(i, k-1, W, V, WM, 0,dbn,S)
                    dbn = backtrace(k, j, W, V, WM, 1,dbn,S)
        # backrack matrix V
        elif x == 1:
            if V[i][j] == eH(i, j,S):
                # hairpin loop
                dbn = dbn[:i] + "(" + dbn[i+1:j] + ")" + dbn[j+1:]
            elif V[i][j] == V[i+1][j-1] + eS(i, j, S):
                # stacking loop
                dbn = dbn[:i] + "(" + dbn[i+1:j] + ")" + dbn[j+1:]
                dbn = backtrace(i+1, j-1, W, V, WM, 1,dbn,S)
            else:
                score, ip, jp = min([(V[ip][jp] + eL(i, j, ip, jp, S), ip, jp)
                       for ip in range(i+1, j)
                       for jp in range(ip+1, j)])
                if V[i][j] == score:
                    # interior loop or bulge
                    dbn = dbn[:i] + "(" + dbn[i+1:j] + ")" + dbn[j+1:]
                    dbn = backtrace(ip, jp, W, V, WM, 1,dbn,S)
                else:
                    # multiloop
                    score, k = min([(WM[i+1][k] + WM[k+1][j-1] + a, k)
                      for k in range(i+1, j)])
                    dbn = dbn[:i] + "(" + dbn[i+1:j] + ")" + dbn[j+1:]
                    dbn = backtrace(i+1, k, W, V, WM, 2,dbn,S)
                    dbn = backtrace(k+1, j-1, W, V, WM, 2,dbn,S)

        # backtrack matrix WM
        else:
            if WM[i][j] == WM[i][j-1] + c:
                # j unpaired
                dbn = backtrace(i, j-1, W, V, WM, 2,dbn,S)
            elif WM[i][j] == WM[i+1][j] + c:
                # i unpaired
                dbn = backtrace(i+1, j, W, V, WM, 2,dbn,S)
            elif WM[i][j] == V[i][j] + b:
                # closed
                dbn = backtrace(i,j,W,V,WM,1,dbn,S)
            else:
                # unclosed
                score, k = min([(WM[i][k] + WM[k+1][j], k)
                      for k in range(i+1, j)])
                dbn = backtrace(i,k,W,V,WM,2,dbn,S)
                dbn = backtrace(k+1,j,W,V,WM,2,dbn,S)

                
    return dbn


def main(ff,dbnf,lim):
    # file = "bpRNA_CRW_296.fasta"
    seq = parse(ff)

    if(lim > -1):
        print("limiting length to",lim,"...")
        S = seq[:lim]
    else:
        S = seq

    starttime = timeit.default_timer()

    W, V, WM = zuker(S)
    print("Minimum energy score:",W[0][len(S)-1],"\nTime:", timeit.default_timer() - starttime)
    dbn = "-" * len(S)
    dbn = backtrace(0, len(S)-1, W, V, WM, 0,dbn,S)
    dbn_output("test",dbn,S,dbnf)

if __name__ == "__main__":
    # fasta_file = sys.argv[1]
    # dbn_file = sys.argv[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="fasta file", required=True)
    parser.add_argument("-d", "--dbn", help="dbn file", required=True)
    parser.add_argument("-l", "--limit", help="limit length of sequence", default=-1, type=int)
    args = parser.parse_args()
    fasta_file = args.fasta
    dbn_file = args.dbn
    limit = args.limit
    main(fasta_file,dbn_file,limit)