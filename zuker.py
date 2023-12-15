import argparse
import timeit
from free_energy import eH, eS, eL


def parse(filename):
    """Parses input FASTA file for RNA sequence."""
    seq = ""
    with open(filename, "r", encoding="UTF-8") as file:
        lines = [line.rstrip() for line in file]
        assert len(lines) == 2, "File is not in FASTA format."
        seq = lines[1]
        file.close()
    return seq


def dbn_output(name, dbn, rna, dbnf):
    """Writes the RNA sequence and its secondary structure into a dot-bracket file.
    Args:
        name: Name of the bpRNA file
        dbn: Predicted secondary structure of an RNA sequence in dot-bracket notation
        rna: RNA sequence corresponding to the predicted structure
        dbn_file: Name of the .dbn file that is written to

    Returns:
        None.
    """
    with open(dbnf, "w",  encoding="UTF-8") as f:
        f.write("#Name: "+name + "\n")
        f.write("#Length: "+str(len(rna)) + "\n")
        f.write("#PageNumber: NA\n")
        f.write(rna + "\n")
        f.write(dbn)
        f.close()


def zuker(S):
    """Returns the 3 matrices constructed by Zuker's algorithm from an input RNA sequence."""
    n = len(S)

    a = 0  # energy contribution for closing of loop
    b = 0 # energy contribution for paired base (secondary structure) in multiloop
    c = 0 # energy contribution for unpaired base in multiloop

    # initialization (for j - i <= m, W(i,j)=0, V(i,j)=inf, WM(i,j)=inf):

    # minimal energy of general substructure i...j
    W = [[0 for _ in range(n)] for _ in range(n)]
    # minimal energy of closed substructure i...j
    V = [[float("inf") for _ in range(n)] for _ in range(n)]
    # minimal energy of true part of a multi-loop i...j
    WM = [[float("inf") for _ in range(n)] for _ in range(n)]

    for j in range(n):
        for i in reversed(range(j-m)):

            # Recursion equation for V:

            # energies from hairpin loop, stacking loop
            eq1 = min(eH(i, j, S), V[i+1][j-1] + eS(i, j, S))
            # energies from interior/bulge loop
            eq2 = min([V[ip][jp] + eL(i, j, ip, jp, S)
                       for ip in range(i+1, min(j, i+34))
                       for jp in range(max(ip+1, j-34), j)])
            # energies from multiloop
            eq3 = min([WM[i+1][k] + WM[k+1][j-1] + a
                      for k in range(i+1, j)])
            V[i][j] = min(eq1, eq2, eq3)

            # Recursion equation for WM:

            # energies from leaving j unpaired, leaving i unpaired, pairing i and j
            eq4 = min(WM[i][j-1] + c, WM[i+1][j] + c, V[i][j] + b)
            # energies from non-closed multiloop
            eq5 = min([WM[i][k] + WM[k+1][j] for k in range(i+1, j)])
            WM[i][j] = min(eq4, eq5)

    for j in range(n):
        for i in reversed(range(j-m)):

            # Recursion equation for W:

            eq0 = min([W[i][k-1] + V[k][j] if k > 0 else V[k][j]
                      for k in range(i, j-m)])
            # energies from j unpaired, j paired
            W[i][j] = min([W[i][j-1], eq0])

    return W, V, WM


# x = 0 (W), 1 (V), WM (2)
def backtrace(i, j, W, V, WM, x, dbn, S):
    """Reconstructs the secondary structure predicted by Zuker's algorithm.
    Args:
        i, j: Start and end indices of the sequence to be backtraced
        W, V, WM: Matrices computed from Zuker's algorithm
        x: Matrix the backtrace is currently recursing on (0 = W, 1 = V, 2 = WM)
        dbn: The dot-bracket string to be modified
        S: The RNA sequence for which the structure is predicted

    Returns: The dot-bracket string dbn with pairs added.
    """
    a = 0
    b = 0
    c = 0
    m = 3

    if j - m > i:
        # backtrack matrix W
        if x == 0:
            if W[i][j] == W[i][j-1]:
                # j unpaired
                dbn = backtrace(i, j-1, W, V, WM, 0, dbn, S)
            else:
                score, k = min([(W[i][k-1] + V[k][j], k) if k >
                               0 else (V[k][j], k) for k in range(i, j-m)])

                if W[i][j] == score:
                    # j paired
                    dbn = backtrace(i, k-1, W, V, WM, 0, dbn, S)
                    dbn = backtrace(k, j, W, V, WM, 1, dbn, S)
        # backrack matrix V
        elif x == 1:
            if V[i][j] == eH(i, j, S):
                # hairpin loop
                dbn = dbn[:i] + "(" + dbn[i+1:j] + ")" + dbn[j+1:]
            elif V[i][j] == V[i+1][j-1] + eS(i, j, S):
                # stacking loop
                dbn = dbn[:i] + "(" + dbn[i+1:j] + ")" + dbn[j+1:]
                dbn = backtrace(i+1, j-1, W, V, WM, 1, dbn, S)
            else:
                score, ip, jp = min([(V[ip][jp] + eL(i, j, ip, jp, S), ip, jp)
                                     for ip in range(i+1, min(j, i+34))
                                     for jp in range(max(ip+1, j-34), j)])
                if V[i][j] == score:
                    # interior or bulge loop
                    dbn = dbn[:i] + "(" + dbn[i+1:j] + ")" + dbn[j+1:]
                    dbn = backtrace(ip, jp, W, V, WM, 1, dbn, S)
                else:
                    # multiloop
                    score, k = min([(WM[i+1][k] + WM[k+1][j-1] + a, k)
                                    for k in range(i+1, j)])
                    dbn = dbn[:i] + "(" + dbn[i+1:j] + ")" + dbn[j+1:]
                    dbn = backtrace(i+1, k, W, V, WM, 2, dbn, S)
                    dbn = backtrace(k+1, j-1, W, V, WM, 2, dbn, S)

        # backtrack matrix WM
        else:
            if WM[i][j] == WM[i][j-1] + c:
                # j unpaired
                dbn = backtrace(i, j-1, W, V, WM, 2, dbn, S)
            elif WM[i][j] == WM[i+1][j] + c:
                # i unpaired
                dbn = backtrace(i+1, j, W, V, WM, 2, dbn, S)
            elif WM[i][j] == V[i][j] + b:
                # closed
                dbn = backtrace(i, j, W, V, WM, 1, dbn, S)
            else:
                # non-closed
                score, k = min([(WM[i][k] + WM[k+1][j], k)
                                for k in range(i+1, j)])
                dbn = backtrace(i, k, W, V, WM, 2, dbn, S)
                dbn = backtrace(k+1, j, W, V, WM, 2, dbn, S)

    return dbn


def main(ff, dbnf, lim):
    seq = parse(ff)

    if lim > -1:
        print("limiting length to", lim, "...")
        S = seq[:lim]
    else:
        S = seq

    starttime = timeit.default_timer()

    W, V, WM = zuker(S)
    print("Minimum energy score:", W[0][len(S)-1],
          "\nTime:", timeit.default_timer() - starttime)
    dbn = "-" * len(S)
    dbn = backtrace(0, len(S)-1, W, V, WM, 0, dbn, S)
    dbn_output("test", dbn, S, dbnf)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help="fasta file", required=True)
    parser.add_argument("-d", "--dbn", help="dbn file", required=True)
    parser.add_argument(
        "-l", "--limit", help="limit length of sequence", default=-1, type=int)
    args = parser.parse_args()
    fasta_file = args.fasta
    dbn_file = args.dbn
    limit = args.limit
    main(fasta_file, dbn_file, limit)
