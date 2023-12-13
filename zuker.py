from free_energy import *
seq = ""


def parse(filename):
    '''
    input: FASTA file
    output: RNA sequence
    '''
    seq = ""
    with open(filename, "r", encoding="UTF-8") as file:
        lines = [line.rstrip() for line in file]
        assert (len(lines) == 2)
        seq = lines[1]
        file.close()
    return seq

def dbn_output(name,dbn,rna,dbn_file):
    '''
    input: name string, dot-bracket string, sequence string, file path
    '''
    f = open(dbn_file,"w")
    f.write("#Name: "+name + "\n")
    f.write("#Length: "+str(len(rna)) +"\n")
    f.write("#PageNumber: NA\n")
    f.write(rna + "\n")
    f.write(dbn)
    f.close()

def a(i, j):
    '''
    input: start and end position
    output: -1 or 0 depending on whether or not a base pairing
            between the ith and jth nucleotides is allowed or not
    '''
    if seq[i] == 'A' and seq[j] == 'U' \
            or seq[i] == 'U' and seq[j] == 'A' \
            or seq[i] == 'C' and seq[j] == 'G' \
            or seq[i] == 'G' and seq[j] == 'C':
        return -1
    return 0


def dp(seq):
    n = len(seq)

    # page 609
    k = 1
    b = 0
    c = 0

    m = 30
    # bound on maximal interior loop size?

    # minimal energy of general substructure i...j
    w = [[0 for _ in range(n+1)] for _ in range(n+1)]
    # minimal energy of closed substructure i...j
    v = [[float("inf") for _ in range(n+1)] for _ in range(n+1)]
    # minimal energy of true part of a multi-loop i...j
    wm = [[float("inf") for _ in range(n+1)] for _ in range(n+1)]

    # (n+1) x (n+1) matrix of entries (i, j, matrix)
    bt = [[None for _ in range(n+1)] for _ in range(n+1)]

    for i in range(n):
        for j in range(n):
            temp, pos = min([(w[i][k-1] + v[k][j], k)
                            for k in range(i, j-m-1)])
            w[i][j] = min(w[i][j-1], temp)
            bt[i][j] = 

            eq1 = min(eH(i, j, seq), v[i+1][j-1] + eS(i, j, seq))
            eq2 = min([v[i_prime][j_prime] + eL(i, j, i_prime, j_prime, seq)
                       for i_prime in range(i+1, j)
                       for j_prime in range(i_prime+1, j)])
            eq3 = min([wm[i+1][k] + wm[k+1][j-1] + a(i, j)
                      for k in range(i+1, j)])
            v[i][j] = min(eq1, eq2, eq3)

            eq4 = min(wm[i][j-1] + c, wm[i+1][j] + c, v[i][j] + b)
            eq5 = min([wm[i][k] + wm[k+1][j] for k in range(i+1, j)])
            wm[i][j] = min(eq4, eq5)

    return w[1][n], bt


def backtrace(bt):
    pass


def main():
    file = "bpRNA_CRW_296.fasta"
    seq = parse(file)
    score, bt = dp(seq)
    structure = backtrace(bt)


if __name__ == "__main__":
    main()
