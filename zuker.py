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


def eH(i, j):
    '''
    input: start and end position
    output: energy of hairpin loop from sequence i..j
    '''
    pass


def eS(i, j):
    '''
    input: start and end position
    output: energy of stacking loop from sequence i..j
    '''
    pass


def eL(i, j, i_prime, j_prime):
    '''
    input:
    output:
    '''
    pass


def dp(seq):
    n = len(seq)
    # minimal energy of general substructure i...j
    w = [[0 for _ in range(n+1)] for _ in range(n+1)]
    # minimal energy of closed substructure i...j
    v = [[0 for _ in range(n+1)] for _ in range(n+1)]
    # minimal energy of true part of a multi-loop i...j
    wm = [[0 for _ in range(n+1)] for _ in range(n+1)]

    for i in range(n):
        for j in range(n):
            temp = min([w[i][k-1] + v[k][j] for k in range(i, j-m-1)])
            w[i][j] = min(w[i][j-1], temp)

            eq1 = min(eH(i, j), v[i+1][j-1] + eS(i, j))
            eq2 = min([v[i_prime][j_prime] + eL(i, j, i_prime, j_prime)
                       for i_prime in range(i+1, j) for j_prime in range(i_prime+1, j)])
            eq3 = min([wm[i+1][k] + wm[k+1][j-1] + a for k in range(i+1, j)])
            v[i][j] = min(eq1, eq2, eq3)

            eq4 = min(wm[i][j-1] + c, wm[i+1][j] + c, v[i][j] + b)
            eq5 = min([wm[i][k] + wm[k+1][j] for k in range(i+1, j)])
            wm[i][j] = min(eq4, eq5)

    return w[1][n]


def backtrace():
    pass


def main():
    file = "bpRNA_CRW_296.fasta"
    print(parse(file))


if __name__ == "__main__":
    main()
