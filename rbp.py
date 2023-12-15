def deltabb(i,j,ip,jp):
    return max(abs(i-ip),abs(j-jp))
def deltasb(i,j,S):
    # make stack of indices of opening pairs, once hit closing pair pop index and that's ip,jp
    stack = []
    mini = None
    for n in range(len(S)):
        if(S[n] == '('):
            stack.append(n)
        elif(S[n]==')'):
            ip = stack.pop()
            jp = n
            dbb = deltabb(i,j,ip,jp)
            if(mini == None or dbb<mini):
                mini = dbb
    assert(mini != None)
    return mini
def distList(S1,S2):
    stack = []
    dsblist = []
    for n in range(len(S1)):
        if(S1[n]=='('):
            stack.append(n)
        elif(S1[n]==')'):
            i = stack.pop()
            j = n
            dsblist.append(deltasb(i,j,S2))
    return dsblist
def relaxed_base_pair(S1,S2,t):
    '''
    input: 2 dbn sequences, t relaxation parameter
    output: relaxed base pair distance
    ASSUMES: there is at least one pair in each dbn sequence and there are no unclosed parentheses.
    '''

    dsb = distList(S1,S2)
    dsb += distList(S2,S1)
    dsb = sorted(dsb,reverse=True)
    print(dsb)
    M = len(dsb)
    minm = M
    for m in range(M-1):
        if(dsb[m]<= t * m):
            if(m < minm):
                minm = m
                break
    return minm

if __name__ == "__main__":
    # t = 0 corresponds to regular base pair, increase to relax it
    rbp = relaxed_base_pair("(((...(.((()))..)..).))","((...(.().().)...))()()",0)
    print(rbp)