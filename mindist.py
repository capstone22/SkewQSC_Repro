import itertools
from ffield import rref, dual_basis

def min_distance_fast(F, G, cap=None):
    H = dual_basis(F, G)
    if not H:
        return 1
    rho=len(H); n=len(H[0])
    cols=[[H[r][c] for r in range(rho)] for c in range(n)]
    maxw = cap if cap else rho+1
    for w in range(1, maxw+1):
        for S in itertools.combinations(range(n), w):
            sub=[cols[c] for c in S]
            rk=len(rref(F, sub)[1])
            if rk < w:
                return w
    return maxw+1

if __name__=="__main__":
    from ffield import GFTable
    F4=GFTable(2,[1,1,1])
    G=[[2,3,3,3,1,0,0,0],[0,3,2,2,2,1,0,0],[0,0,2,3,3,3,1,0],[0,0,0,3,2,2,2,1]]
    print("main example C1 d1 =",min_distance_fast(F4,G),"(expect 4)")
    Gc=[[1,0,0,0,1,0,0,0],[0,1,0,0,0,1,0,0],[0,0,1,0,0,0,1,0],[0,0,0,1,0,0,0,1]]
    print("cyclic [8,4] d =",min_distance_fast(F4,Gc),"(expect 2)")
