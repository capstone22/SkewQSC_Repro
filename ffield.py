import itertools

class GFTable:
    def __init__(self, p, mod):
        self.p=p; self.r=len(mod)-1; self.q=p**self.r; self.mod=mod
        self.elems=list(range(self.q))
        self._tuple=[self._int2tuple(i) for i in range(self.q)]
        self._idx={t:i for i,t in enumerate(self._tuple)}
        self.ADD=[[0]*self.q for _ in range(self.q)]
        self.MUL=[[0]*self.q for _ in range(self.q)]
        self.NEG=[0]*self.q
        for a in range(self.q):
            ta=self._tuple[a]
            self.NEG[a]=self._idx[tuple((-x)%p for x in ta)]
            for b in range(self.q):
                tb=self._tuple[b]
                self.ADD[a][b]=self._idx[tuple((x+y)%p for x,y in zip(ta,tb))]
                self.MUL[a][b]=self._idx[self._polymul(ta,tb)]
        self.INV=[0]*self.q
        for a in range(1,self.q):
            for b in range(1,self.q):
                if self.MUL[a][b]==1: self.INV[a]=b; break
        self.ONE=self._idx[tuple([1]+[0]*(self.r-1))]
        self.ZERO=0
    def _int2tuple(self,i):
        d=[]; 
        for _ in range(self.r): d.append(i%self.p); i//=self.p
        return tuple(d)
    def _polymul(self,ta,tb):
        p=self.p; r=self.r; res=[0]*(2*r-1)
        for i in range(r):
            if ta[i]==0: continue
            for j in range(r):
                res[i+j]=(res[i+j]+ta[i]*tb[j])%p
        for i in range(2*r-2,r-1,-1):
            if res[i]!=0:
                c=res[i]
                for k in range(r+1):
                    res[i-r+k]=(res[i-r+k]-c*self.mod[k])%p
        return tuple(res[:r])
    def add(self,a,b): return self.ADD[a][b]
    def mul(self,a,b): return self.MUL[a][b]
    def neg(self,a): return self.NEG[a]
    def sub(self,a,b): return self.ADD[a][self.NEG[b]]
    def inv(self,a): return self.INV[a]
    def frob(self,a,power):
        e=pow(self.p,power)
        return self.powe(a,e)
    def powe(self,a,e):
        res=self.ONE; base=a
        while e>0:
            if e&1: res=self.MUL[res][base]
            base=self.MUL[base][base]; e>>=1
        return res

def rref(F, rows):
    M=[r[:] for r in rows]; 
    if not M: return [],[]
    ncol=len(M[0]); piv=[]; r=0
    for c in range(ncol):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]!=F.ZERO: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        inv=F.inv(M[r][c])
        M[r]=[F.mul(inv,x) for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]!=F.ZERO:
                f=M[i][c]
                M[i]=[F.sub(a,F.mul(f,b)) for a,b in zip(M[i],M[r])]
        piv.append(c); r+=1
        if r==len(M): break
    return M[:r], piv

def rank(F,rows): return len(rref(F,rows)[1])

def in_rowspace(F, basis_rref, piv, v):
    w=v[:]
    for row,c in zip(basis_rref,piv):
        if w[c]!=F.ZERO:
            f=w[c]
            w=[F.sub(a,F.mul(f,b)) for a,b in zip(w,row)]
    return all(x==F.ZERO for x in w)

def dual_basis(F, G):
    R,piv=rref(F,G); k=len(piv); 
    n=len(G[0]); free=[c for c in range(n) if c not in piv]
    basis=[]
    pivset={c:i for i,c in enumerate(piv)}
    for fc in free:
        v=[F.ZERO]*n; v[fc]=F.ONE
        for c,i in pivset.items():
            v[c]=F.neg(R[i][fc])
        basis.append(v)
    return basis

def min_distance(F, G, cap=None):
    R,piv=rref(F,G); k=len(piv)
    if k==0: return 0
    n=len(R[0]); q=F.q; best=n+1
    for coeffs in itertools.product(range(q),repeat=k):
        if all(c==0 for c in coeffs): continue
        w=[F.ZERO]*n
        for ci,row in zip(coeffs,R):
            if ci==0: continue
            for j in range(n):
                if row[j]!=F.ZERO:
                    w[j]=F.add(w[j],F.mul(ci,row[j]))
        wt=sum(1 for x in w if x!=F.ZERO)
        if 0<wt<best: best=wt
        if cap and best<=cap: pass
    return best

if __name__=="__main__":
    F4=GFTable(2,[1,1,1])
    print("F4 q=",F4.q,"one=",F4.ONE)
    print("2*2=",F4.mul(2,2),"(expect 3)  2*3=",F4.mul(2,3),"(expect 1)  3*3=",F4.mul(3,3),"(expect 2)")
    print("frob(2,1)=",F4.frob(2,1),"(expect 3)")
    G=[[1,0,1,0],[0,1,0,1]]
    print("rank=",rank(F4,G),"dual=",dual_basis(F4,G),"dmin=",min_distance(F4,G))
