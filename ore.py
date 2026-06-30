import itertools
from ffield import GFTable, rref, rank, dual_basis, in_rowspace, min_distance

class Ore:
    def __init__(self,F,theta_power):
        self.F=F; self.tp=theta_power
        r=F.r; g=__import__('math').gcd(theta_power%r if r else 0, r) if r else 1
        self.m = r//g if r>0 else 1
        self.theta_maps=[]
        for k in range(self.m):
            self.theta_maps.append([F.frob(a, (theta_power*k)%(r if r else 1)) for a in range(F.q)])
    def th(self,a,k=1):
        return self.theta_maps[k%self.m][a]
    def trim(self,p):
        F=self.F
        while len(p)>1 and p[-1]==F.ZERO: p.pop()
        return p
    def deg(self,p):
        p=self.trim(p[:])
        if len(p)==1 and p[0]==self.F.ZERO: return -1
        return len(p)-1
    def add(self,p,q):
        F=self.F; n=max(len(p),len(q)); r=[F.ZERO]*n
        for i in range(len(p)): r[i]=F.add(r[i],p[i])
        for i in range(len(q)): r[i]=F.add(r[i],q[i])
        return self.trim(r)
    def neg(self,p): return [self.F.neg(c) for c in p]
    def mul(self,p,q):
        F=self.F
        if self.deg(p)<0 or self.deg(q)<0: return [F.ZERO]
        r=[F.ZERO]*(len(p)+len(q)-1)
        for i,pi in enumerate(p):
            if pi==F.ZERO: continue
            for j,qj in enumerate(q):
                if qj==F.ZERO: continue
                r[i+j]=F.add(r[i+j], F.mul(pi, self.th(qj,i)))
        return self.trim(r)
    def right_rem(self,a,f):
        F=self.F; a=self.trim(a[:]); f=self.trim(f[:]); df=self.deg(f)
        assert f[df]==F.ONE, "f must be monic"
        while self.deg(a)>=df and self.deg(a)>=0:
            da=self.deg(a); c=a[da]; e=da-df
            term=[F.ZERO]*(e+1); term[e]=c
            a=self.add(a, self.neg(self.mul(term,f)))
        return self.trim(a)
    def right_divides(self,f,a):
        return self.deg(self.right_rem(a,f))<0
    def x_times(self,r):
        F=self.F
        return self.trim([F.ZERO]+[self.th(c,1) for c in r])
    def skew_index(self,fs):
        F=self.F; one=[F.ONE]; cur=one[:]
        limit=F.q**self.deg(fs)+2
        for step in range(1,limit):
            cur=self.right_rem(self.x_times(cur),fs)
            if cur==one: return step
        return None
    def order_comm(self,f):
        F=self.F; one=[F.ONE]
        if self.deg(f)==0: return None
        cur=self.right_rem([F.ZERO,F.ONE],f); L=1
        x=[F.ZERO,F.ONE]
        while cur!=one:
            cur=self.right_rem(self.mul(cur,x),f); L+=1
            if L>F.q**self.deg(f): return None
        return L

def xn_minus_1(F,n):
    p=[F.ZERO]*(n+1); p[0]=F.neg(F.ONE); p[n]=F.ONE
    return p

def monic_divisors_right(O,target,s):
    F=O.F; out=[]
    for coeffs in itertools.product(range(F.q),repeat=s):
        if coeffs and coeffs[0]==F.ZERO: continue
        f=list(coeffs)+[F.ONE]
        if O.right_divides(f,target): out.append(f)
    return out

def left_mult_x_vec(O,vec,n):
    F=O.F; out=[F.ZERO]*n
    for i,v in enumerate(vec):
        if v==F.ZERO: continue
        j=(i+1)%n
        out[j]=F.add(out[j], O.th(v,1))
    return out

def skew_gen_matrix(O,g1,n):
    F=O.F; dg=O.deg(g1); k=n-dg
    v=[F.ZERO]*n
    for i,c in enumerate(g1): v[i]=c
    rows=[v[:]]
    for _ in range(k-1):
        v=left_mult_x_vec(O,v,n); rows.append(v[:])
    return rows
