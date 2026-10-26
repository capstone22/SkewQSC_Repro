import csv, itertools, json, math, time
from pathlib import Path
from functools import lru_cache
from ffield import GFTable, rref, rank, dual_basis, in_rowspace
from ore import (Ore, xn_minus_1, monic_divisors_right, skew_gen_matrix)
from mindist import min_distance_fast

def poly_div_exact(O,a,b):
    F=O.F
    if not O.right_divides(b,a): return None
    db=O.deg(b); aa=O.trim(a[:]); q=[F.ZERO]*(O.deg(a)-db+1)
    while O.deg(aa)>=db and O.deg(aa)>=0:
        da=O.deg(aa); c=aa[da]; e=da-db; q[e]=c
        term=[F.ZERO]*(e+1); term[e]=c
        aa=O.add(aa,O.neg(O.mul(term,b)))
    return O.trim(q)

def cyclic_opt(F,n,k1,k2):
    Oid=Ore(F,0); target=xn_minus_1(F,n)
    d1=n-k1; d2=n-k2
    best=0
    for g1 in monic_divisors_right(Oid,target,d1):
        for g2 in monic_divisors_right(Oid,g1,d2):
            q=poly_div_exact(Oid,g1,g2)
            if q is None or Oid.deg(q)!=k2-k1 or q[0]==F.ZERO: continue
            o=Oid.order_comm(q)
            if o and o>best: best=o
    return best

def gain_search(F,tp,n):
    O=Ore(F,tp); m=O.m
    if n%m!=0: return []
    target=xn_minus_1(F,n)
    cyc_cache={}
    hits=[]
    for dg1 in range(1,(n+1)//2):
        k1=n-dg1
        if 2*k1-n<=0: continue
        for g1 in monic_divisors_right(O,target,dg1):
            G=skew_gen_matrix(O,g1,n); Rg,pivg=rref(F,G)
            if len(pivg)!=k1: continue
            Db=dual_basis(F,G)
            if not all(in_rowspace(F,Rg,pivg,v) for v in Db): continue
            for dg2 in range(0,dg1):
                k2=n-dg2; ds=dg1-dg2
                g2list=monic_divisors_right(O,g1,dg2) if dg2>0 else [[F.ONE]]
                for g2 in g2list:
                    fs=poly_div_exact(O,g1,g2)
                    if fs is None or O.deg(fs)!=ds or fs[0]==F.ZERO: continue
                    if O.mul(fs,g2)!=O.trim(g1[:]): continue
                    ind=O.skew_index(fs)
                    if ind is None: continue
                    key=(n,k1,k2)
                    if key not in cyc_cache: cyc_cache[key]=cyclic_opt(F,n,k1,k2)
                    cyc=cyc_cache[key]
                    if cyc>0 and ind>cyc:
                        hits.append(dict(n=n,k1=k1,k2=k2,ds=ds,ind=ind,cyc=cyc,
                                         gain=ind/cyc,kL=2*k1-n,g1=g1[:],g2=g2[:],fs=fs[:]))
    return hits

if __name__=="__main__":
    fields={'F4':(GFTable(2,[1,1,1]),1,16),
            'F8':(GFTable(2,[1,1,0,1]),1,12),
            'F9':(GFTable(3,[1,0,1]),1,10)}
    allhits=[]
    for name,(F,tp,Nmax) in fields.items():
        m=Ore(F,tp).m; t0=time.time()
        for n in range(m,Nmax+1):
            if n%m: continue
            h=gain_search(F,tp,n)
            for r in h:
                r['field']=name
                allhits.append(r)
        print(f"{name}: done in {time.time()-t0:.1f}s, hits so far={len(allhits)}")
    print(f"\nTOTAL dual-containing CSS hits with k_L>0 and strict gain: {len(allhits)}")
    for name,(F,tp,_) in fields.items():
        O=Ore(F,tp)
        for r in [x for x in allhits if x['field']==name]:
            G1=skew_gen_matrix(O,r['g1'],r['n']); G2=skew_gen_matrix(O,r['g2'],r['n'])
            d1=min_distance_fast(F,G1); d2=min_distance_fast(F,G2)
            r['d1']=d1; r['d2']=d2
    import pickle
    with Path(__file__).with_name('hits.pkl').open('wb') as stream:
        pickle.dump(allhits, stream)
    columns = ('field', 'n', 'k1', 'k2', 'ds', 'ind', 'cyc', 'gain',
               'kL', 'd1', 'd2', 'g1', 'g2', 'fs')
    with Path(__file__).with_name('positive_hits.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for hit in allhits:
            writer.writerow({key: json.dumps(hit[key]) if key in ('g1', 'g2', 'fs')
                             else hit[key] for key in columns})
    print(f"Wrote {len(allhits)} search hits to hits.pkl and positive_hits.csv")
