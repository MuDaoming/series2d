"""Euclidean reconstruction of r*D, learning missing poles as well as powers."""
from reference_space import RatFunc

def make_preconditioned_fit(L,p,center,holdout=8):
 cache={}
 def fit(D,times,values):
  xs=tuple((t-center)%p for t in times);n=len(xs)-holdout
  if xs not in cache:
   M=L.one
   for x in xs[:n]:M*=L.poly([-x%p,1])
   bases=[]
   for x in xs[:n]:
    b=M//L.poly([-x%p,1]);bases.append(b*pow(int(b(x)),-1,p))
   cache[xs]=(M,bases)
  M,bases=cache[xs];ys=[v*int(D(x))%p for x,v in zip(xs,values)]
  a=L.poly([0])
  for b,v in zip(bases,ys):a+=b*v
  r0,r1,s0,s1=M,a,L.poly([0]),L.one
  while r1:
   if r1.degree()+s1.degree()<n:
    if all(int(s1(xs[k])) and int(r1(xs[k]))==ys[k]*int(s1(xs[k]))%p for k in range(n,len(xs))):
     f=RatFunc.make(list(r1),list(s1*D),p)
     return f
   q,rem=divmod(r0,r1);r0,r1=r1,rem;s0,s1=s1,s0-q*s1
  return None
 return fit
