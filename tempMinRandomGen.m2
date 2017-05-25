randomMinGeneratingSet = method()
randomMinGeneratingSet(ZZ, ZZ, RR) := (n,D,p) -> (
    return randomMinGeneratingSet(n,D,D:p)
 )

randomMinGeneratingSet(ZZ, ZZ, Sequence) := (n,D,p) -> (
 x:=symbol x;
 R:=QQ[x_1..x_n];
 B:={};
 currentRing=R;
 apply(D, d-> (
  chosen:=select(flatten entries basis(d+1,currentRing), m->random(0.0,1.0)<=p_d);
  B=flatten append(B, chosen/(i->sub(i, R)));
  currentRing = currentRing/promote(ideal(chosen),currentRing)
  )
);
 return(B)
 )
TEST ///
-- Check no terms are chosen for a probability of 0
n=4;
assert ({}=randomMinGeneratingSet(n,5,0.0))
assert ({}=randomMinGeneratingSet(n,4,n:0.0))
///
TEST ///
-- Check all possible monomials of degree d are outputted with probability 1.0
n=4;
D=5;
assert (# flatten entries basis (1, QQ[x_1..x_n])==#randomMinGeneratingSet(n,D,1.0))
assert (# flatten entries basis (2, QQ[x_1..x_n])==#randomMinGeneratingSet(n,D,(0.0,1.0,1.0,1.0,1.0))
assert (# flatten entries basis (3, QQ[x_1..x_n])==#randomMinGeneratingSet(n,D,(0.0,0.0,1.0,1.0,1.0))
assert (# flatten entries basis (4, QQ[x_1..x_n])==#randomMinGeneratingSet(n,D,(0.0,0.0,0.0,1.0,1.0))
assert (# flatten entries basis (5, QQ[x_1..x_n])==#randomMinGeneratingSet(n,D,(0.0,0.0,0.0,0.0,1.0))
///
TEST ///
-- Check all monomials of degree d are generated with probability 1.0
L=randomMinGeneratingSet(3,3,1.0);
R=ring(L#0)
assert(set L===set {R_0, R_1, R_2})
L=randomMinGeneratingSet(3,3,(0.0,1.0,1.0))
assert(set L===set {R_0^2,R_0*R_1,R_1^2,R_0*R_2,R_1*R_2,R_2^2})
L=randomMinGeneratingSet(3,3,(0.0,0.0,1.0))
assert(set L===set {R_0^3,R_0^2*R_1,R_0^2*R_2,R_0*R_1*R_2,R_1^3,R_0*R_1^2,R_1^2*R_2,R_0*R_2^2,R_1*R_2^2,R_2^3})
///
TEST ///
-- Check max degree of monomial less than or equal to D
n=10;
D=5;
assert(D==max(apply((randomGeneratingSets(n,D,(0.0,0.0,0.0,0.0,1.0))#0,m->first degree m)))
///
TEST ///
-- Check min degree of monomial greater than or equal to 1
n=10;
D=5;
assert(D==min(apply((randomGeneratingSets(n,D,1.0)#0,m->first degree m)))
///
