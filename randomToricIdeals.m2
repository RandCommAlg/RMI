-- back in business Jan 2018:
-- hot to get toric ideal I_A from a matrix A: 
loadPackage"FourTiTwo"
help toricMarkov
A= transpose matrix flatten apply(allLaurentMonomials(3,2) , m->exponents m)
n=numcols A
R=ZZ/101[x_1..x_n]
toricMarkov (A,R)
--dim oo
--degree oo


-- SP:
-- here is an example of how to get all Laurent monomials 
-- with a given 1-norm. We start with example and move to 
-- a function that works on n variables and 1-norm D. 

{*

R=QQ[a,b,c,x,y,z]
I=ideal(a*x-1,b*y-1,c*z-1)
 F=R/I

-- example:
 basis(2,F)

L = QQ[a,b,c, MonomialOrder=>Lex,Inverses=>true]
  M= matrix{{a,b,c, a^(-1),b^(-1),c^(-1)}}
phi=map(L, F ,M ) 

-- example : 
B  =  flatten entries basis(2,F)
apply(B,b->phi b)

-- the above example gives you all monomials with 1-norm 2. 
*}

-- let's create a funciton for the above: 
allLaurentMonomials = method();
allLaurentMonomials(ZZ,ZZ) := (n,D) -> (
-- input n 
-- input D
R:=ZZ/101[x_1..x_n,a_1..a_n];
I:=ideal apply(toList(1..n),i-> a_i*x_i-1);
F:=R/I;
L:=QQ[x_1..x_n, MonomialOrder=>Lex,Inverses=>true];
phi := map( L , F ,     matrix{join(toList(x_1..x_n), apply(toList(1..n),i->x_i^(-1)) ) } ); 
B  :=  flatten flatten apply(toList(1..D),d->entries basis(d,F));
apply(B,b->phi b)
)
allLaurentMonomials(3,2) 

allLaurentMonomials(ZZ,ZZ,ZZ) := (n,L,U) -> (
-- input n
-- input L<0
-- input U>0
R:=ZZ/101[x_1..x_n,a_1..a_n, Degrees=>join(toList(n:{1,0}), toList(n:{0,-1}))];
I:=ideal apply(toList(1..n), i->a_i*x_i-1);
F:=R/I;
K:=QQ[x_1..x_n, MonomialOrder=>Lex,Inverses=>true];
phi:= map(K,F, matrix{join(toList(x_1..x_n), apply(toList(1..n),i->x_i^(-1)))});
B:= delete(sub(1,F), flatten flatten flatten apply(toList(0..U), i->apply(toList(L..0),j->entries basis({i,j},F))));
apply(B, b->phi b)
)
allLaurentMonomials(3,-2,2) 

randomMonomialSet = method();
randomMonomialSet(ZZ,ZZ,ZZ) := (n,D,M) -> (
-- fixed M model with L1 norm monomial generationg model
allMonomials = allLaurentMonomials(n,D);
B = take(random(allMonomials),M)
)
randomMonomialSet(3,2,8)

randomMonomialSet(ZZ,ZZ,RR) := (n,D,p) -> (
-- ER model with L1 norm monomial generating model
allMonomials = allLaurentMonomials(n,D);
B = select(allMonomials, m->random(0.0,1.0)<=p)
)
randomMonomialSet(3,2,.2)

randomMonomialSet(ZZ,ZZ,ZZ,ZZ) := (n,L,U,M) -> (
-- fixed M model with positive degree sum/negative degree sum monomial generating model
allMonomials = allLaurentMonomials(n,L,U);
B = take(random(allMonomials),M)
)
randomMonomialSet(3,-2,1,5)

randomMonomialSet(ZZ,ZZ,ZZ,RR) := (n,L,U,p) -> (
-- ER model with positive degree sum/negative degree sum monomial generating model
allMonomials = allLaurentMonomials(n,L,U);
B= select(allMonomials, m->random(0.0,1.0)<=p)
)
randomMonomialSet(3,-1,2,.2)

randomMonomialSet(ZZ,ZZ,List) := (n,D,pOrM) -> (
-- start of graded model
if all(pOrM,q->instance(q,ZZ)) then (
        allMonomials = sort values partition(m-> first degree m, allLaurentMonomials(n,D));
        B = flatten apply(toList(1..(2*D+1)), d->take(random(allMonomials_(d-1)), pOrM_(d-1)))
    )
else if all(pOrM,q->instance(q,RR)) then (
        allMonomials = sort values partition(m-> first degree m, allLaurentMonomials(n,D));
        B = flatten apply(toList(1..(2*D+1)), d->select(allMonomials_(d-1),m->random(0.0,1.0)<=pOrM_(d-1)))
    )
)
randomMonomialSet(3,2,{2,1,0,1,2})
randomMonomialSet(3,2,{.5,.2,.0,.2,.5})

randomMonomialSet(ZZ,ZZ,ZZ,List) := (n,L,U,pOrM) -> (
-- start of graded model
if all(pOrM,q->instance(q,ZZ)) then (
        allMonomials = sort values partition(m-> first degree m, allLaurentMonomials(n,L,U));
        B = flatten apply(toList(1..(U-L+1)), d->take(random(allMonomials_(d-1)), pOrM_(d-1)))
    )
else if all(pOrM,q->instance(q,RR)) then (
        allMonomials = sort values partition(m-> first degree m, allLaurentMonomials(n,L,U));
        B = flatten apply(toList(1..(U-L+1)), d->select(allMonomials_(d-1),m->random(0.0,1.0)<=pOrM_(d-1)))
    )
)
randomMonomialSet(3,-1,2,{1,0,1,2})
randomMonomialSet(3,-1,2,{.2,.0,.2,.5})





