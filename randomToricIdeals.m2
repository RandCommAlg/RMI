
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
allLaurentMonomials=method()

allLaurentMonomials = (n,D) -> (
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

allLaurentMonomials = (n,L,U) -> (
R:=ZZ/101[x_1..x_n,a_1..a_n, Degrees=>join(toList(n:{1,0}), toList(n:{0,-1}))];
I:=ideal apply(toList(1..n), i->a_i*x_i-1);
F:=R/I;
K:=QQ[x_1..x_n, MonomialOrder=>Lex,Inverses=>true];
phi:= map(K,F, matrix{join(toList(x_1..x_n), apply(toList(1..n),i->x_i^(-1)))});
B:= delete(sub(1,F), flatten flatten flatten apply(toList(0..U), i->apply(toList(L..0),j->entries basis({i,j},F))));
apply(B, b->phi b)
)

allLaurentMonomials(3,-2,2) 

randomLMonomialSet = (n,D,M) -> (
-- fixed M model
allMonomials = allLaurentMonomials(n,D);
B = take(random(allMonomials),M)
)

randomLMonomialSet = (n,D,p) -> (
-- ER model
allMonomials = allLaurentMonomials(n,D);
B = select(allMonomials, m->random(0.0,1.0)<=p)
)


