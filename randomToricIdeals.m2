
-- SP:
-- here is an example of how to get all Laurent monomials 
-- with a given 1-norm. We start with example and move to 
-- a function that works on n variables and 1-norm D. 

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


-- let's create a funciton for the above: 
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
