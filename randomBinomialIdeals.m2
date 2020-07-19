------------------------------------------
----- this part: Jun 2020 code for generating binomial ideals as differences of randomly generated monomials: 

--- Notes to self.
-- I guess one can try something even more straightfoward so we never create the big list
-- of monomials (do we care?): 
-- mutliply randomly many random variables to get a monomial and repeat to get the 2nd one.
-- is that silly? it was just easier for me to code as what i did below. 


---------------------------------------------------------------------------------------------
--      randomBinomials
---------------------------------------------------------------------------------------------
-- make a list "mons" of all mons of degree d in ring R; 
-- then select two random indices, indPlus and indMinus,
-- and make a binomial out of mons_indPlus - mons_indMinus,
-- and add k such binomials to a list. 
-- 
-- if Homogeneous=>true then: 
-- OUTPUT: list of k homogeneous binomials (of format M1-M2) of degree d in ring R.
-- if Homogeneous=>false then: 
-- OUTPUT: list of k (non-homogeneous) binomials (of format M1-M2) of degree at most d in ring R.
---------------------------------------------------------------------------------------------
randomBinomials = method(TypicalValue => List, Options=>{Homogeneous=>false})
randomBinomials(PolynomialRing,ZZ,ZZ) := List => o -> (R,maxDegree,k) -> ( 
    if o.Homogeneous then mons = flatten entries basis(maxDegree,R) else (
	mons = {};
	scan(0..maxDegree, d-> mons = append(mons, flatten entries basis(d,R)));
	mons = flatten mons
	);
    mons = drop(mons,1); -- no 0.  TO DO: enable this later; it may produce monomials not true binomials.
    binomials = {};
    scan(k, i-> (
	    indPlus = random(0,#mons-1);
	    indMinus = random(0,#mons-1);
	    binomials = append(binomials, mons_indPlus-mons_indMinus);
	    )
    	);
    flatten binomials
    )
-- but the list may include zeros. 
---------------------------------------------------------------------------------------------

-- see code at end of file for examples on how this runs. 
-- here's a short example of usage:
{*
    randomBinomials(QQ[x,y,z],3,4)
    randomBinomials(QQ[x,y,z],3,4,Homogeneous=>true)
    *}

writeData = method()
-- outline of method taken from RandomMonomialIdeals.m2 writeSample :) 
writeData (List, String, String) := (s, parameters, filename) -> (
    if fileExists filename then (
	stderr << "error: file with this name already exists. Rename the old file and try again." << endl;
 	--stderr << "warning: overwrting file." << endl;
	--    removeFile filename;
	) 
    else filename << parameters << endl << s << endl << close;
)
writeData (Array, String, String) := (s, parameters, filename) -> (
    if fileExists filename then (
	stderr << "error: file with this name already exists. Rename the old file and try again." << endl;
 	--stderr << "warning: overwrting file." << endl;
	--    removeFile filename;
	) 
    else filename << parameters << endl << s << endl << close;
)




-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------
----- below: Jan 2018 code for generating I_A from a random A, where A is a mtx of exponents of randomly generated Laurent monomials.


-- let's create a function for the above: 
allLaurentMonomials = method(TypicalValue=>List);
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

randomMonomialSet = method();
randomMonomialSet(ZZ,ZZ,ZZ) := (n,D,M) -> (
-- fixed M model with L1 norm monomial generationg model
allMonomials = allLaurentMonomials(n,D);
B = take(random(allMonomials),M)
)

randomMonomialSet(ZZ,ZZ,RR) := (n,D,p) -> (
-- ER model with L1 norm monomial generating model
allMonomials = allLaurentMonomials(n,D);
B = select(allMonomials, m->random(0.0,1.0)<=p)
)

randomMonomialSet(ZZ,ZZ,ZZ,ZZ) := (n,L,U,M) -> (
-- fixed M model with positive degree sum/negative degree sum monomial generating model
allMonomials = allLaurentMonomials(n,L,U);
B = take(random(allMonomials),M)
)

randomMonomialSet(ZZ,ZZ,ZZ,RR) := (n,L,U,p) -> (
-- ER model with positive degree sum/negative degree sum monomial generating model
allMonomials = allLaurentMonomials(n,L,U);
B= select(allMonomials, m->random(0.0,1.0)<=p)
)

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




--##########################################################################--
end   -- terminate reading ...
--##########################################################################--


--***************************************************************************--
-- Some examples and how to run, etc. 
--***************************************************************************--

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

{* -- this is a multi-line comment --

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


allLaurentMonomials(3,2) 
allLaurentMonomials(3,-2,2) 
-- fixed M model with L1 norm monomial generationg model: 
randomMonomialSet(3,2,8)
-- ER model with L1 norm monomial generating model:
randomMonomialSet(3,2,.2)
-- fixed M model with positive degree sum/negative degree sum monomial generating model:
randomMonomialSet(3,-2,1,5)
-- ER model with positive degree sum/negative degree sum monomial generating model:
randomMonomialSet(3,-1,2,.2)


-- graded? :  [appears to not be done yet??] 
randomMonomialSet(3,2,{2,1,0,1,2})
randomMonomialSet(3,2,{.5,.2,.0,.2,.5})
randomMonomialSet(3,-1,2,{1,0,1,2})
randomMonomialSet(3,-1,2,{.2,.0,.2,.5})



---------------------------------------------------------------------------------------------
-- EXAMPLE for trying things out for randomBinomials:
-- Just trying an example before I wrote the method: 
---------------------------------------------------------------------------------------------

R=QQ[a,b,c,x,y,z]

randomBinomials(R,2,3,Homogeneous=>true)
-- now we can maybe use things like bettiStats or whatever to see how many gens we got. 
-- but not after we hav some more general code. Right now we're getting homogeneous binomials of 
-- fixed degree. 

-- want to change degree on the fly? use Homogeneous=>false. 



--degree check:
betti ideal  binomials
--these are *not* homogeneous. it now makes sense to run a command like degreeStats (from RMI) :) 
loadPackage"RandomMonomialIdeals"
----degStats (binomials,ShowTally => true)--this breaks so need to make each bin into its ideal (stupid):
--bin = apply(binomials,b->ideal(b)); 

--degStats (bin,ShowTally => true)
s=new Sample from {SampleSize => k, Data=>binomials} --,ModelName=>"foo",Parameters=>k}
binomialDegrees = apply(s.Data, degree@@ideal)



---------------------
---->>> I DON'T THINK I NEED THIS OPTION BELOW: to get it, simply run the Homogeneous=>true
-- in the method above as many times as you want for each desired degree. :) 
---------------------
-- Note: maybe we want homogeneous binomials in degrees from 0 to maxDegree? Let's try that too: 
maxDegree = 4;
mons = {};
scan(0..maxDegree, d-> mons = append(mons, flatten entries basis(d,R)));
-- don't flatten! this keeps the distinct degrees in distinct sublists :) 
k = 3; -- how many binomials do you want? 
binomials = {};
scan(maxDegree, i-> (
	currentDegree = random(1,maxDegree);
	indPlus = random(0,#mons_currentDegree-1);
	indMinus = random(0,#mons_currentDegree-1);
	binomials = append(binomials, mons_currentDegree_indPlus-mons_currentDegree_indMinus)
	)
    )
binomials = flatten binomials -- <<<<<, THIS IS WHAT WE WANT TO SAVE. 
--degree check:
betti ideal  binomials -- i got 1 quadric nad 3 quartics. 
s=new Sample from {SampleSize => k, Data=>binomials} --,ModelName=>"foo",Parameters=>k}
peek s
-- but this is OK:
binomialDegrees = apply(s.Data, degree@@ideal)