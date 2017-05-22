--**************************--
-- -*- coding: utf-8 -*-
newPackage(
	"RandomMonomialIdeals",
    	Version => "1.0", 
    	Date => "May 5, 2017",
    	Authors => {
	    {
		Name => "Sonja Petrovic", 
		Email => "sonja.petrovic@iit.edu", 
		HomePage => "http://math.iit.edu/~spetrov1/"
	    },
	    {
		Name => "Despina Stasi", 
		Email => "stasdes@iit.edu", 
		HomePage => "http://math.iit.edu/~stasdes/"
	    },	
	    {
		Name => "Dane Wilburne", 
		Email => "dwilburn@hawk.iit.edu", 
		HomePage => "http://mypages.iit.edu/~dwilburn/"
	    },	
	    {
		Name => "Tanner Zielinski", 
		Email => "tzielin1@hawk.iit.edu", 
		HomePage => "https://www.linkedin.com/in/tannerzielinski/"
	    },	
	    {
		Name => "Daniel Kosmas", 
		Email => "dkosmas@hawk.iit.edu", 
		HomePage => "https://www.linkedin.com/in/daniel-kosmas-03160988/"
	    },
	    {
		Name => "Parker Joncus", 
		Email => "pjoncus@hawk.iit.edu", 
		HomePage => ""
	    },
	    {
		Name => "Richard Osborn", 
		Email => "rosborn@hawk.iit.edu", 
		HomePage => ""
	    },
	    {
	    	Name => "Monica Yun", 
	    	Email => "myun1@hawk.iit.edu", 
	    	HomePage => ""
	    },
	    {
	    	Name => "Genevieve Hummel", 
	    	Email => "ghummel1@hawk.iit.edu", 
	    	HomePage => ""
	    }
          -- {Name=> "Contributing authors and collaborators: add any acknowledgements here", 
	  -- Email=> "",
	  -- HomePage=>""}      
	},
    	Headline => "A package for generating Erdos-Renyi-type random monomial ideals",
    	DebuggingMode => false,
	Reload => true 
    	)

export {
    "randomGeneratingSets"
    }

--***************************************--
--  Exported methods 	     	     	 --
--***************************************--

randomGeneratingSets = method(TypicalValue => List)
randomGeneratingSets (ZZ,ZZ,RR,ZZ) := List =>  (n,D,p,N) -> (
    x :=symbol x;
    R := QQ[x_1..x_n];

    allMonomials := flatten flatten apply(toList(1..D),d->entries basis(d,R));
    -- this generates a list of all possible monomials of degree <= D in n variables
    -- go through list allMonomials, and for each monomial m in the list, select a number in Unif(0,1);
    -- if that number <= p, then include the monomial m in a generating set B
    -- since iid~Unif(0,1), this is same as keeping each possible monomial w/ probability p.
    -- we need a sample of size N of sets of monomials like these, so we repeat this process N times:
    B := apply(N,i-> select(allMonomials, m-> random(0.0,1.0)<=p) );
    -- we need 0 ideals stored in the appropriate ring; hence do this: 
    apply(#B,i-> if B_i==={} then B=replace(i,{0_R},B));
    return(B)
)
randomGeneratingSets (ZZ,ZZ,ZZ,ZZ) := List =>  (n,D,M,N) -> (
    x :=symbol x;
    R := QQ[x_1..x_n];
    --this generates a list of all possible monomials of degree <=D in n variables
    --randomizes the list of all monomials and selects the first M as the a generating set
    --this is repeated N times to get a sample size of N sets of monomials 
    allMonomials := toList(flatten flatten apply(toList(1..D),d->entries basis(d,R)));
    B :=  apply(N,i-> take(random(allMonomials), {0,M-1}) );
    return(B)
)

randomGeneratingSets (ZZ,ZZ,List,ZZ) := List =>  (n,D,p,N) -> (
    x :=symbol x;
    R := QQ[x_1..x_n];
    B := apply(N,i-> flatten apply(toList(1..D),d-> select(flatten entries basis(d,R),m-> random(0.0,1.0)<=p_(d-1))));
    apply(#B,i-> if B_i==={} then B=replace(i,{0_R},B));
    return(B)
)

--**********************************--
--  Internal methods	    	    --
--**********************************--



--******************************************--
-- DOCUMENTATION     	       	    	    -- 
--******************************************--
beginDocumentation()

doc ///
 Key
  RandomMonomialIdeals
 Headline
  A package for generating Erdos-Renyi-type random monomial ideals
 Description
  Text
   {\em RandomMonomialIdeals} is a  package that... 
  -- Caveat
  -- Still trying to figure this out. [REMOVE ME]
///

doc ///
 Key
  randomGeneratingSets
  (randomGeneratingSets,ZZ,ZZ,RR,ZZ)
  (randomGeneratingSets,ZZ,ZZ,ZZ,ZZ)
  (randomGeneratingSets,ZZ,ZZ,List,ZZ)
 Headline
  randomly generates lists of monomials, up to a given degree
 Usage
  randomGeneratingSets (ZZ,ZZ,RR,ZZ)
  randomGeneratingSets(ZZ,ZZ,ZZ,ZZ)
  randomGeneratingSets (ZZ,ZZ,List,ZZ)
 Inputs
  n: ZZ
    number of variables
  D: ZZ
    maximum degree
  p: RR
     or @ofClass List@
     , probability to select a monomial
  M: ZZ
     number of monomials in each generating set
  N: ZZ
    number of sets generated
 Outputs
  B: List
   random generating sets of monomials
 Description
  Text
   randomGeneratingSets creates $N$ random sets of monomials of degree $d$, $1\leq d\leq D$, in $n$ variables. 
   If $p$ is a real number, it generates each of these sets according to the Erdos-Renyi-type model: 
   from the list of all monomials of degree $1,\dots,D$ in $n$ variables, it selects each one, independently, with probability $p$. 
  Example
   randomGeneratingSets(2,3,0.2,10)
   randomGeneratingSets(3,2,0.6,4)
  Text
   Note that this model does not generate the monomial $1$: 
  Example
   randomGeneratingSets(3,2,1.0,1)
  Text 
   If $M$ is an integer, then randomGeneratingSets creates $N$ random sets of monomials of size $M$:
   randomly select $M$ monomials from the list of all monomials of degree $1,\dots,D$ in $n$ variables.
  Example
   randomGeneratingSets(2,3,3,1)
  Text
   Note that the degree 1 monomials were not generated, and each set has $M$ monomials.
  Text
   If $M$ is bigger than the total number of monomials in $n$ variables of degree at most $D$, then the method will simply return all those monomials (and not $M$ of them).
  Example
   randomGeneratingSets(2,2,10,1)
  Text
   returns 5 monomials in a generating set, and not 10, since there do not exist 10 to choose from.
  Text 
   If $p=p_1,\dots,p_D$ is a list of real numbers of length $D$, then randomGeneratingSets generates the sets utilizing the graded Erdos-Renyi-type model:
   select each monomial of degree $1\le d\le D$, independently, with probability $p_d$.
  Example
   randomGeneratingSets(2,3,{0.0,1.0,1.0},1)
  Text
   Note that the degree 1 monomials were not generated.
///



--******************************************--
-- TESTS     	     	       	    	    -- 
--******************************************--

--************************--
--  randomGeneratingSets  --
--************************--

TEST ///
    -- Check there are N samples
    N=10;
    n=3; D=2; p=0.5;
    assert (N==#randomGeneratingSets(n,D,p,N))
    N=13;
    n=5; D=3; p={0.5,0.25,0.3};
    assert (N==#randomGeneratingSets(n,D,p,N))
    N=10;
    n=3; D=2; M=10;
    assert (N==#randomGeneratingSets(n,D,M,N))
///

TEST ///
    -- Check no terms are chosen for a probability of 0
    assert (0==(randomGeneratingSets(5,5,0.0,1))#0#0)
    assert (0==(randomGeneratingSets(5,4,toList(4:0.0),1))#0#0)
///

TEST ///
    -- Check all possible values are outputted with a probability of 1
    D=3; n=4;
    assert (product(toList((D+1)..D+n))/n!-1==#(randomGeneratingSets(n,D,1.0,1))#0)
    assert (product(toList((D+1)..D+n))/n!-1==#(randomGeneratingSets(n,D,{1.0,1.0,1.0},1))#0)
    D=2; n=6;
    assert (product(toList((D+1)..D+n))/n!-1==#(randomGeneratingSets(n,D,1.0,1))#0)
    assert (product(toList((D+1)..D+n))/n!-1==#(randomGeneratingSets(n,D,{1.0,1.0},1))#0)
///

TEST ///
    -- Check every monomial is generated
    L=(randomGeneratingSets(2,3,1.0,1))#0
    R=ring(L#0)
    assert(set L===set {R_0,R_1,R_0^2,R_0*R_1,R_1^2,R_0^3,R_0^2*R_1,R_0*R_1^2,R_1^3})
    L=(randomGeneratingSets(3,3,{0.0,1.0,0.0},1))#0
    R=ring(L#0)
    assert(set L===set {R_0^2,R_0*R_1,R_1^2,R_0*R_2,R_1*R_2,R_2^2})
    L=(randomGeneratingSets(2,3,9,1))#0
    R=ring(L#0)
    assert(set L===set {R_0,R_1,R_0^2,R_0*R_1,R_1^2,R_0^3,R_0^2*R_1,R_0*R_1^2,R_1^3})
///

TEST ///
    -- Check max degree of monomial less than or equal to D
    n=10; D=5;
    assert(D==max(apply((randomGeneratingSets(n,D,1.0,1))#0,m->first degree m)))
    assert(D==max(apply((randomGeneratingSets(n,D,toList(D:1.0),1))#0,m->first degree m)))
    n=4; D=7;
    assert(D==max(apply((randomGeneratingSets(n,D,1.0,1))#0,m->first degree m)))
    assert(D==max(apply((randomGeneratingSets(n,D,toList(D:1.0),1))#0,m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(D==max(apply((randomGeneratingSets(n,D,M,1))#0,m->first degree m)))
    n=4; D=7;
    assert(D==max(apply((randomGeneratingSets(n,D,1.0,1))#0,m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(D==max(apply((randomGeneratingSets(n,D,M,1))#0,m->first degree m)))
///

TEST ///
    -- Check min degree of monomial greater than or equal to 1
    n=8; D=6;
    assert(1==min(apply((randomGeneratingSets(n,D,1.0,1))#0,m->first degree m)))
    assert(1==min(apply((randomGeneratingSets(n,D,toList(D:1.0),1))#0,m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(1==min(apply((randomGeneratingSets(n,D,M,1))#0,m->first degree m)))
    n=3; D=5;
    assert(1==min(apply((randomGeneratingSets(n,D,1.0,1))#0,m->first degree m)))
    assert(1==min(apply((randomGeneratingSets(n,D,toList(D:1.0),1))#0,m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(1==min(apply((randomGeneratingSets(n,D,M,1))#0,m->first degree m)))
///

end

You can write anything you want down here.  I like to keep examples
as I’m developing here.  Clean it up before submitting for
publication.  If you don't want to do that, you can omit the "end"
above.
