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
    "randomGeneratingSets",
    "randomGeneratingSet",
    "idealsFromGeneratingSets",
    "Coefficients",
    "VariableName",
    "Strategy",
    "IncludeZeroIdeals"
    }

--***************************************--
--  Exported methods 	     	     	 --
--***************************************--

randomGeneratingSets = method(TypicalValue => List, Options => {Coefficients => QQ,
	                                                        VariableName => "x",
								Strategy => "ER"})
randomGeneratingSets (ZZ,ZZ,RR,ZZ) := List => o -> (n,D,p,N) -> (
    if p<0.0 or 1.0<p then error "p expected to be a real number between 0.0 and 1.0";
    randomGeneratingSets(n,D,toList(D:p),N,
	                 Coefficients=>o.Coefficients,
			 VariableName=>o.VariableName,
			 Strategy=>o.Strategy)
)

randomGeneratingSets (ZZ,ZZ,ZZ,ZZ) := List => o -> (n,D,M,N) -> (
    if N<1 then stderr << "warning: N expected to be a positive integer" << endl;
    apply(N,i-> randomGeneratingSet(n,D,M,
	                            Coefficients=>o.Coefficients,
				    VariableName=>o.VariableName,
				    Strategy=>o.Strategy))
)

randomGeneratingSets (ZZ,ZZ,List,ZZ) := List => o -> (n,D,p,N) -> (
    if N<1 then stderr << "warning: N expected to be a positive integer" << endl;
    apply(N,i-> randomGeneratingSet(n,D,p,
	                            Coefficients=>o.Coefficients,
				    VariableName=>o.VariableName,
				    Strategy=>o.Strategy))
)

randomGeneratingSet = method(TypicalValue => List, Options => {Coefficients => QQ,
	                                                       VariableName => "x",
							       Strategy => "ER"})
randomGeneratingSet (ZZ,ZZ,RR) := List => o -> (n,D,p) -> (
    if p<0.0 or 1.0<p then error "p expected to be a real number between 0.0 and 1.0";
    randomGeneratingSet(n,D,toList(D:p),
	                Coefficients=>o.Coefficients,
			VariableName=>o.VariableName,
			Strategy=>o.Strategy)
)

randomGeneratingSet (ZZ,ZZ,ZZ) := List => o -> (n,D,M) -> (
    if M<0 then stderr << "warning: M expected to be a nonnegative integer" << endl;
    if o.Strategy === "Minimal" then error "Minimal not implemented for fixed size ER model";
    x := toSymbol o.VariableName;
    R := o.Coefficients[x_1..x_n];
    allMonomials := flatten flatten apply(toList(1..D),d->entries basis(d,R));
    C := take(random(allMonomials), M);
    if C==={} then {0_R} else C
)

randomGeneratingSet (ZZ,ZZ,List) := List => o -> (n,D,p) -> (
    if n<1 then error "n expected to be a positive integer";
    if #p != D then error "p expected to be a list of length D";
    if not all(p, q->instance(q, ZZ)) and not all(p, q->instance(q,RR)) then error "p must be a list of all integers or all real numbers";
    x := toSymbol o.VariableName;
    R := o.Coefficients[x_1..x_n];
    B := {};
    if all(p,q->instance(q,ZZ)) then (
        if o.Strategy === "Minimal" then error "Minimal not implemented for fixed size ER model";
        B = flatten apply(toList(1..D), d->take(random(flatten entries basis(d,R)), p_(d-1)));
	)
    else if all(p,q->instance(q,RR)) then (
        if any(p,q-> q<0.0 or 1.0<q) then error "p expected to be a list of real numbers between 0.0 and 1.0";
        if o.Strategy === "Minimal" then (
            currentRing := R;
            apply(D, d->(
                chosen := select(flatten entries basis(d+1, currentRing), m->random(0.0,1.0)<=p_d);
                B = flatten append(B, chosen/(i->sub(i, R)));
                currentRing = currentRing/promote(ideal(chosen), currentRing)
            )))
        else
            B = flatten apply(toList(1..D),d-> select(flatten entries basis(d,R),m-> random(0.0,1.0)<=p_(d-1)));
	);
    if B==={} then {0_R} else B
)


--creates a list of monomialIdeal objects from a list of monomial generating sets 
idealsFromGeneratingSets =  method(TypicalValue => List, Options => {IncludeZeroIdeals => false})
-- ^^ change this to by default NOT write to file; and if option " SaveToFile=> true " then do write to file.
-- see branch @25 for this fix. 
idealsFromGeneratingSets (List,RR,ZZ,String) := o -> (B,p,D,basefilename) -> (
	-- ^^ we can decide if we want p,D,basefilename to be optionalinputs that are put together in a sequence 
	-- i.e., do (p,D,baseFileName) as input. 
	-- maybe the filename should be optional and make it "temp" for default. 
    N := # B;
    n := numgens ring ideal B#0; -- ring of the first monomial in the first gen set
    -- see branch @25 for the file writing: 
    --    fileNameExt := concatenate("_for_params_n",toString(n),"_p",toString(p),"_D",toString(D),"_N",toString(N));
    --    filename := concatenate(basefilename,"randomIdeals",fileNameExt,".txt");
    ideals := {};
    for i from 0 to #B-1 do {
	ideals = B / (b-> monomialIdeal b);
	--	filename << toString B#i << endl; 
	};
    --    filename<<close;
    (nonzeroIdeals,numberOfZeroIdeals) := extractNonzeroIdeals(ideals);
    print(concatenate("There are ", toString(#B)," ideals in this sample."));
    print(concatenate("Of those, ", toString numberOfZeroIdeals, " were the zero ideal."));
    if o.IncludeZeroIdeals then return ideals else return (nonzeroIdeals,numberOfZeroIdeals); 
)


--**********************************--
--  Internal methods	    	    --
--**********************************--

toSymbol = (p) -> (
     if instance(p,Symbol)
         then p
     else if instance(p,String)
         then getSymbol p
     else
         error ("expected a string or symbol, but got: ", toString p))

-- Internal method that takes as input list of ideals and splits out the zero ideals, counting them:
    -- input list of ideals 
    -- output a sequence (list of non-zero ideals from the list , the number of zero ideals in the list)
-- (not exported, therefore no need to document) 
extractNonzeroIdeals = ( ideals ) -> (
    nonzeroIdeals := select(ideals,i->i != 0);
    numberOfZeroIdeals := # ideals - # nonzeroIdeals;
    -- numberOfZeroIdeals = # positions(B,b-> b#0==0); -- sinze 0 is only included if the ideal = ideal{}, this is safe too
    return(nonzeroIdeals,numberOfZeroIdeals)
    )
-- we may not need the next one for any of the methods in this file; we'll be able to determine this soon. keep for now.
-- Internal method that takes as input list of generating sets and splits out the zero ideals, counting them:
    -- input list of generating sets
    -- output a sequence (list of non-zero ideals from the list , the number of zero ideals in the list)
-- (not exported, therefore no need to document) 
extractNonzeroIdealsFromGens = ( generatingSets ) -> (
    nonzeroIdeals := select(generatingSets,i-> i#0 != 0_(ring i#0)); --ideal(0)*ring(i));
    numberOfZeroIdeals := # generatingSets - # nonzeroIdeals;
    -- numberOfZeroIdeals = # positions(B,b-> b#0==0); -- sinze 0 is only included if the ideal = ideal{}, this is safe too
    return(nonzeroIdeals,numberOfZeroIdeals)
    )

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
  randomGeneratingSets(ZZ,ZZ,RR,ZZ)
  randomGeneratingSets(ZZ,ZZ,ZZ,ZZ)
  randomGeneratingSets(ZZ,ZZ,List,ZZ)
 Inputs
  n: ZZ
    number of variables
  D: ZZ
    maximum degree
  p: RR
     or @ofClass List@
     , probability to select a monomial
  M: ZZ
     , number of monomials in each generating set
  N: ZZ
    number of sets generated
 Outputs
  B: List
   random generating sets of monomials
 Description
  Text
   randomGeneratingSets creates $N$ random sets of monomials of degree $d$, $1\leq d\leq D$, in $n$ variables. 
   It does so by calling @TO randomGeneratingSet$ $N$ times. 
  SeeAlso
   randomGeneratingSet
///

doc ///
 Key
  randomGeneratingSet
  (randomGeneratingSet,ZZ,ZZ,RR)
  (randomGeneratingSet,ZZ,ZZ,ZZ)
  (randomGeneratingSet,ZZ,ZZ,List)
 Headline
  randomly generates a list of monomials, up to a given degree
 Usage
  randomGeneratingSet(ZZ,ZZ,RR)
  randomGeneratingSet(ZZ,ZZ,ZZ)
  randomGeneratingSet(ZZ,ZZ,List)
 Inputs
  n: ZZ
    number of variables
  D: ZZ
    maximum degree
  p: RR
     , the probability of selecting a monomial, 
     or @ofClass List@
     of real numbers whose i-th entry is the probability of selecing a monomial of degree i
  M: ZZ
     , number of monomials in each set, 
     or @ofClass List@
     of integers whose i-th entry is the number of monomials of degree i in each set
 Outputs
  B: List
   random set of monomials
 Description
  Text
   randomGeneratingSet creates a list of monomials, up to a given degree $d$, $1\leq d\leq D$, in $n$ variables. 
   If $p$ is a real number, it generates the set according to the Erdos-Renyi-type model:
   from the list of all monomials of degree $1,\dots,D$ in $n$ variables, it selects each one, independently, with probability $p$.
  Example
   n=2; D=3; p=0.2;
   randomGeneratingSet(n,D,p)
   randomGeneratingSet(3,2,0.6)
  Text
   Note that this model does not generate the monomial $1$:
  Example
   randomGeneratingSet(3,2,1.0)
  Text
   If $M$ is an integer, then randomGeneratingSet creates a list of monomials of size $M$:
   randomly select $M$ monomials from the list of all monomials of degree $1,\dots,D$ in $n$ variables.
  Example
   n=10; D=5; M=4;
   randomGeneratingSet(n,D,M)
  Text
   Note that it returns a set with $M = 4$ monomials.
  Text
   If $M$ is bigger than the total number of monomials in $n$ variables of degree at most $D$, then the method will simply return all those monomials (and not $M$ of them). For example:
  Example
   randomGeneratingSet(2,2,10)
  Text
   returns 5 monomials in a generating set, and not 10, since there are fewer than 10 monomials to choose from.
  Text
   If $p=p_1,\dots,p_D$ is a list of real numbers of length $D$, then randomGeneratingSet generates the set utilizing the graded Erdos-Renyi-type model:
   select each monomial of degree $1\le d\le D$, independently, with probability $p_d$.
  Example
   p={0.0, 1.0, 1.0};
   randomGeneratingSet(2,3,p)
  Text
   Note that the degree-1 monomials were not generated, since the first probability vector entry is 0.
  Text
   If $M=M_1,\dots,M_D$ is a list of integers of length $D$, then randomGeneratingSet creates a list of monomials, where $M_d$ monomials are of degree $d$.
  Example
   M={2,1,1};
   randomGeneratingSet(2,3,M)
  Text
   Observe that there are two degree-1 monomials, one degree-2 monomial, and one degree-3 monomial.
  SeeAlso
   randomGeneratingSets
///


doc ///
  Key
    Coefficients
    [randomGeneratingSet, Coefficients]
    [randomGeneratingSets, Coefficients]
  Headline
    optional input to choose the coefficient ring of the generated polynomials
  Description
    Text
      Put {\tt Coefficients => r} for a choice of ring r as an argument in
      the function @TO randomGeneratingSet@ or @TO randomGeneratingSets@. 
    Example 
      n=2; D=3; p=0.2;
      randomGeneratingSet(n,D,p)
      ring ideal oo
      randomGeneratingSet(n,D,p,Coefficients=>ZZ/101)
      ring ideal oo
  SeeAlso
    randomGeneratingSet
    randomGeneratingSets
///

doc ///
  Key
    VariableName
    [randomGeneratingSet, VariableName]
    [randomGeneratingSets, VariableName]
  Headline
    optional input to choose the variable name for the generated polynomials
  Description
    Text
      Put {\tt VariableName => x} for a choice of string or symbol x as an argument in
      the function @TO randomGeneratingSet@ or @TO randomGeneratingSets@
    Example 
      n=2; D=3; p=0.2;
      randomGeneratingSet(n,D,p)
      randomGeneratingSet(n,D,p,VariableName => y)
  SeeAlso
    randomGeneratingSet
    randomGeneratingSets
///

doc ///
  Key
    Strategy
    [randomGeneratingSet, Strategy]
    [randomGeneratingSets, Strategy]
  Headline
    optional input to choose the strategy for generating the monomial set
  Description
    Text
      Put {\tt Strategy => "ER"} or {\tt Strategy => "Minimal"} as an argument in the function @TO randomGeneratingSet@ or @TO randomGeneratingSets@. 
      "ER" draws random sets of monomials from the ER-type distribution B(n,D,p), while "Minimal" saves computation time by using quotient rings to exclude any non-minimal generators from the list.
  SeeAlso
    randomGeneratingSet
    randomGeneratingSets
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
    N=7;
    n=4; D=3; M={3,3,3};
    assert (N==#randomGeneratingSets(n,D,M,N))
///

TEST ///
    -- Check multiple samples agree
    n=4; D=3;
    L = randomGeneratingSets(n,D,1.0,3);
    R = ring(L#0#0);
    L = apply(L,l-> apply(l,m-> sub(m,R)));
    assert (set L#0===set L#1)
    assert (set L#0===set L#2)
    assert (set L#1===set L#2)
///

--***********************--
--  randomGeneratingSet  --
--***********************--

TEST ///
    -- Check no terms are chosen for a probability of 0
    assert (0==(randomGeneratingSet(5,5,0.0))#0)
    assert (0==(randomGeneratingSet(5,4,toList(4:0.0)))#0)
    assert (0==(randomGeneratingSet(5,4,0.0, Strategy=>"Minimal"))#0)
    assert (0==(randomGeneratingSet(5,4,toList(4:0.0), Strategy=>"Minimal"))#0)
    assert (0==(randomGeneratingSet(5,4,0))#0)
    assert (0==(randomGeneratingSet(5,4,toList(4:0)))#0)
///

TEST ///
    -- Check all possible values are outputted with a probability of 1
    n=4; D=3;
    assert (product(toList((D+1)..D+n))/n!-1==#randomGeneratingSet(n,D,1.0))
    assert (product(toList((D+1)..D+n))/n!-1==#randomGeneratingSet(n,D,{1.0,1.0,1.0}))
    n=6; D=2;
    assert (product(toList((D+1)..D+n))/n!-1==#randomGeneratingSet(n,D,1.0))
    assert (product(toList((D+1)..D+n))/n!-1==#randomGeneratingSet(n,D,{1.0,1.0}))
    n=4;D=5;
    assert (# flatten entries basis (1, QQ[x_1..x_n])==#randomGeneratingSet(n,D,1.0, Strategy=>"Minimal"))
    assert (# flatten entries basis (2, QQ[x_1..x_n])==#randomGeneratingSet(n,D,{0.0,1.0,1.0,1.0,1.0}, Strategy=>"Minimal"))
    assert (# flatten entries basis (3, QQ[x_1..x_n])==#randomGeneratingSet(n,D,{0.0,0.0,1.0,1.0,1.0}, Strategy=>"Minimal"))
    assert (# flatten entries basis (4, QQ[x_1..x_n])==#randomGeneratingSet(n,D,{0.0,0.0,0.0,1.0,1.0}, Strategy=>"Minimal"))
    assert (# flatten entries basis (5, QQ[x_1..x_n])==#randomGeneratingSet(n,D,{0.0,0.0,0.0,0.0,1.0}, Strategy=>"Minimal"))
///

TEST ///
    -- Check every monomial is generated
    L=randomGeneratingSet(2,3,1.0)
    R=ring(L#0)
    assert(set L===set {R_0,R_1,R_0^2,R_0*R_1,R_1^2,R_0^3,R_0^2*R_1,R_0*R_1^2,R_1^3})
    L=randomGeneratingSet(2,3,9)
    R=ring(L#0)
    assert(set L===set {R_0,R_1,R_0^2,R_0*R_1,R_1^2,R_0^3,R_0^2*R_1,R_0*R_1^2,R_1^3})
    L=randomGeneratingSet(3,3,{0.0,1.0,0.0})
    R=ring(L#0)
    assert(set L===set {R_0^2,R_0*R_1,R_1^2,R_0*R_2,R_1*R_2,R_2^2})
    L=randomGeneratingSet(3,3,1.0, Strategy=>"Minimal");
    R=ring(L#0);
    assert(set L===set {R_0, R_1, R_2})
    L=randomGeneratingSet(2,3,{2,3,4})
    R=ring(L#0)
    assert(set L===set {R_0,R_1,R_0^2,R_0*R_1,R_1^2,R_0^3,R_0^2*R_1,R_0*R_1^2,R_1^3})
    L=randomGeneratingSet(3,3,{0.0,1.0,1.0}, Strategy=>"Minimal");
    R=ring(L#0);
    assert(set L===set {R_0^2,R_0*R_1,R_1^2,R_0*R_2,R_1*R_2,R_2^2})
    L=randomGeneratingSet(3,3,{0.0,0.0,1.0}, Strategy=>"Minimal");
    R=ring(L#0);
    assert(set L===set {R_0^3,R_0^2*R_1,R_0^2*R_2,R_0*R_1*R_2,R_1^3,R_0*R_1^2,R_1^2*R_2,R_0*R_2^2,R_1*R_2^2,R_2^3})
///

TEST ///
    -- Check max degree of monomial less than or equal to D
    n=10; D=5;
    assert(D==max(apply(randomGeneratingSet(n,D,1.0),m->first degree m)))
    assert(D==max(apply(randomGeneratingSet(n,D,toList(D:1.0)),m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(D==max(apply(randomGeneratingSet(n,D,M),m->first degree m)))
    assert(D==max(apply((randomGeneratingSet(n,D,{0.0,0.0,0.0,0.0,1.0}, Strategy=>"Minimal"),m->first degree m))))
    n=4; D=7;
    assert(D==max(apply(randomGeneratingSet(n,D,1.0),m->first degree m)))
    assert(D==max(apply(randomGeneratingSet(n,D,toList(D:1.0)),m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(D==max(apply(randomGeneratingSet(n,D,M),m->first degree m)))
    assert(D==max(apply(randomGeneratingSet(n,D,{1,1,1,1,1,1,1}), m->first degree m)))
///

TEST ///
    -- Check min degree of monomial greater than or equal to 1
    n=8; D=6;
    assert(1==min(apply(randomGeneratingSet(n,D,1.0),m->first degree m)))
    assert(1==min(apply(randomGeneratingSet(n,D,toList(D:1.0)),m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(1==min(apply(randomGeneratingSet(n,D,M),m->first degree m)))
    n=3; D=5;
    assert(1==min(apply(randomGeneratingSet(n,D,1.0),m->first degree m)))
    assert(1==min(apply(randomGeneratingSet(n,D,toList(D:1.0)),m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(1==min(apply(randomGeneratingSet(n,D,M),m->first degree m)))
    n=10; D=5;
    assert(1==min(apply((randomGeneratingSet(n,D,1.0, Strategy=>"Minimal"),m->first degree m))))
    assert(1==min(apply((randomGeneratingSet(n,D,toList(D:1.0), Strategy=>"Minimal"),m->first degree m))))
    assert(1==min(apply(randomGeneratingSet(n,D,toList(D:1)), m->first degree m)))
///

end

You can write anything you want down here.  I like to keep examples
as I’m developing here.  Clean it up before submitting for
publication.  If you don't want to do that, you can omit the "end"
above.
