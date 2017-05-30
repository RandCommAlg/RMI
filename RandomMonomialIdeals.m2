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
    "randomMonomialIdeals",
    "Coefficient",
    "VariableName",
    "Strategy"
    }

--***************************************--
--  Exported methods 	     	     	 --
--***************************************--

ER = getSymbol "ER"
Minimal = getSymbol "Minimal"

randomGeneratingSets = method(TypicalValue => List, Options => {Coefficient => QQ,
	                                                        VariableName => "x",
								Strategy => ER})
randomGeneratingSets (ZZ,ZZ,RR,ZZ) := List => o -> (n,D,p,N) -> (
    if p<0.0 or 1.0<p then error "p expected to be a real number between 0.0 and 1.0";
    randomGeneratingSets(n,D,toList(D:p),N,
	                 Coefficient=>o.Coefficient,
			 VariableName=>o.VariableName,
			 Strategy=>o.Strategy)
)

randomGeneratingSets (ZZ,ZZ,ZZ,ZZ) := List => o -> (n,D,M,N) -> (
    if N<1 then stderr << "warning: N expected to be a positive integer" << endl;
    apply(N,i-> randomGeneratingSet(n,D,M,
	                            Coefficient=>o.Coefficient,
				    VariableName=>o.VariableName,
				    Strategy=>o.Strategy))
)

randomGeneratingSets (ZZ,ZZ,List,ZZ) := List => o -> (n,D,p,N) -> (
    if N<1 then stderr << "warning: N expected to be a positive integer" << endl;
    apply(N,i-> randomGeneratingSet(n,D,p,
	                            Coefficient=>o.Coefficient,
				    VariableName=>o.VariableName,
				    Strategy=>o.Strategy))
)

randomGeneratingSet = method(TypicalValue => List, Options => {Coefficient => QQ,
	                                                       VariableName => "x",
							       Strategy => ER})
randomGeneratingSet (ZZ,ZZ,RR) := List => o -> (n,D,p) -> (
    if p<0.0 or 1.0<p then error "p expected to be a real number between 0.0 and 1.0";
    randomGeneratingSet(n,D,toList(D:p),
	                Coefficient=>o.Coefficient,
			VariableName=>o.VariableName,
			Strategy=>o.Strategy)
)

randomGeneratingSet (ZZ,ZZ,ZZ) := List => o -> (n,D,M) -> (
    if M<0 then stderr << "warning: M expected to be a nonnegative integer" << endl;
    if o.Strategy === Minimal then error "Minimal not implemented for fixed size ER model";
    x := toSymbol o.VariableName;
    R := o.Coefficient[x_1..x_n];
    allMonomials := flatten flatten apply(toList(1..D),d->entries basis(d,R));
    take(random(allMonomials), M)
)

randomGeneratingSet (ZZ,ZZ,List) := List => o -> (n,D,p) -> (
    if n<1 then error "n expected to be a positive integer";
    if #p != D then error "p expected to be a list of length D";
    if any(p,q-> q<0.0 or 1.0<q) then error "p expected to be a list of real numbers between 0.0 and 1.0";
    x := toSymbol o.VariableName;
    R := o.Coefficient[x_1..x_n];
    B := {};
    if o.Strategy === Minimal then (
        currentRing := R;
        apply(D, d->(
            chosen := select(flatten entries basis(d+1, currentRing), m->random(0.0,1.0)<=p_d);
            B = flatten append(B, chosen/(i->sub(i, R)));
            currentRing = currentRing/promote(ideal(chosen), currentRing)
        )))
    else
        B = flatten apply(toList(1..D),d-> select(flatten entries basis(d,R),m-> random(0.0,1.0)<=p_(d-1)));
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

 randomMonomialIdeals = method(TypicalValue => List, Options => {Coefficient => QQ, VariableName => "x", IncludeZeroIdeals => false})
			
 randomMonomialIdeals (ZZ,ZZ,List,ZZ) := List => o -> (n,D,p,N) -> (
 	B:=randomGeneratingSets(n,D,p,N,Coefficient=>o.Coefficient,VariableName=>o.VariableName,Strategy=>Minimal);
	idealsFromGeneratingSets(B,p_0,D,"temporary")
	-- idealsFromGeneratingSets currently doesn't have input options for p being a list
	-- "temporary" needed since idealsFromGeneratingSets needs a String input
)
 randomMonomialIdeals (ZZ,ZZ,RR,ZZ) := List => o -> (n,D,p,N) -> (
 	B:=randomGeneratingSets(n,D,p,N,Coefficient=>o.Coefficient,VariableName=>o.VariableName,Strategy=>Minimal);
	idealsFromGeneratingSets(B,p,D,"temporary",IncludeZeroIdeals=>o.IncludeZeroIdeals)
)
 randomMonomialIdeals (ZZ,ZZ,ZZ,ZZ) := List => o -> (n,D,M,N) -> (
 	B:=randomGeneratingSets(n,D,M,N,Coefficient=>o.Coefficient,VariableName=>o.VariableName);
	idealsFromGeneratingSets(B,.5,D,"temporary",IncludeZeroIdeals=>o.IncludeZeroIdeals)
	-- need .5, since idealsFromGeneratingSets needs a RR input
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
   n=2; D=3; p=0.2; N=10;
   randomGeneratingSets(n,D,p,N)
   randomGeneratingSets(3,2,0.6,4)
  Text
   Note that this model does not generate the monomial $1$: 
  Example
   randomGeneratingSets(3,2,1.0,1)
  Text 
   If $M$ is an integer, then randomGeneratingSets creates $N$ random sets of monomials of size $M$:
   randomly select $M$ monomials from the list of all monomials of degree $1,\dots,D$ in $n$ variables.
  Example
   n=10; D=5; M=4; N=3;
   randomGeneratingSets(n,D,M,N)
  Text
   Note that each set has $M = 4$ monomials.
  Text
   If $M$ is bigger than the total number of monomials in $n$ variables of degree at most $D$, then the method will simply return all those monomials (and not $M$ of them). For example: 
  Example
   randomGeneratingSets(2,2,10,1)
  Text
   returns 5 monomials in a generating set, and not 10, since there do not exist 10 to choose from.
  Text 
   If $p=p_1,\dots,p_D$ is a list of real numbers of length $D$, then randomGeneratingSets generates the sets utilizing the graded Erdos-Renyi-type model:
   select each monomial of degree $1\le d\le D$, independently, with probability $p_d$.
  Example
   p={0.0, 1.0, 1.0}; 
   randomGeneratingSets(2,3,p,1)
  Text
   Note that the degree-1 monomials were not generated, since the first probability vector entry is 0.
///

doc ///
 Key
  randomMonomialIdeals
  (randomMonomialIdeals,ZZ,ZZ,RR,ZZ)
  (randomMonomialIdeals,ZZ,ZZ,ZZ,ZZ)
  (randomMonomialIdeals,ZZ,ZZ,List,ZZ)
 Headline
  randomly generates mononial ideals, with each monomial up to a given degree
 Usage
  randomMonomialIdeals(ZZ,ZZ,RR,ZZ)
  randomMonomialIdeals(ZZ,ZZ,ZZ,ZZ)
  randomMonomialIdeals(ZZ,ZZ,List,ZZ)
 Inputs
  n: ZZ
    number of variables
  D: ZZ
    maximum degree
  p: RR
     or @ofClass List@
     , probability to select a monomial
  M: ZZ
     maximum number of monomials in each generating set for the ideal
  N: ZZ
    number of ideals generated
 Outputs
  B: List
   random monomial ideals
 Description
  Text
   randomMonomialIdeals creates $N$ random monomial ideals, with each monomial having degree $d$, $1\leq d\leq D$, in $n$ variables. 
   If $p$ is a real number, it generates each of these ideals according to the Erdos-Renyi-type model: 
   from the list of all monomials of degree $1,\dots,D$ in $n$ variables, it selects each one, independently, with probability $p$. 
  Example
   n=2; D=3; p=0.2; N=10;
   randomMonomialIdeals(n,D,p,N)
   randomMonomialIdeals(5,3,0.4,4)
  Text
   Note that this model does not generate the monomial $1$: 
  Example
   randomMonomialIdeals(3,2,1.0,1)
  Text 
   If $M$ is an integer, then randomMonomialIdeals creates $N$ random monomial ideals of size at most $M$:
   randomly select $M$ monomials from the list of all monomials of degree $1,\dots,D$ in $n$ variables, then generate its ideal.
  Example
   n=8; D=4; M=7; N=3;
   randomMonomialIdeals(n,D,M,N)
  Text
   Note that each generating set of each ideal has at most $M = 7$ monomials. If one monomial divides another monomial that was generated, it will not be in the generating set
  Text 
   If $p=p_1,\dots,p_D$ is a list of real numbers of length $D$, then randomMonomialIdeals generates the generating sets utilizing the graded Erdos-Renyi-type model:
   select each monomial of degree $1\le d\le D$, independently, with probability $p_d$.
  Example
   p={0.0, 1.0, 1.0}; 
   randomGeneratingSets(2,3,p,1)
  Text
   Note that the degree-1 monomials were not generated, since the first probability vector entry is 0.
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
    assert (0==(randomGeneratingSet(5,4,0.0, Strategy=>Minimal))#0)
    assert (0==(randomGeneratingSet(5,4,toList(4:0.0), Strategy=>Minimal))#0)
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
    assert (# flatten entries basis (1, QQ[x_1..x_n])==#randomGeneratingSet(n,D,1.0, Strategy=>Minimal))
    assert (# flatten entries basis (2, QQ[x_1..x_n])==#randomGeneratingSet(n,D,{0.0,1.0,1.0,1.0,1.0}, Strategy=>Minimal))
    assert (# flatten entries basis (3, QQ[x_1..x_n])==#randomGeneratingSet(n,D,{0.0,0.0,1.0,1.0,1.0}, Strategy=>Minimal))
    assert (# flatten entries basis (4, QQ[x_1..x_n])==#randomGeneratingSet(n,D,{0.0,0.0,0.0,1.0,1.0}, Strategy=>Minimal))
    assert (# flatten entries basis (5, QQ[x_1..x_n])==#randomGeneratingSet(n,D,{0.0,0.0,0.0,0.0,1.0}, Strategy=>Minimal))
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
    L=randomGeneratingSet(3,3,1.0, Strategy=>Minimal);
    R=ring(L#0);
    assert(set L===set {R_0, R_1, R_2})
    L=randomGeneratingSet(3,3,{0.0,1.0,1.0}, Strategy=>Minimal);
    R=ring(L#0);
    assert(set L===set {R_0^2,R_0*R_1,R_1^2,R_0*R_2,R_1*R_2,R_2^2})
    L=randomGeneratingSet(3,3,{0.0,0.0,1.0}, Strategy=>Minimal);
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
    assert(D==max(apply((randomGeneratingSet(n,D,{0.0,0.0,0.0,0.0,1.0}, Strategy=>Minimal),m->first degree m))))
    n=4; D=7;
    assert(D==max(apply(randomGeneratingSet(n,D,1.0),m->first degree m)))
    assert(D==max(apply(randomGeneratingSet(n,D,toList(D:1.0)),m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(D==max(apply(randomGeneratingSet(n,D,M),m->first degree m)))
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
    assert(1==min(apply((randomGeneratingSet(n,D,1.0, Strategy=>Minimal),m->first degree m))))
    assert(1==min(apply((randomGeneratingSet(n,D,toList(D:1.0), Strategy=>Minimal),m->first degree m))))
///

--************************--
--  randomMonomialIdeals  --
--************************--

TEST ///
  -- check the number of ideals
  n=5; D=5; p=.6; N=3;
  B = flatten randomMonomialIdeals(n,D,p,N);
  assert ((N+1)===#B)
///

TEST ///
  -- check the number of monomials in the generating set of the ideal
  n=4; D=6; M=7; N=1;
  B = flatten randomMonomialIdeals(n,D,p,N);
  assert (M>=numgens B_0)
///
end

You can write anything you want down here.  I like to keep examples
as I’m developing here.  Clean it up before submitting for
publication.  If you don't want to do that, you can omit the "end"
above.
