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
needsPackage "Depth";

needsPackage "BoijSoederberg"

export {
    "randomMonomialSets",
    "randomMonomialSet",
    "idealsFromGeneratingSets",
    "randomMonomialIdeals",
    "Coefficients",
    "VariableName",
    "mingenStats",
    "IncludeZeroIdeals",
    "dimStats",
    "regStats",
    "CMStats",
    "borelFixedStats",
    "ShowTally",
    "degStats",
    "bettiStats",
    "SaveBettis"
    "pdimStats"
}

--***************************************--
--  Exported methods 	     	     	 --
--***************************************--

randomMonomialSets = method(TypicalValue => List, Options => {Coefficients => QQ,
	                                                        VariableName => "x",
								Strategy => "ER"})
randomMonomialSets (ZZ,ZZ,RR,ZZ) := List => o -> (n,D,p,N) -> (
    if p<0.0 or 1.0<p then error "p expected to be a real number between 0.0 and 1.0";
    randomMonomialSets(n,D,toList(D:p),N,o)
)

randomMonomialSets (PolynomialRing,ZZ,RR,ZZ) := List => o -> (R,D,p,N) -> (
    if p<0.0 or 1.0<p then error "p expected to be a real number between 0.0 and 1.0";
    randomMonomialSets(R,D,toList(D:p),N,o)
)

randomMonomialSets (ZZ,ZZ,ZZ,ZZ) := List => o -> (n,D,M,N) -> (
    if N<1 then stderr << "warning: N expected to be a positive integer" << endl;
    x := toSymbol o.VariableName;
    R := o.Coefficients[x_1..x_n];
    apply(N,i-> randomMonomialSet(R,D,M,o))
)

randomMonomialSets (PolynomialRing,ZZ,ZZ,ZZ) := List => o -> (R,D,M,N) -> (
    if N<1 then stderr << "warning: N expected to be a positive integer" << endl;
    apply(N,i-> randomMonomialSet(R,D,M,o))
)

randomMonomialSets (ZZ,ZZ,List,ZZ) := List => o -> (n,D,pOrM,N) -> (
    if n<1 then error "n expected to be a positive integer";
    if N<1 then stderr << "warning: N expected to be a positive integer" << endl;
    x := toSymbol o.VariableName;
    R := o.Coefficients[x_1..x_n];
    apply(N,i-> randomMonomialSet(R,D,pOrM,o))
)

randomMonomialSets (PolynomialRing,ZZ,List,ZZ) := List => o -> (R,D,pOrM,N) -> (
    if N<1 then stderr << "warning: N expected to be a positive integer" << endl;
    apply(N,i-> randomMonomialSet(R,D,pOrM,o))
)

randomMonomialSet = method(TypicalValue => List, Options => {Coefficients => QQ,
	                                                       VariableName => "x",
							       Strategy => "ER"})
randomMonomialSet (ZZ,ZZ,RR) := List => o -> (n,D,p) -> (
    if p<0.0 or 1.0<p then error "p expected to be a real number between 0.0 and 1.0";
    randomMonomialSet(n,D,toList(D:p),o)
)

randomMonomialSet (PolynomialRing,ZZ,RR) := List => o -> (R,D,p) -> (
    if p<0.0 or 1.0<p then error "p expected to be a real number between 0.0 and 1.0";
    randomMonomialSet(R,D,toList(D:p),o)
)

randomMonomialSet (ZZ,ZZ,ZZ) := List => o -> (n,D,M) -> (
    if n<1 then error "n expected to be a positive integer";
    x := toSymbol o.VariableName;
    R := o.Coefficients[x_1..x_n];
    randomMonomialSet(R,D,M)
)

randomMonomialSet (PolynomialRing,ZZ,ZZ) := List => o -> (R,D,M) -> (
    if M<0 then stderr << "warning: M expected to be a nonnegative integer" << endl;
    if o.Strategy === "Minimal" then error "Minimal not yet implemented for fixed size ER model";
    allMonomials := flatten flatten apply(toList(1..D),d->entries basis(d,R));
    C := take(random(allMonomials), M);
    if C==={} then {0_R} else C
)

randomMonomialSet (ZZ,ZZ,List) := List => o -> (n,D,pOrM) -> (
    if n<1 then error "n expected to be a positive integer";
    x := toSymbol o.VariableName;
    R := o.Coefficients[x_1..x_n];
    randomMonomialSet(R,D,pOrM,o)
)

randomMonomialSet (PolynomialRing,ZZ,List) := List => o -> (R,D,pOrM) -> (
    if #pOrM != D then error "pOrM expected to be a list of length D";
    if not all(pOrM, q->instance(q, ZZ)) and not all(pOrM, q->instance(q,RR)) then error "pOrM must be a list of all integers or all real numbers";
    B := {};
    if all(pOrM,q->instance(q,ZZ)) then (
        if o.Strategy === "Minimal" then error "Minimal not implemented for fixed size ER model";
        B = flatten apply(toList(1..D), d->take(random(flatten entries basis(d,R)), pOrM_(d-1)));
    )
    else if all(pOrM,q->instance(q,RR)) then (
        if any(pOrM,q-> q<0.0 or 1.0<q) then error "pOrM expected to be a list of real numbers between 0.0 and 1.0";
        if o.Strategy === "Minimal" then (
            currentRing := R;
            apply(D, d->(
                chosen := select(flatten entries basis(d+1, currentRing), m->random(0.0,1.0)<=pOrM_d);
                B = flatten append(B, chosen/(i->sub(i, R)));
                currentRing = currentRing/promote(ideal(chosen), currentRing)
            )))
        else
            B = flatten apply(toList(1..D),d-> select(flatten entries basis(d,R),m-> random(0.0,1.0)<=pOrM_(d-1)));
	);
    B = apply(B,m->sub(m,R));
    if B==={} then {0_R} else B
)

bettiStats = method(TypicalValue =>Sequence, Options =>{IncludeZeroIdeals=>true, SaveBettis => false})
bettiStats List :=  o-> (ideals) -> ( 
    if o.SaveBettis then (
	basefilename:="DELETEme"; fileNameExt:="GetRidOfThis";
	filename1 := concatenate(basefilename,"bettis",fileNameExt)
	);
    -- sum of the betti tables and betti shapes: 
    N := #ideals; Z:=0;
    if not o.IncludeZeroIdeals then (ideals,Z) = extractNonzeroIdeals(ideals);
    beta := new BettiTally; 
    betaShapes := new BettiTally;
    pure := 0; -- count pure Betti tables. THIS IS NOT CURRENTLY RETURNED BY THE METHOD. DECIDE. 
    -- add up all the betti tables: 
    apply(#ideals,i->( 
        resi := betti res ideals_i;
	if isPure resi then pure = pure +1; -- IS there a shorter way to do this? why can't i run "value(true)" and get 1? :-/
        if o.SaveBettis then filename1 << resi << endl;
    	-- add Betti tables as we compute them, so as to not store an entire sample of Betti tables in memory:
        beta = beta + resi;
  	-- let's only keep a 1 in all spots where ther was a non-zero Betti number: 
	beta1mtx := matrix(resi);
	Rtemp := (ring ideals_i)^1/ideals_i;
	beta1shape := new BettiTally from mat2betti  matrix pack(1+pdim(Rtemp), apply(flatten entries beta1mtx, i-> if i>0 then i=1 else i=0));
	betaShapes = betaShapes + beta1shape
	)
    );
    if o.SaveBettis then filename1 << close;
    -- compute the average Betti table shape: 
    bShape := mat2betti(1/#ideals*(sub(matrix(betaShapes), RR)));
    -- compute the average Betti table:
    b := mat2betti(1/#ideals*(sub(matrix(beta), RR)));
    return (bShape,b)--,pure)
    )
    
degStats = method(TypicalValue =>Sequence, Options =>{ShowTally => false})
degStats List :=  o-> (listOfIdeals) -> (
    N := #listOfIdeals;
    deg := 0;
    degHistogram:={};
    apply(#listOfIdeals, i->( 
        degi := degree listOfIdeals_i;
        degHistogram = append(degHistogram, degi)
	)
    );
    ret:=();
    avg:=sub(1/N*(sum degHistogram), RR);
    Ex2:=sub(1/N*(sum apply(elements(tally degHistogram), i->i^2)), RR);
    var:= Ex2 - avg^2;
    stdDev:= var^(1/2);
    if o.ShowTally
    	then(ret=(avg, stdDev,tally degHistogram); return ret;);
    print "Average Degree:" expression(sub(1/N*(sum degHistogram), RR));
    ret = (avg, stdDev)
)

--creates a list of monomialIdeal objects from a list of monomial generating sets 
idealsFromGeneratingSets =  method(TypicalValue => List, Options => {IncludeZeroIdeals => true})
idealsFromGeneratingSets(List):= o -> (B) -> (
    N := # B;
    n := numgens ring ideal B#0; -- ring of the first monomial in the first gen set
    ideals := {};
    for i from 0 to #B-1 do {
	ideals = B / (b-> monomialIdeal b);
	};
    (nonzeroIdeals,numberOfZeroIdeals) := extractNonzeroIdeals(ideals);
    print(concatenate("There are ", toString(#B)," ideals in this sample."));
    print(concatenate("Of those, ", toString numberOfZeroIdeals, " were the zero ideal."));
    if o.IncludeZeroIdeals then return ideals else return (nonzeroIdeals,numberOfZeroIdeals); 
)

dimStats = method(TypicalValue => Sequence, Options => {ShowTally => false})
dimStats List := o-> (listOfIdeals) -> (
    N := #listOfIdeals;
    dims:=0;
    dimsHistogram:={};
    apply(#listOfIdeals,i->( 
        dimi := dim listOfIdeals_i;
    dimsHistogram = append(dimsHistogram, dimi)
    )
    );
    ret:= ();
    avg:=sub(1/N*(sum dimsHistogram), RR);
    Ex2:=sub(1/N*(sum apply(elements(tally dimsHistogram), i->i^2)), RR);
    var:= Ex2 - avg^2;
    stdDev:= var^(1/2);
    if o.ShowTally 
         then(ret = (avg, stdDev, tally dimsHistogram), return ret;);
    print "Average Krull dimension:" expression(sub(1/N*(sum dimsHistogram), RR));
    ret = (avg, stdDev)
)

regStats = method(TypicalValue => Sequence, Options => {ShowTally => false})
regStats List := o-> (listOfIdeals) -> (
    N:=#listOfIdeals;
    ideals := extractNonzeroIdeals(listOfIdeals);
    ideals = ideals_0;
    reg := 0;
    ret := ();
    regHistogram:={};
    if set {} === set ideals then (
	regHistogram = N:-infinity;
	stdDev := 0;
	if o.ShowTally then(
	    print "All ideals in this list are the zero ideal";
	    ret=(-infinity, 0, tally regHistogram); 
	    return ret;
	    );
	print "All ideals in this list are the zero ideal";
	ret = (-infinity, 0)
    )
    else (
	apply(#ideals,i->( 
              regi := regularity ideals_i;
              regHistogram = append(regHistogram, regi)
	     ));
             avg := sub(1/#ideals*(sum regHistogram), RR);
    	     Ex2 := sub((1/(#ideals))*(sum apply(elements(tally regHistogram), i->i^2)), RR);
    	     var := Ex2-avg^2;
    	     stdDev = var^(1/2);
    	     if o.ShowTally
    	        then(ret=(avg, stdDev,tally regHistogram); return ret;);
	     print(concatenate(toString(N-#ideals), " zero ideals were extracted from this sample"));
    	     ret = (avg, stdDev)
         )
    
)

 randomMonomialIdeals = method(TypicalValue => List, Options => {Coefficients => QQ, VariableName => "x", IncludeZeroIdeals => true})
			
 randomMonomialIdeals (ZZ,ZZ,List,ZZ) := List => o -> (n,D,pOrM,N) -> (
        B:={};
        if all(pOrM,q->instance(q,RR)) then 
	    B=randomMonomialSets(n,D,pOrM,N,Coefficients=>o.Coefficients,VariableName=>o.VariableName,Strategy=>"Minimal")
	else if all(pOrM,q->instance(q,ZZ)) then 
	    B=randomMonomialSets(n,D,pOrM,N,Coefficients=>o.Coefficients,VariableName=>o.VariableName);
	idealsFromGeneratingSets(B,IncludeZeroIdeals=>o.IncludeZeroIdeals)
)
 randomMonomialIdeals (ZZ,ZZ,RR,ZZ) := List => o -> (n,D,p,N) -> (
 	B:=randomMonomialSets(n,D,p,N,Coefficients=>o.Coefficients,VariableName=>o.VariableName,Strategy=>"Minimal");
	idealsFromGeneratingSets(B,IncludeZeroIdeals=>o.IncludeZeroIdeals)
)
 randomMonomialIdeals (ZZ,ZZ,ZZ,ZZ) := List => o -> (n,D,M,N) -> (
 	B:=randomMonomialSets(n,D,M,N,Coefficients=>o.Coefficients,VariableName=>o.VariableName);
	idealsFromGeneratingSets(B,IncludeZeroIdeals=>o.IncludeZeroIdeals)
)
--checks if each RMI is CM and prints the % CM as a real number
CMStats = method(TypicalValue => RR)
CMStats (List) :=  (listOfIdeals) -> (
    cm := 0;
    N:= #listOfIdeals;
    R := ring(listOfIdeals#0);
    for i from 0 to #listOfIdeals-1 do (
       if isCM(R/listOfIdeals_i) == true then cm = cm + 1 else cm = cm);
     print "Percent Cohen-Macaulay:" expression(sub((cm)/N, RR));
   sub((cm)/N, RR)
)

--checks whether each RMI is Borel-fixed, 
--prints and returns % of Borel-fixed RMIs in sample (real number) 
borelFixedStats = method()
borelFixedStats (List) :=  (ideals) -> (
    bor := 0;
    N:=#ideals;
    for i from 0 to #ideals-1 do ( 
        if isBorel((ideals_i)) == true then bor = bor + 1 else bor = bor);     
    print "Percent Borel-fixed:" expression(sub((bor)/N, RR));
    sub((bor)/N, RR)
)
mingenStats = method(TypicalValue => Sequence, Options => {ShowTally => false})
mingenStats (List) :=  o -> (ideals) -> (
    N:=#ideals;
    ideals = extractNonzeroIdeals(ideals);
    ideals = ideals_0;
    num := 0;
    numgensHist := {};
    m := 0;
    complexityHist := {};
    ret:=();
    if set {} === set ideals then (
        numgensHist = N:-infinity;
	complexityHist = N:-infinity;
	numStdDev := 0;
	comStdDev := 0;
	if o.ShowTally then(ret=(-infinity, 0, tally numgensHist, -infinity, 0, tally complexityHist); return ret;);
	print "Average # of min gens:" expression(-infinity);
	print "Average degree complexity:" expression(-infinity);
	ret = (-infinity, 0, -infinity, 0)
    )
    else (
        apply(#ideals,i->( 
            mingensi := gens gb ideals_i;
            numgensi := numgens source mingensi; 
            mi := max({degrees(mingensi)}#0#1); 
	    numgensHist = append(numgensHist, numgensi); 
	    complexityHist = append(complexityHist, mi#0)
	    )
        );
    numAvg:=sub((1/(#ideals))*(sum numgensHist), RR);
    comAvg:=sub((1/(#ideals))*(sum complexityHist), RR);
    numEx2:=sub((1/(#ideals))*(sum apply(elements(tally numgensHist), i->i^2)), RR);
    comEx2:=sub((1/(#ideals))*(sum apply(elements(tally complexityHist), i->i^2)), RR);
    numVar:= numEx2 - numAvg^2;
    comVar:= comEx2 - comAvg^2;
    numStdDev= numVar^(1/2);
    comStdDev= comVar^(1/2);
    if o.ShowTally 
       then(ret=(numAvg, numStdDev, tally numgensHist, comAvg, comStdDev, tally complexityHist); return ret;); 
    print "Average # of min gens:" expression(sub((1/(#ideals))*(sum numgensHist), RR));
    print "Average degree complexity:" expression(sub((1/(#ideals))*(sum complexityHist), RR));
    ret = (numAvg, numStdDev, comAvg, comStdDev)
    )
)
-- example to run this^^ right now: 
-- L=randomMonomialIdeals(3,4,0.5,2)
-- (mu,reg) = mingenStats(L);
-- mu
-- reg


pdimStats = method(TypicalValue=>Sequence, Options => {ShowTally => false})
pdimStats (List) := o-> (ideals) -> (
    N:=#ideals;
    pdHist:={};
    R:=ring(ideals_0);
    apply(#ideals,i-> 
	(
        pdimi := pdim(R^1/ideals_i);
	pdHist = append(pdHist, pdimi)
	)
    );           
    ret:=();
    avg:=sub(((1/N)*(sum pdHist)),RR);
    Ex2:=sub(((1/N)*(sum apply(elements(tally pdHist), i->i^2))), RR);
    var:= Ex2 - avg^2;
    stdDev:= var^(1/2);
    if o.ShowTally 
         then(ret = (avg, stdDev, tally pdHist), return ret;);
    ret=(avg, stdDev) 
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
  A package for generating Erdos-Renyi-type random monomial ideals and variations.
 Description
  Text
   {\em RandomMonomialIdeals} is a  package that... 
  -- Caveat
  -- Still trying to figure this out. [REMOVE ME]
///

doc ///
 Key
  randomMonomialSets
  (randomMonomialSets,ZZ,ZZ,RR,ZZ)
  (randomMonomialSets,PolynomialRing,ZZ,RR,ZZ)
  (randomMonomialSets,ZZ,ZZ,ZZ,ZZ)
  (randomMonomialSets,PolynomialRing,ZZ,ZZ,ZZ)
  (randomMonomialSets,ZZ,ZZ,List,ZZ)
  (randomMonomialSets,PolynomialRing,ZZ,List,ZZ)
 Headline
  randomly generates lists of monomials in fixed number of variables up to a given degree
 Usage
  randomMonomialSets(ZZ,ZZ,RR,ZZ)
  randomMonomialSets(PolynomialRing,ZZ,RR,ZZ)
  randomMonomialSets(ZZ,ZZ,ZZ,ZZ)
  randomMonomialSets(PolynomialRing,ZZ,ZZ,ZZ)
  randomMonomialSets(ZZ,ZZ,List,ZZ)
  randomMonomialSets(PolynomialRing,ZZ,List,ZZ)
 Inputs
  n: ZZ
    number of variables
  : PolynomialRing
    the ring in which monomial sets are to live if n is not specified
  D: ZZ
    maximum degree
  p: RR
     the probability of selecting a monomial, 
  M: ZZ
     number of monomials in the set, 
  : List 
     of real numbers whose i-th entry is the probability of selecing a monomial of degree i, 
     or of integers whose i-th entry is the number of monomials of degree i in each set
  N: ZZ
    number of sets generated
 Outputs
  B: List
   random generating sets of monomials
 Description
  Text
   randomMonomialSets creates $N$ random sets of monomials of degree $d$, $1\leq d\leq D$, in $n$ variables. 
   It does so by calling @TO randomMonomialSet@ $N$ times.
 SeeAlso
  randomMonomialSet
///

doc ///
 Key
  bettiStats
  (bettiStats,List)
 Headline
  statistics on Betti tables of a sample of monomial ideals
 Usage
  bettiStats(List)
 Inputs
  L: List
   of @TO monomialIdeal@s
 Outputs
  : Sequence
   of BettyTallies, representing the mean Betti table shape and the mean Betti table of the elements in the list {\tt L}.
 Description
  Text
   For a sample of ideals stored as a List, this method computes some basic Betti table statistics of the sample.
   Namely, it computes the average shape of the Betti tables (where 1 is recorded in entry (ij) for each element if beta_{ij} is not zero), 
   and it also computes the average Betti table (think of beta_{ij} as the mean value of beta_{ij} for all ideals in the sample). 
  Example
   R = ZZ/101[a..e];
   L={monomialIdeal"a2b,bc", monomialIdeal"ab,bc3",monomialIdeal"ab,ac,bd,be,ae,cd,ce,a3,b3,c3,d3,e3"}
   (meanBettiShape,meanBetti) = bettiStats L;
   meanBettiShape
   meanBetti
  Text
   For sample size $N$, the average Betti table {\em shape} simply considers nonzero Betti numbers. It is to be interpreted as follows:
   entry (i,j) encodes the following sum of indicators: 
   $\sum_{all ideals} 1_{beta_{ij}>0} / N$; that is,
   the proportion of ideals with a nonzero beta_{ij}.
   Thus an entry of 0.33 means 33% of ideals have a non-zero Betti number there.
  Example
   apply(L,i->betti res i)
   meanBettiShape   
  Text 
   For sample size $N$, the average Betti table is to be interpreted as follows: 
   entry $(i,j)$ encodes  $\sum_{I\in ideals}beta_{ij}(R/I) / N$:
  Example
   apply(L,i->betti res i)
   meanBetti
  Text 
   Note that this method will work on a List of any objects to which @TO betti@ @TO res@ can be applied. 

  Text 
   If {\tt ideals} contains zero ideals, you may wish to exclude them from the Betti statistics. 
   In this case, use the optional input as follows: 
  Example
   L={monomialIdeal (a^2*b,b*c), monomialIdeal(a*b,b*c^3),monomialIdeal 0_R};
   apply(L,i->betti res i)
   bettiStats(L,IncludeZeroIdeals=>false)
   
  Text 
   For now, the method also has an option to save all betti tables to a file (because they take forefer to compute for some ideals; this is how: 
   [** but this is really going to move to option documentation, if we even decide to keep it. **] 
  Example
   bettiStats (L,SaveBettis=>true) -- saves 1 file (for now). 
///

doc ///
 Key
  degStats
  (degStats,List)
 Headline
  statistics on the degrees of a list of monomialIdeals
 Usage
  degStats(List)
 Inputs
  listOfIdeals: List
   of @TO monomialIdeal@s
 Outputs
  : Sequence
   whose first entry is the average degree of a list of monomialIdeals, second entry is the standard deviation of the degree, and third entry (if option turned on) is the degree tally
 Description
  Text
   degStats finds the average and the standard deviation of the degree of R/I for a list of monomialIdeals.
   The degree of each ideal is calculated using the @TO degree@ function.
   It has the optional input of ShowTally.
  Example
   L=randomMonomialSet(3,3,1.0);
   R=ring(L#0);
   listOfIdeals={monomialIdeal(R_0^5*R_1^2,R_2),monomialIdeal(R_0,R_1,R_2),monomialIdeal(R_0^3*R_1^5,R_1^4*R_2,R_0^2*R_2^3)};
   degStats(listOfIdeals)
  Text
   The following expamples use the existing functions @TO randomMonomialSets@ and @TO idealsFromGeneratingSets@ or @TO randomMonomialIdeals@ to automatically generate a list of ideals, rather than creating the list manually:
  Example
   listOfIdeals = idealsFromGeneratingSets(randomMonomialSets(4,3,1.0,3));
   degStats(listOfIdeals)
  Example
   listOfIdeals = randomMonomialIdeals(4,3,1.0,3);
   degStats(listOfIdeals)
  Text
   Note that this function can be run with a list of any objects to which @TO degree@ can be applied.
///

doc ///
 Key
  randomMonomialIdeals
  (randomMonomialIdeals,ZZ,ZZ,RR,ZZ)
  (randomMonomialIdeals,ZZ,ZZ,ZZ,ZZ)
  (randomMonomialIdeals,ZZ,ZZ,List,ZZ)
 Headline
  randomly generates monomial ideals in fixed number of variables, with each monomial up to a given degree
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
     probability to select a monomial in the ER model, 
     or @ofClass List@ of probabilities of selecting monomials in each degree for the graded ER model
  M: ZZ
     maximum number of monomials in each generating set for the ideal
  N: ZZ
    number of ideals generated
 Outputs
  : List
   list of randomly generated @TO monomialIdeal@, and the number of zero ideals removed, if any
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
   randomly select $M$ monomials from the list of all monomials of degree $1,\dots,D$ in $n$ variables, then generate the ideal from this set.
  Example
   n=8; D=4; M=7; N=3;
   randomMonomialIdeals(n,D,M,N)
  Text
   Note that each generating set of each ideal has at most $M = 7$ monomials. If one monomial divides another monomial that was generated, it will not be in the generating set.
  
   The input of type @TO List@ controls the number of monomials in the generating set of each degree for the graded ER model.
   Specifically, this input is either a list of real numbers between 0 and 1, inclusive, whose i-th entry is 
   the probability of including a monomial of degree i in the monomial set, or it is a list of nonnegative 
   integers whose i-th entry is the number of monomials of each degree to include in the monomial set. Consider the following two examples: 
   If $p=p_1,\dots,p_D$ is a list of real numbers of length $D$, then randomMonomialSet generates the set utilizing the graded Erdos-Renyi-type model:
   select each monomial of degree $1\le d\le D$, independently, with probability $p_d$.
  Example
   randomMonomialIdeals(2,3,{0.0, 1.0, 1.0},1)
  Text
   Note that the degree-1 monomials were not generated, since the first probability vector entry is 0.
 
  
   If $M=M_1,\dots,M_D$ is a list of integers of length $D$, then randomMonomialIdeal creates a list of @TO MonomialIdeal@, where at most $M_d$ monomials are of degree $d$.
  Example
   randomMonomialIdeals(3,3,{1,1,1},1)
  Text
   Observe that there are at most one degree-1 monomials, one degree-2 monomial, and one degree-3 monomial.
 Caveat
  Since the method returns a list of @TO MonomialIdeal@s, only the minimal generating set will be displayed. 
  In contrast, @TO randomMonomialSet@ will display the full (not necessarily minimal) generating set produced by the model.
 SeeAlso
   randomMonomialSets
   idealsFromGeneratingSets
///

doc ///
 Key
  randomMonomialSet
  (randomMonomialSet,ZZ,ZZ,RR)
  (randomMonomialSet,PolynomialRing,ZZ,RR)
  (randomMonomialSet,ZZ,ZZ,ZZ)
  (randomMonomialSet,PolynomialRing,ZZ,ZZ)
  (randomMonomialSet,ZZ,ZZ,List)
  (randomMonomialSet,PolynomialRing,ZZ,List)
 Headline
  randomly generates a list of monomials in fixed number of variables up to a given degree
 Usage
  randomMonomialSet(ZZ,ZZ,RR)
  randomMonomialSet(PolynomialRing,ZZ,RR)
  randomMonomialSet(ZZ,ZZ,ZZ)
  randomMonomialSet(PolynomialRing,ZZ,ZZ)
  randomMonomialSet(ZZ,ZZ,List)
  randomMonomialSet(PolynomialRing,ZZ,List)
 Inputs
  n: ZZ
    number of variables
  : PolynomialRing
    the ring in which monomial sets are to live if n is not specified
  D: ZZ
    maximum degree
  p: RR
     the probability of selecting a monomial, 
  M: ZZ
     number of monomials in the set, 
  : List 
     of real numbers whose i-th entry is the probability of selecing a monomial of degree i, 
     or of integers whose i-th entry is the number of monomials of degree i in each set
 Outputs
  : List
   random set of monomials
 Description
  Text
   randomMonomialSet creates a list of monomials, up to a given degree $d$, $1\leq d\leq D$, in $n$ variables. 
   If $p$ is a real number, it generates the set according to the Erdos-Renyi-type model, that is, based on a Binomial distribution:
   from the list of all monomials of degree $1,\dots,D$ in $n$ variables, it selects each one, independently, with probability $p$.
  Example
   n=2; D=3; p=0.2;
   randomMonomialSet(n,D,p)
   randomMonomialSet(3,2,0.6)
  Text
   Note that this model does not generate the monomial $1$:
  Example
   randomMonomialSet(3,2,1.0)
  Text
   If $M$ is an integer, then randomMonomialSet creates a list of monomials of size $M$:
   randomly select $M$ monomials from the list of all monomials of degree $1,\dots,D$ in $n$ variables.
  Example
   n=10; D=5; M=4;
   randomMonomialSet(n,D,M)
  Text
   Note that it returns a set with $M = 4$ monomials.
  Text
   If $M$ is greater than the total number of monomials in $n$ variables of degree at most $D$, then the method will simply return all those monomials (and not $M$ of them). For example:
  Example
   randomMonomialSet(2,2,10)
  Text
   returns 5 monomials in a generating set, and not 10, since there are fewer than 10 monomials to choose from. 

   The input of type @TO List@ controls the number of monomials in the generating set of each degree for the graded ER model.
   Specifically, this input is either a list of real numbers between 0 and 1, inclusive, whose i-th entry is 
   the probability of including a monomial of degree i in the monomial set, or it is a list of nonnegative 
   integers whose i-th entry is the number of monomials of each degree to include in the monomial set. Consider the following two examples: 
   If $p=p_1,\dots,p_D$ is a list of real numbers of length $D$, then randomMonomialSet generates the set utilizing the graded Erdos-Renyi-type model:
   select each monomial of degree $1\le d\le D$, independently, with probability $p_d$.
  Example
   randomMonomialSet(2,3,{0.0, 1.0, 1.0})
  Text
   Note that the degree-1 monomials were not generated, since the first probability vector entry is 0.
  Text
   If $M=M_1,\dots,M_D$ is a list of integers of length $D$, then randomMonomialSet creates a list of monomials, where $M_d$ monomials are of degree $d$.
  Example
   randomMonomialSet(2,3,{2,1,1})
  Text
   Observe that there are two degree-1 monomials, one degree-2 monomial, and one degree-3 monomial.
   
   Sometimes we are already working in a specific ring and would like the random sets of monomials to live in the same ring:
  Example
   D=3;p=.5; R=ZZ/101[a,b,c];
   randomMonomialSet(R,D,p)
   ring oo_0
 SeeAlso
   randomMonomialSets
///

doc ///
 Key
  idealsFromGeneratingSets
  (idealsFromGeneratingSets, List)
 Headline
  creates ideals from sets of monomials
 Usage
  idealsFromGeneratingSets(List)
 Inputs
  B: List
    of sets of monomials
 Outputs
  : List
    of @TO monomialIdeal@s
 Description
  Text
   idealsFromGeneratingSets takes a list of sets of monomials and converts each set into a monomial ideal. It counts how many sets are given, and how many sets are converted to the zero ideal.
  Example
   n=4; D=2; p=1.0; N=3;
   B=randomMonomialSets(n,D,p,N); B/print
   idealsFromGeneratingSets(B)
 SeeAlso
  randomMonomialIdeals
///

doc ///
 Key
  mingenStats
  (mingenStats, List)
 Headline
  statistics on the minimal generators of a list of monomialIdeals: number and degree complexity 
 Usage
  mingenStats(List)
 Inputs
  ideals: List
    of @TO monomialIdeal@s
 Outputs
  : Sequence
    with the following entries: the average number of minimal generators, the standard deviation of the number of minimial generators, the average degree complexity, and the stadard deviation of the degree complexity. 
    If ShowTally is turned on, then the output sequence also includes the tallies of the two numbers following their standard deviation. 
 Description
  Text
   mingenStats removes zero ideals from the list of ideals, then calculates the average and the standard deviation for the number of minimal generators and degree complexity of the list of nonzero ideals.
  Example
   n=4; D=3; p={0.0,1.0,0.0}; N=3;
   B=randomMonomialIdeals(n,D,p,N);
   mingenStats(B)
  Text
   If the list given is a list of all zero ideals, mingenStats returns -infinity for the mean number of minimal generators and mean degree complexity.
  Example
   B=randomMonomialIdeals(3,3,0.0,1);
   mingenStats(B)
  Text
   Note that this function can be called on a list of @TO Ideal@ objects instead.
 Caveat
  mingenStats removes zero ideals from the list of ideals before computing the two values.
///

doc ///
  Key
    Coefficients
    [randomMonomialSet, Coefficients]
    [randomMonomialSets, Coefficients]
    [randomMonomialIdeals, Coefficients]
  Headline
    optional input to choose the coefficients of the ambient polynomial ring
  Description
    Text
      Put {\tt Coefficients => r} for a choice of field r as an argument in
      the function @TO randomMonomialSet@ or @TO randomMonomialSets@. 
    Example 
      n=2; D=3; p=0.2;
      randomMonomialSet(n,D,p)
      ring ideal oo
      randomMonomialSet(n,D,p,Coefficients=>ZZ/101)
      ring ideal oo
  SeeAlso
    randomMonomialSet
    randomMonomialSets
    randomMonomialIdeals
///

doc ///
  Key
    VariableName
    [randomMonomialSet, VariableName]
    [randomMonomialSets, VariableName]
    [randomMonomialIdeals, VariableName]
  Headline
    optional input to choose the variable name for the generated polynomials
  Description
    Text
      Put {\tt VariableName => x} for a choice of string or symbol x as an argument in
      the function @TO randomMonomialSet@ or @TO randomMonomialSets@
    Example 
      n=2; D=3; p=0.2;
      randomMonomialSet(n,D,p)
      randomMonomialSet(n,D,p,VariableName => y)
///

doc ///
  Key
    [randomMonomialSet, Strategy]
    [randomMonomialSets, Strategy]
  Headline
    optional input to choose the strategy for generating the monomial set
  Description
    Text
      Put {\tt Strategy => "ER"} or {\tt Strategy => "Minimal"} as an argument in the function @TO randomMonomialSet@ or @TO randomMonomialSets@. 
      "ER" draws random sets of monomials from the ER-type distribution B(n,D,p), while "Minimal" saves computation time by using quotient rings to exclude any non-minimal generators from the list.
///

doc ///
 Key
   IncludeZeroIdeals
   [idealsFromGeneratingSets, IncludeZeroIdeals]
   [randomMonomialIdeals, IncludeZeroIdeals]
   [bettiStats, IncludeZeroIdeals]
 Headline
   optional input to choose whether or not zero ideals should be included
 Description
   Text
     When the option is used with the method @TO randomMonomialIdeals@, if {\tt IncludeZeroIdeals => true} (the default), then zero ideals will be included in the list of random monomial ideals. 
     If {\tt IncludeZeroIdeals => false}, then any zero ideals produced will be excluded, along with the number of them. 
   Example
     n=2;D=2;p=0.0;N=1;
     ideals = randomMonomialIdeals(n,D,p,N)
   Text
     The 0 listed is the zero ideal: 
   Example
     ideals_0
   Text
     In the example below, in contrast, the list of ideals returned is empty since the single zero ideal generated is excluded:
   Example
     randomMonomialIdeals(n,D,p,N,IncludeZeroIdeals=>false)
   Text
     When the option is used with the method @TO bettiStats@, ....
     **** Need to move here the doc for [bettiStats, IncludeZeroIdeals] from the  bettiStats doc node. :) 
 SeeAlso
   randomMonomialIdeals
   bettiStats
///
doc ///
 Key
  dimStats
  (dimStats,List)
 Headline
  statistics on the Krull dimension of a list of monomialIdeals 
 Usage
  dimStats(List)
 
 Inputs
  listOfIdeals: List
    of @TO monomialIdeal@s
  
 Outputs
  : Sequence 
   whose first entry is the average Krull dimension of a list of monomialIdeals, the second entry is the standard deviation of the Krull dimension, and third entry (if option turned on) is the Krull dimension tally
 Description
  Text
   dimStats finds the average and standard deviaation of the Krull dimension for a list of monomialIdeals.   
  Example
    L=randomMonomialSet(3,3,1.0);
    R=ring(L#0);
    listOfIdeals = {monomialIdeal(R_0^3,R_1,R_2^2), monomialIdeal(R_0^3, R_1, R_0*R_2)};
    dimStats(listOfIdeals)
  Text
   The following examples use the existing functions @TO randomMonomialSets@ and @TO idealsFromGeneratingSets@ or @TO randomMonomialIdeals@ to automatically generate a list of ideals, rather than creating the list manually:
  Example
   listOfIdeals = idealsFromGeneratingSets(randomMonomialSets(4,3,1.0,3));
   dimStats(listOfIdeals)
  Example
   listOfIdeals = randomMonomialIdeals(4,3,1.0,3);
   dimStats(listOfIdeals)
  Text
   Note that this function can be run with a list of any objects to which @TO dim@ can be applied. 
///


doc ///
 Key
   ShowTally
   [dimStats, ShowTally]
   [mingenStats, ShowTally]
   [degStats, ShowTally]
   [regStats, ShowTally]
   [pdimStats, ShowTally]
 Headline
   optional input to choose if the tally is to be returned 
 Description
   Text
     If {\tt ShowTally => false} (the default value), then only the 2 basic statistics - mean and standard deviation - of the function will be returned. 
     If {\tt ShowTally => true}, then both the statistics and the tally will be returned. 

   Example
     n=3;D=3;p=0.0;N=3;
     ideals = randomMonomialIdeals(n,D,p,N);
     dimStats(ideals)
     mingenStats(ideals)
     degStats(ideals)
     pdimStats(ideals)
   Text
     In the example above, only the statistics are outputted since by default {\tt ShowTally => false}. 
   Text
    In order to view the tally, ShowTally must be set to true ({\tt ShowTally => true}) when the function is called: 
   Example
     L=randomMonomialSet(3,3,1.0);
     -- QUESTION: what is the purpose of using L here?! 
     R=ring(L#0);
     listOfIdeals = {monomialIdeal(R_0^3,R_1,R_2^2), monomialIdeal(R_0^3, R_1, R_0*R_2)};
     dimStats(listOfIdeals,ShowTally=>true)
     mingenStats(listOfIdeals,ShowTally=>true)
     degStats(listOfIdeals,ShowTally=>true)
     regStats(listOfIdeals,ShowTally=>true)
     pdimStats(listOfIdeals,ShowTally=>true)
 SeeAlso
   dimStats
   mingenStats
   degStats
   regStats
///

doc ///
 Key
  pdimStats
  (pdimStats,List)
 Headline
  statistics on projective dimension of a list of monomial ideals
 Usage
  pdimStats(List)
 Inputs
  ideals: List
    of @TO monomialIdeal@s
 Outputs
  : Sequence 
   whose first entry is the mean projective dimension, the second entry is the standard deviation of the projective dimension, and third entry (if option turned on) is the projective dimension tally for quotient rings of ideals in the list {\tt ideals}.
 Description
  Text
   pdimStats finds the mean and standard deviation of the projective dimension of elements in the list: 
  Example
   R=ZZ/101[a,b,c]
   ideals = {monomialIdeal(a^3,b,c^2), monomialIdeal(a^3,b,a*c)}
   pdimStats(ideals)
  Text
   pdimStats will also output the projective dimension Tally using the optional input ShowTally
  Example
   R=ZZ/101[a,b,c]
   ideals = {monomialIdeal(a,c),monomialIdeal(b),monomialIdeal(a^2*b,b^2)}
   pdimStats(ideals, ShowTally=>true)
  Text
   The following example uses the existing function @TO randomMonomialIdeals@ to automatically generate a list of ideals, rather than creating the list manually:
  Example
   listOfIdeals = randomMonomialIdeals(4,3,1.0,3)
   pdimStats(listOfIdeals)
  Text
   Note that this function can be run with a list of @TO ideal@s as well. 
 SeeAlso
   ShowTally
///

doc ///
 Key
  regStats
  (regStats, List)
 Headline
  statistics on the regularities of a list of monomialIdeals
 Usage
  regStats(List)
 Inputs
  : List
   of @TO monomialIdeal@s
 Outputs
  : Sequence
   whose first entry is the mean regularity of a list of monomialIdeals, second entry is the standard deviation of the regularities, and third entry (if option is turned on) is the regularity tally.
 Description
  Text
   regStats removes zero ideals from the list of ideals, then calculates the average and the standard deviation of the regularity of the list of nonzero ideals.
  Example
   n=4; D=3; p={0.0,1.0,0.0}; N=2;
   B=randomMonomialIdeals(n,D,p,N)
   regStats(B)
  Text
   If the list given is a list of all zero ideals, regStats returns -infinity for the mean regularity.
  Example
   B=randomMonomialIdeals(3,3,0.0,1)
   regStats(B)
  Text
   Note that this function can be called on a list of @TO Ideal@ objects instead.
 Caveat
  regStats removes zero ideals from the list of ideals before computing the two values.
 SeeAlso
  ShowTally
 ///
 
doc ///
 Key
  CMStats
  (CMStats,List)
 Headline
  percentage of monomialIdeals in the given list whose quotient ring is Cohen-Macaulay
 Usage
  CMStats(List)
 Inputs
  ideals: List
    of @TO monomialIdeal@s
 Outputs
  : RR
   the percentage of Cohen-Macaulay ideals in the list
 Description
  Text
   CMStats simply checks whether the coordinate ring of each ideal in the given sample is arithmetically Cohen-Macaulay, and returns the percentage that are.
  Example
    L=randomMonomialSet(3,3,1.0);
    R=ring(L#0);
    ideals = {monomialIdeal(R_0^3,R_1,R_2^2), monomialIdeal(R_0^3, R_1, R_0*R_2)};
    CMStats(ideals)
  Text
    Note that the method can be run on a list of @TO Ideal@s, too.
///

doc ///
 Key
  borelFixedStats
  (borelFixedStats ,List)
 Headline
  percentage of Borel-fixed monomialIdeals in the given list
 Usage
  borelFixedStats(List)
 
 Inputs
  listOfIdeals: List
    of @TO monomialIdeal@s
 Outputs
  : RR
   the percentage of Borel-fixed monomialIdeals in the list
 Description
  Text
   borelFixedStats takes a list of monomialIdeals and returns the percentage of Borel-fixed ideals in the list of monomialIdeals as a real number  
  Example
    L=randomMonomialSet(3,3,1.0);
    R=ring(L#0);
    listOfIdeals = {monomialIdeal(R_0^3), monomialIdeal(R_0^3, R_1, R_0*R_2)};
    borelFixedStats(listOfIdeals)
///

--******************************************--
-- TESTS     	     	       	    	    -- 
--******************************************--

--************************--
--  randomMonomialSets  --
--************************--

TEST ///
    -- Check there are N samples
    N=10;
    n=3; D=2; p=0.5;
    assert (N==#randomMonomialSets(n,D,p,N))
    N=13;
    n=5; D=3; p={0.5,0.25,0.3};
    assert (N==#randomMonomialSets(n,D,p,N))
    N=10;
    n=3; D=2; M=10;
    assert (N==#randomMonomialSets(n,D,M,N))
    N=7;
    n=4; D=3; M={3,3,3};
    assert (N==#randomMonomialSets(n,D,M,N))
///

TEST ///
    -- Check multiple samples agree
    n=4; D=3;
    L = randomMonomialSets(n,D,1.0,3);
    assert (set L#0===set L#1)
    assert (set L#0===set L#2)
    
///

TEST ///
    --Check monomials are in the same ring
    n = 4; D = 3;
    L = randomMonomialSets(n,D,1.0,3);
    assert(ring(L#0#0)===ring(L#1#0))
    assert(ring(L#1#1)===ring(L#1#2))
    assert(ring(L#2#0)===ring(L#1#2))
///

--***********************--
--  randomMonomialSet  --
--***********************--

TEST ///
    --Check monomials are in the same ring
    n = 4; D = 3;
    L = randomMonomialSet(n,D,1.0);
    assert(ring(L#0)===ring(L#1))
    assert(ring(L#2)===ring(L#3))
///

TEST ///
    -- Check no terms are chosen for a probability of 0

    assert (0==(randomMonomialSet(5,5,0.0))#0)
    assert (0==(randomMonomialSet(5,4,toList(4:0.0)))#0)
    assert (0==(randomMonomialSet(5,4,0.0, Strategy=>"Minimal"))#0)
    assert (0==(randomMonomialSet(5,4,toList(4:0.0), Strategy=>"Minimal"))#0)
    assert (0==(randomMonomialSet(5,4,0))#0)
    assert (0==(randomMonomialSet(5,4,toList(4:0)))#0)
///

TEST ///
    -- Check all possible values are outputted with a probability of 1
    n=4; D=3;
    assert (product(toList((D+1)..D+n))/n!-1==#randomMonomialSet(n,D,1.0))
    assert (product(toList((D+1)..D+n))/n!-1==#randomMonomialSet(n,D,{1.0,1.0,1.0}))
    n=6; D=2;
    assert (product(toList((D+1)..D+n))/n!-1==#randomMonomialSet(n,D,1.0))
    assert (product(toList((D+1)..D+n))/n!-1==#randomMonomialSet(n,D,{1.0,1.0}))
    n=4;D=5;
    assert (# flatten entries basis (1, QQ[x_1..x_n])==#randomMonomialSet(n,D,1.0, Strategy=>"Minimal"))
    assert (# flatten entries basis (2, QQ[x_1..x_n])==#randomMonomialSet(n,D,{0.0,1.0,1.0,1.0,1.0}, Strategy=>"Minimal"))
    assert (# flatten entries basis (3, QQ[x_1..x_n])==#randomMonomialSet(n,D,{0.0,0.0,1.0,1.0,1.0}, Strategy=>"Minimal"))
    assert (# flatten entries basis (4, QQ[x_1..x_n])==#randomMonomialSet(n,D,{0.0,0.0,0.0,1.0,1.0}, Strategy=>"Minimal"))
    assert (# flatten entries basis (5, QQ[x_1..x_n])==#randomMonomialSet(n,D,{0.0,0.0,0.0,0.0,1.0}, Strategy=>"Minimal"))
///

TEST ///
    -- Check every monomial is generated
    L=randomMonomialSet(2,3,1.0)
    R=ring(L#0)
    assert(set L===set {R_0,R_1,R_0^2,R_0*R_1,R_1^2,R_0^3,R_0^2*R_1,R_0*R_1^2,R_1^3})
    L=randomMonomialSet(2,3,9)
    R=ring(L#0)
    assert(set L===set {R_0,R_1,R_0^2,R_0*R_1,R_1^2,R_0^3,R_0^2*R_1,R_0*R_1^2,R_1^3})
    L=randomMonomialSet(3,3,{0.0,1.0,0.0})
    R=ring(L#0)
    assert(set L===set {R_0^2,R_0*R_1,R_1^2,R_0*R_2,R_1*R_2,R_2^2})
    L=randomMonomialSet(3,3,1.0, Strategy=>"Minimal");
    R=ring(L#0);
    assert(set L===set {R_0, R_1, R_2})
    L=randomMonomialSet(2,3,{2,3,4})
    R=ring(L#0)
    assert(set L===set {R_0,R_1,R_0^2,R_0*R_1,R_1^2,R_0^3,R_0^2*R_1,R_0*R_1^2,R_1^3})
    L=randomMonomialSet(3,3,{0.0,1.0,1.0}, Strategy=>"Minimal");
    R=ring(L#0);
    assert(set L===set {R_0^2,R_0*R_1,R_1^2,R_0*R_2,R_1*R_2,R_2^2})
    L=randomMonomialSet(3,3,{0.0,0.0,1.0}, Strategy=>"Minimal");
    R=ring(L#0);
    assert(set L===set {R_0^3,R_0^2*R_1,R_0^2*R_2,R_0*R_1*R_2,R_1^3,R_0*R_1^2,R_1^2*R_2,R_0*R_2^2,R_1*R_2^2,R_2^3})
///

TEST ///
    -- Check max degree of monomial less than or equal to D
    n=10; D=5;
    assert(D==max(apply(randomMonomialSet(n,D,1.0),m->first degree m)))
    assert(D==max(apply(randomMonomialSet(n,D,toList(D:1.0)),m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(D==max(apply(randomMonomialSet(n,D,M),m->first degree m)))
    assert(D==max(apply((randomMonomialSet(n,D,{0.0,0.0,0.0,0.0,1.0}, Strategy=>"Minimal"),m->first degree m))))
    n=4; D=7;
    assert(D==max(apply(randomMonomialSet(n,D,1.0),m->first degree m)))
    assert(D==max(apply(randomMonomialSet(n,D,toList(D:1.0)),m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(D==max(apply(randomMonomialSet(n,D,M),m->first degree m)))
    assert(D==max(apply(randomMonomialSet(n,D,toList(D:1)), m->first degree m)))
///

TEST ///
    -- Check min degree of monomial greater than or equal to 1
    n=8; D=6;
    assert(1==min(apply(randomMonomialSet(n,D,1.0),m->first degree m)))
    assert(1==min(apply(randomMonomialSet(n,D,toList(D:1.0)),m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(1==min(apply(randomMonomialSet(n,D,M),m->first degree m)))
    n=3; D=5;
    assert(1==min(apply(randomMonomialSet(n,D,1.0),m->first degree m)))
    assert(1==min(apply(randomMonomialSet(n,D,toList(D:1.0)),m->first degree m)))
    M=lift(product(toList((D+1)..(D+n)))/n!-1,ZZ);
    assert(1==min(apply(randomMonomialSet(n,D,M),m->first degree m)))
    n=10; D=5;
    assert(1==min(apply((randomMonomialSet(n,D,1.0, Strategy=>"Minimal"),m->first degree m))))
    assert(1==min(apply((randomMonomialSet(n,D,toList(D:1.0), Strategy=>"Minimal"),m->first degree m))))
    assert(1==min(apply(randomMonomialSet(n,D,toList(D:1)), m->first degree m)))
///


--*************************--
--  bettiStats  --
--*************************--
TEST///
   R = ZZ/101[a..c];
   L={monomialIdeal (a^2*b,b*c), monomialIdeal(a*b,b*c^3)};
   (meanBettiShape,meanBetti) = bettiStats L;
   -- mean Betti table:
   b=new BettiTally from { (0,{0},0) => 2, (1,{2},2) => 2, (1,{3},3) => 1, (2,{4},4) => 1, (1,{4},4) => 1, (2,{5},5) =>1 }
   assert(1/2*sub(matrix lift(2*meanBetti,ZZ),RR) ==  1/2*sub(matrix b,RR))
   -- mean Betti shape: 
   b=new BettiTally from { (0,{0},0) => 1, (1,{2},2) => 1, (1,{3},3) => 0.5, (2,{4},4) => 0.5, (1,{4},4) => 0.5, (2,{5},5) =>0.5 }
   assert(1/2*sub(matrix lift(2*meanBettiShape,ZZ),RR) ==  1/2*sub(matrix lift(2*b,ZZ),RR))
///


--*************************--
--  degStats  --
--*************************--
TEST///
   --check for p = 0 the average degree should be 1
   listOfIdeals = idealsFromGeneratingSets(randomMonomialSets(3,4,0.0,6));
   assert(1==(degStats(listOfIdeals))_0)
   assert(0==(degStats(listOfIdeals))_1)
   listOfIdeals = idealsFromGeneratingSets(randomMonomialSets(7,2,0,3));
   assert(1==(degStats(listOfIdeals))_0)
   assert(0==(degStats(listOfIdeals))_1)
   --check for p = 1 the average degree is 1
   listOfIdeals = idealsFromGeneratingSets(randomMonomialSets(3,4,1.0,6));
   assert(1==(degStats(listOfIdeals))_0)
   assert(0==(degStats(listOfIdeals))_1)
   --Check average is correct for set monomials
   L=randomMonomialSet(3,3,1.0);
   R=ring(L#0);
   listOfIdeals={monomialIdeal(R_0^3,R_1,R_2^2),monomialIdeal(R_0^3,R_1,R_0*R_2)};
   assert(3.5==(degStats(listOfIdeals,ShowTally=>true))_0)
   assert(2.5==(degStats(listOfIdeals,ShowTally=>true))_1)
   assert(2==sum(values(degStats(listOfIdeals, ShowTally=>true))_2))
   listOfIdeals={monomialIdeal(0_R),monomialIdeal(R_2^2)};
   assert(1.5==(degStats(listOfIdeals, ShowTally=>true))_0)
   assert(0.5==(degStats(listOfIdeals, ShowTally=>true))_1)
   assert(2==sum(values(degStats(listOfIdeals,ShowTally=>true))_2))
   listOfIdeals={monomialIdeal(R_0),monomialIdeal(R_0^2*R_2),monomialIdeal(R_0*R_1^2,R_1^3,R_1*R_2,R_0*R_2^2)};
   assert(sub(8/3,RR)==(degStats(listOfIdeals,ShowTally=>true))_0)
   assert((sub(14/9,RR))^(1/2)==(degStats(listOfIdeals,ShowTally=>true))_1)
   assert(3==sum(values(degStats(listOfIdeals,ShowTally=>true))_2))
 
///

--************************--
--  dimStats  --
--************************--
TEST ///
    --check for p = 0 the average krull dimension is n
    listOfIdeals = idealsFromGeneratingSets(randomMonomialSets(3,4,0.0,6));
    assert(3==(dimStats(listOfIdeals))_0)
    assert(0==(dimStats(listOfIdeals))_1)
    listOfIdeals = idealsFromGeneratingSets(randomMonomialSets(7,2,0,3));
    assert(7==(dimStats(listOfIdeals))_0)
    assert(0==(dimStats(listOfIdeals))_1)
    --check for p = 1 the average krull dimension is 0
    listOfIdeals = idealsFromGeneratingSets(randomMonomialSets(3,4,1.0,6));
    assert(0==(dimStats(listOfIdeals))_0)
    assert(0==(dimStats(listOfIdeals))_1)
    --check for set monomials
    L=randomMonomialSet(3,3,1.0);
    R=ring(L#0);
    listOfIdeals = {monomialIdeal(R_0^3,R_1,R_2^2), monomialIdeal(R_0^3, R_1, R_0*R_2)};
    assert(.5==(dimStats(listOfIdeals, ShowTally=>true))_0)
    assert(.5==(dimStats(listOfIdeals, ShowTally=>true))_1)
    assert(2==sum( values (dimStats(listOfIdeals, ShowTally=>true))_2))
    listOfIdeals = {monomialIdeal (0_R), monomialIdeal (R_2^2)};
    assert(2.5== (dimStats(listOfIdeals,ShowTally=>true))_0)
    assert(.5==(dimStats(listOfIdeals, ShowTally=>true))_1)
    assert(2==sum( values (dimStats(listOfIdeals, ShowTally=>true))_2))
    listOfIdeals = {monomialIdeal (R_0), monomialIdeal (R_0^2*R_2), monomialIdeal(R_0*R_1^2,R_1^3,R_1*R_2,R_0*R_2^2)};
    assert(sub(5/3,RR)==(dimStats(listOfIdeals,ShowTally=>true))_0)
    assert((3-(sub(5/3,RR))^2)^(1/2)==(dimStats(listOfIdeals,ShowTally=>true))_1)
    assert(3==sum( values (dimStats(listOfIdeals, ShowTally=>true))_2))
///

--************************--
--  randomMonomialIdeals  --
--************************--

TEST ///
  -- check the number of ideals
  n=5; D=5; p=.6; N=3;
  B = flatten randomMonomialIdeals(n,D,p,N,IncludeZeroIdeals=>false);
  assert (N===(#B-1+last(B))) -- B will be a sequence of nonzero ideals and the number of zero ideals in entry last(B)
  C = randomMonomialIdeals(n,D,p,N,IncludeZeroIdeals=>true);
  assert (N===#C)
///

TEST ///
  -- check the number of monomials in the generating set of the ideal
  n=4; D=6; M=7; N=1;
  B = flatten randomMonomialIdeals(n,D,M,N);
  assert (M>=numgens B_0)
///

--************--
--  regStats  --
--************--
TEST ///
  -- check average regularity
  n=3; D=5; N=4; p=1.0;
  B=randomMonomialIdeals(n,D,p,N);
  assert((1,0)==regStats(B))
  p={0,1,0,0,0};
  B=randomMonomialIdeals(n,D,p,N);
  assert((2,0)==regStats(B))
  p={0,0,1,0,0};
  B=randomMonomialIdeals(n,D,p,N);
  assert((3,0)==regStats(B))
  p={0,0,0,1,0};
  B=randomMonomialIdeals(n,D,p,N);
  assert((4,0)==regStats(B))
  p={0,0,0,0,1};
  B=randomMonomialIdeals(n,D,p,N);
  assert((5,0)==regStats(B))
  p=0;
  B=randomMonomialIdeals(n,D,p,N);
  assert((-infinity,0)==regStats(B))
///

TEST ///
  -- check all stats
  L=randomMonomialSet(3,3,1.0);
  R=ring(L#0);
  listOfIdeals={monomialIdeal(R_1,R_2^2),monomialIdeal(R_0^3,R_1,R_0*R_2)};
  A = regStats(listOfIdeals, ShowTally=>true);
  assert(2.5===A_0)
  assert(0.5===A_1)
  assert(2==sum(values(A_2)))
///
--***************--
--    CMStats    --
--***************--

TEST ///
 L=randomMonomialSet(5,1,1.0); R=ring(L#0);
 listOfIdeals = {monomialIdeal(0_R)};
 assert(1==CMStats(listOfIdeals))
 listOfIdeals = {monomialIdeal(R_0*R_1, R_2*R_0)};
 assert(0==CMStats(listOfIdeals))
 listOfIdeals = {monomialIdeal(0_R), monomialIdeal(R_0*R_1, R_2*R_0)};
 assert(.5==CMStats(listOfIdeals))
 listOfIdeals = {monomialIdeal(0_R), monomialIdeal(R_0*R_1, R_2*R_0), monomialIdeal(R_0)};
 assert(sub(2/3,RR)==CMStats(listOfIdeals))
///

--********************--
--  borrelFixedStats  --
--********************--

TEST ///
L=randomMonomialSet(5,1,1.0); R=ring(L#0);
listOfIdeals = {monomialIdeal(0_R)};
assert(1==borelFixedStats(listOfIdeals))
listOfIdeals = {monomialIdeal(R_0*R_1)};
assert(0==borelFixedStats(listOfIdeals))
listOfIdeals = {monomialIdeal(R_0), monomialIdeal(R_0*R_1)};
assert(.5==borelFixedStats(listOfIdeals))
listOfIdeals = {monomialIdeal(0_R), monomialIdeal(R_0*R_1, R_2*R_0), monomialIdeal(R_0)};
assert(sub(2/3,RR)==borelFixedStats(listOfIdeals))
///

--***************--
--  mingenStats  --
--***************--

TEST ///
  -- check average number of minimum generators
  n=4; D=3; p=1.0; N=3;
  B = randomMonomialIdeals(n,D,p,N);
  C = mingenStats(B);
  assert (sub(n,RR)===C_0)
  assert (sub(0,RR)===C_1)
  p={0.0,1.0,0.0};
  D = randomMonomialIdeals(n,D,p,N);
  E = mingenStats(D);
  assert (sub(10,RR)===E_0)
  assert (sub(0,RR)===E_1)
///

TEST ///
  -- check average degree complexity
  n=3; D=5; p=1.0; N=5;
  B = randomMonomialIdeals(n,D,p,N);
  C = mingenStats(B);
  assert(sub(1,RR)===C_2)
  assert(sub(0,RR)===C_3)
  p={0.0,0.0,0.0,0.0,1.0};
  D = randomMonomialIdeals(n,D,p,N);
  E = mingenStats(D);
  assert(sub(5,RR)===E_2)
  assert(sub(0,RR)===E_3)
///

TEST ///
  L=randomMonomialSet(3,3,1.0);
  R=ring(L#0);
  listOfIdeals={monomialIdeal(R_1,R_2^2),monomialIdeal(R_0^3,R_1,R_0*R_2)};
  A = mingenStats(listOfIdeals, ShowTally=>true);
  assert(2.5===A_0)
  assert(0.5===A_1)
  assert(2==sum(values(A_2)))
  assert(2.5===A_3)
  assert(0.5===A_4)
  assert(2==sum(values(A_5)))
///

--***************--
--   pdimStats   --
--***************--

TEST ///
  L=randomMonomialSet(3,3,1.0);
  R=ring(L#0);
  listOfIdeals={monomialIdeal(0_R)};
  assert(sub(0,RR)==(pdimStats(listOfIdeals))_0)
  assert(sub(0,RR)==(pdimStats(listOfIdeals))_1)
  listOfIdeals={monomialIdeal(R_0,R_1,R_2)};
  assert(sub(3,RR)==(pdimStats(listOfIdeals))_0)
  assert(sub(0,RR)==(pdimStats(listOfIdeals))_1)
  listOfIdeals={monomialIdeal(0_R),monomialIdeal(R_0*R_1^2,R_1^3,R_2)};
  assert(sub(1.5,RR)==(pdimStats(listOfIdeals))_0)
  assert(sub(1.5,RR)==(pdimStats(listOfIdeals))_1)
  listOfIdeals={monomialIdeal(R_0^2*R_1,R_2)};
  assert(sub(2,RR)==(pdimStats(listOfIdeals))_0)
  assert(sub(0,RR)==(pdimStats(listOfIdeals))_1)
  listOfIdeals={monomialIdeal(R_0,R_2),monomialIdeal(0_R),monomialIdeal(R_0^2*R_1,R_1^2)};
  assert(sub(4/3,RR)==(pdimStats(listOfIdeals))_0)
  assert(sub(((8/3)-(16/9))^(1/2),RR)==(pdimStats(listOfIdeals))_1)
///

end

--****************************--
--  idealsFromGeneratingSets  --
--****************************--

TEST///
  -- check the number of ideals
  n=5; D=5; p=.6; N=3;
  B = flatten idealsFromGeneratingSets(randomMonomialSets(n,D,p,N),IncludeZeroIdeals=>false);
  assert (N===(#B-1+last(B))) -- B will be a sequence of nonzero ideals and the number of zero ideals in entry last(B)
  C = idealsFromGeneratingSets(randomMonomialSets(n,D,p,N),IncludeZeroIdeals=>true);
  assert (N===#C)
///

TEST ///
  --check that all elements are MonomialIdeal
  n=3;D=3;p=1.0;N=5;
  B=idealsFromGeneratingSets(randomMonomialSets(n,D,p,N));
  assert (all(B,b->instance(b,MonomialIdeal)))
  C=idealsFromGeneratingSets(randomMonomialSets(n,D,p,N),IncludeZeroIdeals=>false);
  assert (all(C#0,c->instance(c,MonomialIdeal)))
///

end

restart;
uninstallPackage"RandomMonomialIdeals";
installPackage"RandomMonomialIdeals";
viewHelp bettiStats

check RandomMonomialIdeals 
viewHelp RandomMonomialIdeals
