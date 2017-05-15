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
	    }
	    {
	    	Name => "Monica Yun", 
	    	Email => "myun1@hawk.iit.edu", 
	    	HomePage => ""
	    }
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
    "firstFunction",
    "randomGeneratingSets"
    }

--***************************************--
--  Function provided by FirstPackage.m2 --
--***************************************--

firstFunction = method(TypicalValue => String)
firstFunction ZZ := String => n -> if n == 1 then "Hello World!" else "D'oh!"

--****************************************************************--
--  Methods written by D.Wilburne and S.Petrovic for RMI, 2016-17 --
--****************************************************************--


--**********************************--
--  Methods that need documentation --
--**********************************--
randomGeneratingSets = method(TypicalValue => List)
randomGeneratingSets (ZZ,ZZ,RR,ZZ) := List =>  (n,D,p,N) -> (
    x :=symbol x;
    R := QQ[x_1..x_n];
    allMonomials := flatten flatten apply(1..D,d->entries basis(d,R));
    -- this generates a list of all possible monomials of degree <= D in n variables
    -- go through list allMonomials, and for each monomial m in the list, select a number in Unif(0,1);
    -- if that number <= p, then include the monomial m in a generating set B
    -- since iid~Unif(0,1), this is same as keeping each possible monomial w/ probability p.
    -- we need a sample of size N of sets of monomials like these, so we repeat this process N times:
    B := apply(N,i-> select(allMonomials, m-> random(0.0,1.0)<=p) )
    --the result:
    -- B = list of random monomial ideal generating sets.
)

--**********************************--
--  Methods that need reformatting  --
--**********************************--

-- TO BE ADDED BY S.P. 


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
  Caveat
   Still trying to figure this out. [REMOVE ME]
///

doc ///
 Key
  (firstFunction,ZZ)
  firstFunction
  Headline
   a silly first function
  Usage
   firstFunction n
  Inputs
   n:
  Outputs
   :
   a silly string, depending on the value of {\tt n}
  Description
   Text
    Here we show an example.
   Example
    firstFunction 1
    firstFunction 0
 ///
 
doc ///
 Key
  (randomGeneratingSets,ZZ,ZZ,RR,ZZ)
  randomGeneratingSets
 Headline
  randomly generates lists of monomials, up to a given dimension
 Usage
  randomGeneratingSets (ZZ,ZZ,RR,ZZ)
 Inputs
  n: ZZ
  D: ZZ
  p: RR
  N: ZZ
 Outputs
  :
  a list of generating sets of monomials
 Description
  Text
   given number of variables {n}, degree bound {D}, and probability parameter {p}, {N} generating sets of monomials are randomly generated
  Example
   randomGeneratingSets(2,3,0.2,10)
   randomGeneratingSets(3,4,1,1)
   randomGeneratingSets(5,2,0.6,4)
  Text
   Explain some more.
 SeeAlso
  firstFunction
///


--******************************************--
-- TESTS     	     	       	    	    -- 
--******************************************--


TEST ///
    assert ( firstFunction 2 == "D'oh!" )
///

end

You can write anything you want down here.  I like to keep examples
as I’m developing here.  Clean it up before submitting for
publication.  If you don't want to do that, you can omit the "end"
above.

