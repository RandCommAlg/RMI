-------------------------------------------------------------------------------------
-- THIS FILE IS LOADED BY script.m2 -- 
-------------------------------------------------------------------------------------

-- n (num vars) is predefined in script.m2
-- D (max deg) is predefined in script.m2
-- N (sample size) is predefined in script.m2
-- p (probability for ER model) is predefined in sript.m2

-------------------------------------------------------------------------------------
-- Random Monomial Ideals 
-- Sept. 2016
-- Dane Wilburne, edit by SP Oct.2016

-- Creating and analyzing samples of ER-type random monomial ideals
-------------------------------------------------------------------------------------

 load("RMIs.m2")  -- get all the functions loaded ! 
 
-------------------------------------------------------------------------------------
-------------------------------------Analyze RMIs ----------------------------
-------------------------------------------------------------------------------------

R = QQ[x_1..x_n] 
allMonomials = drop(flatten flatten apply(D+1,d->entries basis(d,R)),1);


-- go through list allMonomials, and for each monomial m in the list, select a number in Unif(0,1); 
-- if that number <= p, then include the monomial m in a generating set: 
--B=select(allMonomials, m-> random(0.0,1.0)<=p) 
-- since iid~Unif(0,1), this is same as keeping each possible monomial w/ probability p. 
--In addition, we need a sample of size N of sets of monomials like these, so here we go: 
B=apply(N,i-> 
    select(allMonomials, m-> random(0.0,1.0)<=p) 
    );
--the result:
-- B = list of random monomial ideal generating sets. 

--make a RMI from each random gen set, store as list called `ideals', and Z is number of zero ideals:
(ideals,Z) = makeRandIdeals(B,p,D,basefilename); 


--compute everything you ever wanted to know about your sample of RMIs: 
if N==Z then (
    print("The model only generated zero ideals in this sample! Stopping simulation analysis now.") 
    ) else randomIdealData(ideals,N,Z,p,n,D,basefilename) 



