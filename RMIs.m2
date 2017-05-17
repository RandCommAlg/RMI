-->> THE BASIC METHODS FOR PACKAGE RandomMonomialIdeals.m2 <<-- 

-------------------------------------------------------------------------------------
-- THIS FILE IS LOADED BY RMIscript.m2 -- 
-------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------
-- Random Monomial Ideals (M2 Side)
-- Sept. 2016
-- Dane Wilburne
-- SP edited 18.9.2016. to save different file name extentsions encoding values of n,p,D,N for each run.
-- SP edited 8.10.2016. to add testEisenbudGoto function.
-- SP edited 28.11.2016. to make shorter simulations for the paper (2 functions ending with "ForPaper") so we can try to get those 2 tables...

--This is the M2 side of the code for creating and analyzing samples of ER-type random monomial ideals

--This code defines a set of functions that, from   N randomly generated sets of monomials (see RMIscript.m2), does 2 things:
--1) computes a set of invariants (see below) for the random monomial ideal corresponding to each set, 
--2) outputs a collection of text files with the computed invariants, and prints samples averages/summaries.
--The invariants computed by this code are: regularity, projective dimension, Krull dimension, degree, number of minimal generators, degree complexity, Cohen-Macaulayness, Borel-fixedness, and Betti numbers
-------------------------------------------------------------------------------------


-------------------------------------------------------------------------------------
------------------------------ Load gens sets, create ideals ------------------------
-------------------------------------------------------------------------------------

-- Here is a method that we may want to use if we end up re-reading some of the simulated ideals:
--reads .txt file of gen sets into M2
readBs = filename -> lines(get(filename))


--creates ideals from random gen sets
--saves them to file `randomIdeals' - with an extension encoding values of n,p,D,N. 
--see RMIscript.m2 for short script to run this function
--makeRandIdeals = (B,p,D,basefilename) -> ( 
idealsFromGeneratingSets =  method(TypicalValue => List, Options => {includeZeroIdeals => false})
idealsFromGeneratingSets (List,RR,ZZ,String) := o -> (B,p,D,basefilename) -> (
    --inputs:
	--B = list of random generating sets
	--p (probability for ER model) is predefined in sript.m2 
	--D (max deg) is predefined in script.m2
    fileNameExt = concatenate("_for_params_n",toString(n),"_p",toString(p),"_D",toString(D),"_N",toString(N));
    filename=concatenate(basefilename,"randomIdeals",fileNameExt,".txt");
    for i from 0 to #B-1 do {
	ideals = B / (b-> ideal b);
	filename << toString B#i << endl; 
    	if not o.includeZeroIdeals then (nonzeroIdeals,numberOfZeroIdeals) = extractNonzeroIdeals(ideals); 
	};
    filename<<close;
    print(concatenate("There are ", toString(#B)," ideals in this sample."));
    print(concatenate("Of those, ", toString numZeroIdeals, " were the zero ideal."));
    if o.includeZeroIdeals then return ideals else return (nonzeroIdeals,numberOfZeroIdeals); 
)


-- Internal method that takes as input list of ideals and splits out the zero ideals, counting them:
    -- input list of ideals 
    -- output a sequence (list of non-zero ideals from the list , the number of zero ideals in the list)
-- (not exported, therefore no need to document) 
extractNonzeroIdeals = ( ideals ) -> (
    nonzeroIdeals = select(ideals,i->i != 0);
    numberOfZeroIdeals = # ideals - # nonzeroIdeals;
    -- numberOfZeroIdeals = # positions(B,b-> b#0==0); -- sinze 0 is only included if the ideal = ideal{}, this is safe too
    return(nonzeroIdeals,numberOfZeroIdeals)
    )


-- we may not need the next one for any of the methods in this file; we'll be able to determine this soon. keep for now.
-- Internal method that takes as input list of generating sets and splits out the zero ideals, counting them:
    -- input list of generating sets
    -- output a sequence (list of non-zero ideals from the list , the number of zero ideals in the list)
-- (not exported, therefore no need to document) 
extractNonzeroIdealsFromGens = ( generatingSets ) -> (
    nonzeroIdeals = select(generatingSets,i-> i#0 != 0_(ring i#0)); --ideal(0)*ring(i));
    numberOfZeroIdeals = # generatingSets - # nonzeroIdeals;
    -- numberOfZeroIdeals = # positions(B,b-> b#0==0); -- sinze 0 is only included if the ideal = ideal{}, this is safe too
    return(nonzeroIdeals,numberOfZeroIdeals)
    )


-------------------------------------------------------------------------------------
------------------------------Combine all of the functions below to one function-----
-------------------------------------------------------------------------------------

--This method calls each of the methods defined below and produces all of the associated output
--see RMIscript.m2 for short script to run this function
randomIdealData =     (ideals,N,Z,p,n,D,basefilename) -> (
--randomIdealData = method()--method can't take more than 4 arguments, what do we do, m2!! 
--randomIdealData (List,ZZ,ZZ,RR,ZZ,ZZ,String) :=     (ideals,N,Z,p,n,D,basefilename) -> (
        --inputs:
	--list of ideals
	--N, total sample size 
	--Z, number of zero ideals in sample 
	--p, ER model probability
	--n, num vars
	--D, max degree for ER model
    fileNameExt = concatenate("_for_params_n",toString(n),"_p",toString(p),"_D",toString(D),"_N",toString(N));
    print concatenate("All files stored in filenames with extension that describe parameters as follows: ", fileNameExt,"."); 
    averageRegularity = avgReg(ideals,N,basefilename,fileNameExt); 
    averagePdimension = avgPdim(ideals,N,basefilename,fileNameExt); 
    averageKrullDimension = avgDim(ideals,N,Z,basefilename,fileNameExt); 
    averageDegree = avgDeg(ideals, N,basefilename,fileNameExt); 
    (averageNumGens,averageDegreeComplexity) = avgMinGens(ideals,N,basefilename,fileNameExt);
    percentCM = cohMac(ideals,N,Z);
    percentBorel = borelFixed(ideals,N,Z);
    (bAvg, bAvgShape) = bettis(ideals,N,Z,basefilename,fileNameExt);
    summary(N,averageRegularity,averagePdimension,
	averageKrullDimension,averageDegree, averageNumGens,
	averageDegreeComplexity,percentCM,percentBorel,bAvg,
	bAvgShape,basefilename, fileNameExt)
    )

--This method calls SELECTED ONES of the methods defined below and produces all of the associated output
--see RMIscript.m2 for short script to run this function
randomIdealDataForPaper =     (ideals,N,Z,p,n,D,basefilename) -> (
--randomIdealData = method()--method can't take more than 4 arguments, what do we do, m2!! 
--randomIdealData (List,ZZ,ZZ,RR,ZZ,ZZ,String) :=     (ideals,N,Z,p,n,D,basefilename) -> (
        --inputs:
	--list of ideals
	--N, total sample size 
	--Z, number of zero ideals in sample 
	--p, ER model probability
	--n, num vars
	--D, max degree for ER model
    fileNameExt = concatenate("_for_params_n",toString(n),"_p",toString(p),"_D",toString(D),"_N",toString(N));
    print concatenate("All files stored in filenames with extension that describe parameters as follows: ", fileNameExt,"."); 
    averageRegularity = avgReg(ideals,N,basefilename,fileNameExt); 
    averagePdimension = avgPdim(ideals,N,basefilename,fileNameExt); 
    percentCM = cohMac(ideals,N,Z);
    summaryForPaper(N,averageRegularity,averagePdimension,percentCM,basefilename, fileNameExt)
    )




-------------------------------------------------------------------------------------
--------------------------------------Betti #s---------------------------------------
-------------------------------------------------------------------------------------

--computes bettis of each ideal, saves to file `bettis'  - with an extension encoding values of n,p,D,N. 
--computes avg. betti table, prints and saves to file `avgBetti' - with an extension encoding values of n,p,D,N. 
--returns avg betti table and avg betti table shape.
bettis =    (ideals,N,Z,basefilename,fileNameExt) -> (
    --
    -- NOTE --- THIS (AND OTHERS?) METHOD BREAKS IF THERE ARE NO ZERO IDEALS GENERATED. (!!) 
    -- THIS IS A SIMPLE BUG.
    -- FIX IT.
    -- 
  -- sum of the betti tables and betti shapes: 
    beta = new BettiTally; 
    betaShapes = new BettiTally;
    filename1 = concatenate(basefilename,"bettis",fileNameExt);
    -- add up all the betti tables: 
    apply(#ideals,i->( 
        resi = betti res ideals_i;
        filename1 << resi << endl;
        beta = beta + resi;
  	-- let's only keep a 1 in all spots where ther was a non-zero Betti number: 
	beta1mtx = matrix(resi);
	Rtemp=(ring ideals_i)^1/ideals_i;
	beta1shape = new BettiTally from mat2betti  matrix pack(1+pdim(Rtemp), apply(flatten entries beta1mtx, i-> if i>0 then i=1 else i=0));
	betaShapes = betaShapes + beta1shape
	)
    );
    -- compute the average Betti table:
    b = mat2betti(1/#ideals*(sub(matrix(beta), RR)));
    -- compute the average Betti table shape: 
    bShape = mat2betti(1/#ideals*(sub(matrix(betaShapes), RR)));
    -- Now, add the betti tables of all the zero ideals:
    if Z>0 then ( 
	betaWithZeroIdeals=new BettiTally from beta;
	betaShapesWithZeroIdeals=new BettiTally from betaShapes;
	apply(Z,i-> ( -- note that Z = N - #ideals -- 
    	    -- compute the average Betti table:
	    betaWithZeroIdeals=betaWithZeroIdeals + betti res (ideal(0)*ring ideals_0);
    	    -- compute the average Betti table shape: 
	    betaShapeWithZeroIdeals=betaShapesWithZeroIdeals + betti res (ideal(0)*ring ideals_0);
	    )    
    	);
    -- compute the average Betti table including the count of zero ideals:
    bWithZeroIdeals = mat2betti(1/N*(sub(matrix(betaWithZeroIdeals), RR)));
    -- compute the average Betti table SHAPE including the count of zero ideals:
    bShapeWithZeroIdeals = mat2betti(1/N*(sub(matrix(betaShapeWithZeroIdeals), RR)));
    );
    filename1 << close;
  -- averages of the betti tables:     
    filename2 = concatenate(basefilename,"avgBetti",fileNameExt);
    if Z>0 then ( 
    filename2 << "SUM OF BETTI TABLES (not including Betti tabls of any zero ideals generated)" << endl;
    filename2 << beta << endl;
    filename2 << "SUM OF BETTI TABLES (with zero ideals)" << endl;
    filename2 << betaWithZeroIdeals << endl;
    filename2 << "AVERAGE BETTI NUMBERS (without zero ideals)" << endl;
    filename2 << b << endl;
    filename2 << "AVERAGE BETTI NUMBERS (with zero ideals)" << endl;
    filename2 << bWithZeroIdeals << endl;
    ) else ( -- no zero ideals so don't worry about that part: 
    filename2 << "SUM OF BETTI TABLES" << endl;
    filename2 << beta << endl;
    filename2 << "AVERAGE BETTI NUMBERS" << endl;
    filename2 << b << endl;
    );
    filename2 << close;
    print "Average Betti numbers:" expression(b);
    print "Interpretation: entry (i,j) in average Betti table encodes 
      sum_{all ideals}beta_{ij} / (sample size).";
    print "Note: the average Betti table does not include any zero ideals generated.";
  -- average betti table SHAPE: 
    filename3 = basefilename|"avgBettiShape"|fileNameExt;
    if Z>0 then ( 
    filename3 << "AVERAGE BETTI TABLE SHAPE (without zero ideals): an entry of 0.2 means 20% of ideals have a non-zero Betti number there" << endl;
    filename3 << bShape << endl;
    filename3 << "AVERAGE BETTI TABLE SHAPE (with zero ideals): an entry of 0.2 means 20% of ideals have a non-zero Betti number there" << endl;
    filename3 << bShapeWithZeroIdeals << endl;
    filename3<< "Interpretation: entry (i,j) in average Betti table SHAPE encodes 
      sum_{all ideals} 1_{beta_{ij}>0} / (sample size); that is,
      the proportion of ideals with a nonzero beta_{ij}."<<endl;
    ) else ( -- no zero ideals so don't worry about that part: 
    filename3 << "AVERAGE BETTI TABLE SHAPE: an entry of 0.2 means 20% of ideals have a non-zero Betti number there" << endl;
    filename3 << bShape << endl;
    filename3<< "Interpretation: entry (i,j) in average Betti table SHAPE encodes 
      sum_{all ideals} 1_{beta_{ij}>0} / (sample size); that is,
      the proportion of ideals with a nonzero beta_{ij}."<<endl;
    );
    filename3 << close;
    print "Average Betti table shape:" expression(bShape);
    print "Interpretation: entry (i,j) in average Betti table SHAPE encodes 
      sum_{all ideals} 1_{beta_{ij}>0} / (sample size); that is,
      the proportion of ideals with a nonzero beta_{ij}.";
    print "Note: the average Betti table shape does not include any zero ideals generated.";
    return (b,bShape)
    )

-------------------------------------------------------------------------------------
--------------------------------------Regularity-------------------------------------
-------------------------------------------------------------------------------------

--computes regularity of each RMI, prints, returns and saves to file `regularity'  - with an extension encoding values of n,p,D,N. 
-- also saves a distribution and a TALLY (i.e. histogram) of all regularities computed at the end of that file! 
avgReg = method()
avgReg (List,ZZ,String,String) :=   (ideals,N,basefilename,fileNameExt) -> (
    reg = 0;
    regHistogram={};
    filename = basefilename|"regularity"|fileNameExt;
    fileHist = basefilename|"regularityHistogram"|fileNameExt;
    apply(#ideals,i->( 
        regi = regularity ideals_i;
        reg = reg + regi;
        regHistogram = append( regHistogram, regi);
        filename << regi << endl
	)
    );
    filename << close;
    fileHist << values tally regHistogram << endl; 
    fileHist << tally regHistogram;
    fileHist<<close; 
    print "Average regularity (of non-0 ideals):" expression(sub((1/#ideals)*reg, RR));
    sub((1/#ideals)*reg, RR)
    )
{*
-- HEY HOW ABOUT MORE THAN JUST AVERAGES? 
--computes regularity of each RMI, prints, returns and saves to file `regularity'  - with an extension encoding values of n,p,D,N. 
regularityHistogram = method()
regularityHistogram (List,ZZ,String,String) :=   (ideals,N,basefilename,fileNameExt) -> (
    reg = {};
    filename = concatenate(basefilename,"regularity",fileNameExt);
    apply(#ideals,i->( 
        regi = regularity ideals_i;
        reg = append( reg, regi);
        filename << regi << endl
	)
    );
    filename << close;
    print "Average regularity (of non-0 ideals):" expression(sub((1/#ideals)*sum(reg), RR));
    print tally reg;
    )
*}
-------------------------------------------------------------------------------------
----------------------------------Krull Dimension------------------------------------
-------------------------------------------------------------------------------------

--computes Krull dim of each RMI, saves to file `dimension' - with an extension encoding values of n,p,D,N. 
--prints and returns the avg. Krull dim (real number) 
--also saves the histogram of dimensions
avgDim =   (ideals,N,Z,basefilename,fileNameExt) -> (
    dims = (numgens ring ideals_0)*Z; --since zero ideals fill the space but were not included in ideals
    dimsHistogram=toList(Z:numgens ring ideals_0);
    filename = basefilename|"dimension"|fileNameExt;
    fileHist = basefilename|"dimensionHistogram"|fileNameExt;
    apply(#ideals,i->( 
        dimi = dim ideals_i;
	filename << dimi << endl;
        dims = dims + dimi;
	dimsHistogram = append(dimsHistogram, dimi)
	)
    );
    filename << close;
    fileHist << values tally dimsHistogram << endl; 
    fileHist << tally dimsHistogram;
    fileHist<<close; 
    print "Average Krull dimension:" expression(sub(1/N*dims, RR));
    sub(1/N*dims, RR)
    )

-------------------------------------------------------------------------------------
-------------------------Minimal Generators/Degree Complexity------------------------
-------------------------------------------------------------------------------------

--computes # of min gens of each RMI, saves to file `mingens' - with an extension encoding values of n,p,D,N. 
--prints and returns the avg. # of min gens
--computes deg. complexity of each RMI, saves to file `degreecomplexity' - with an extension encoding values of n,p,D,N. 
--prints and returns avg. deg. complexity
--also saves the histograms of these invariants.
avgMinGens = method()
avgMinGens (List,ZZ,String,String) :=   (B,N,basefilename,fileNameExt) -> (
    num = 0;
    numgensHist={};
    m = 0;
    complexityHist={};
    filename1 = concatenate(basefilename,"mingens",fileNameExt);
    file1hist = concatenate(basefilename,"mingensHistogram",fileNameExt);
    filename2 = concatenate(basefilename,"degreecomplexity",fileNameExt);
    file2hist = concatenate(basefilename,"degreecomplexityHistogram",fileNameExt);
    apply(#ideals,i->( 
        mingensi = gens gb ideals_i;
        numgensi = numgens source mingensi;
        mi = max({degrees(mingensi)}#0#1);
        m = m + mi#0;
        filename1 << mingensi << endl;
        filename2 << mi#0 <<endl;
        num = num + numgensi;
	numgensHist = append(numgensHist, numgensi); 
	complexityHist = append(complexityHist, mi#0) -- (??) THE DEGREE COMPLEXITY TALLY IS WRONG. All other tallys are correct. so why this one?
	)
    );
    filename1 << close;
    filename2 << close;
    file1hist << values tally numgensHist << endl; 
    file1hist << tally numgensHist;
    file1hist<<close; 
    file2hist << values tally complexityHist << endl; 
    file2hist << tally complexityHist;
    file2hist<<close; 
    print "Average # of min gens:" expression(sub((1/N)*num, RR));
    print "Average degree complexity:" expression(sub(1/N*m, RR));
    (sub((1/N)*num, RR), sub(1/N*m, RR))
    )

-------------------------------------------------------------------------------------
---------------------------------Cohen-Macaulayness----------------------------------
-------------------------------------------------------------------------------------

--checks whether each RMI is CM, prints and returns (real number) % of CM RMIs in sample
cohMac = method()
cohMac (List,ZZ,ZZ) :=  (ideals,N,Z) -> (
    cm = 0;
    for i from 0 to #ideals-1 do (
        if isCM(R/ideals_i) == true then cm = cm + 1 else cm = cm);
    print "Percent Cohen-Macaulay:" expression(sub((cm+Z)/N, RR));
    sub((cm+Z)/N, RR)
    )

-------------------------------------------------------------------------------------
---------------------------------Borel-Fixedness-------------------------------------
-------------------------------------------------------------------------------------

--checks whether each RMI is Borel-fixed, 
--prints and returns % of Borel-fixed RMIs in sample (real number) 
borelFixed = method()
borelFixed (List,ZZ,ZZ) :=  (ideals,N,Z) -> (
    bor = 0;
    for i from 0 to #ideals-1 do ( 
        if isBorel(monomialIdeal(ideals_i)) == true then bor = bor + 1 else bor = bor);            
    print "Percent Borel-fixed:" expression(sub((bor+Z)/N, RR));
    sub((bor+Z)/N, RR) -- this is the returned value.
    )

-------------------------------------------------------------------------------------
-------------------------------------Degree------------------------------------------
-------------------------------------------------------------------------------------

--computes degree of R/I for each RMI, saves degrees to file “degree” - with an extension encoding values of n,p,D,N. 
--prints and returns avg. degree (real number)
avgDeg = method()
avgDeg (List,ZZ,String,String) :=  (ideals,N,basefilename,fileNameExt) -> (
    deg = 0;
    degHist={};
    filename = concatenate(basefilename,"degree",fileNameExt);
    fileHist = concatenate(basefilename,"degreeHistogram",fileNameExt);
    apply(#ideals,i-> ( 
        degi = degree ideals_i;
        filename << degi << endl;
        deg = deg + degi;
	degHist = append(degHist, degi)
	)
    );
    filename << close;
    fileHist << values tally degHist << endl; 
    fileHist << tally degHist;
    fileHist<<close; 
    print "Average degree:" expression(sub(1/N*deg, RR));
    sub(1/N*deg, RR)
    )

-------------------------------------------------------------------------------------
------------------------------Projective Dimension-----------------------------------
-------------------------------------------------------------------------------------

--computes proj. dim. of each RMI, saves to file “projdims” - with an extension encoding values of n,p,D,N. 
--prints and returns avg. proj dim (real number) and their histogram
avgPdim = method()
avgPdim (List,ZZ,String,String) :=  (ideals,N,basefilename,fileNameExt) -> (
    pd = 0;
    pdHist={};
    filename = concatenate(basefilename,"projdims",fileNameExt);
    fileHist = concatenate(basefilename,"projdimsHistogram",fileNameExt);
--    for i from 0 to #B-1 do ( -- wrong index set for i. #B\neq #ideals. 
    apply(#ideals,i-> 
	(
        pdimi = pdim(R^1/ideals_i);
        filename << pdimi << endl;
        pd = pd + pdimi;
	pdHist = append(pdHist, pdimi)
	)
    );           
    filename << close;
    fileHist << values tally pdHist << endl; 
    fileHist << tally pdHist;
    fileHist<<close; 
    print "Average projective dimension:" expression(sub(1/N*pd, RR));
    sub(1/N*pd, RR)
    )

-------------------------------------------------------------------------------------
----------------------------------- Summary -----------------------------------------
-------------------------------------------------------------------------------------

--write a summary of each sample to file
-- SHORT SUMMARY
summaryForPaper = (N,reg,pd,cm,basefilename, fileNameExt) -> (
    filename = concatenate(basefilename,"summary", fileNameExt);
	filename << "Average regularity (of non-0 ideals):" << endl;
	filename << reg << endl;
	filename << "Average projective dimension:" << endl;
	filename << pd << endl;
	filename << "Percent Cohen-Macaulay:" << endl;
	filename << cm << endl;
	filename << close
	)

-- LONG SUMMARY
--summary = (ideals, fileNameExt) -> (
    -- this was bad because we are not passing it any values, only ideals, and those are not at all used.
    -- so we cannot use this method in post-simulation mode...
-- so I edited the input:
--    summary = (N,averageRegularity,averagePdimension,averageKrullDimension,averageDegree,
--	averageNumGens,averageDegreeComplexity,percentCM,percentBorel,bAvg,bAvgShape,basefilename, fileNameExt) -> (
summary = (N,reg,pd,dims,deg,num,m,cm,bor,b,bShape,basefilename, fileNameExt) -> (
    filename = concatenate(basefilename,"summary", fileNameExt);
	filename << "Average regularity (of non-0 ideals):" << endl;
	filename << reg << endl;
	filename << "Average projective dimension:" << endl;
	filename << pd << endl;
	filename << "Average Krull dimension:" << endl;
	filename << dims << endl;
	filename << "Average degree:" << endl;
	filename << deg << endl;
	filename << "Average # of min gens:" << endl;
	filename << num << endl;
	filename << "Average degree complexity:" << endl;
	filename << m << endl;
	filename << "Percent Cohen-Macaulay:" << endl;
	filename << cm << endl;
	filename << "Percent Borel-fixed:" << endl;
	filename << bor << endl;
	filename << "Average Betti numbers:" << endl;
	filename << b << endl;
	filename << "Average Betti table shape:" << endl;
	filename << bShape << endl;
	filename << close
	)