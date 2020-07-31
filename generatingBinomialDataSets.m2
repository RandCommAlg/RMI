-- ******************************************************************
-- README  
-- ******************************************************************
-- Random binomial data is written to a FILE. The filename for this data set is of
-- the format "RandomBinomialDataSet.3vars.deg4.sampleSize10". 
-- The format of the file is: 
-- 1st line contiains parameters for the data generated,
-- 2nd line is a list of lists of binomial exponents, as in this:
bins = randomBinomials(QQ[a,b,c],2, 3,Homogeneous=>false)
expos = apply(bins,b->exponents b) 
-- if you do not like all the {}s used, we can easily remove them one level at a time:
flatten expos
flatten flatten expos  -- etc. 
-- We can also easily put new lines after each set of binomials, or whatever. 
--
-- END OF README. 
-- ******************************************************************



-- ******************************************************************
-- PRESETS -- this is setup by hand for each data-generating run -- 
-- set up what kind of data you want 
-- ******************************************************************
numVars = 3
maxDegree = 10 
homogeneous = true --  or false
binomialsInEachSample = 5 -- how many binomials in each sample 
sampleSize = 15  -- how many samples of the above many binomials each.

-- ******************************************************************
-- PREREQUISITES 
-- ******************************************************************
load"randomBinomialIdeals.m2"


-- ******************************************************************
-- GENERATE GB training DATA 
-- run this section of code once to get BOTH
-- the random binomial ideals saved as a text file,
-- and the size of each minimal GB saved in a different text file. 
-- ******************************************************************
filename = concatenate("RandomBinomialDataSet.",toString numVars,"vars.deg",toString maxDegree,".sampleSize",toString sampleSize,".",toString currentTime(),".txt");
parameters = concatenate("numVars = ",toString numVars,", maxDegree = ",toString maxDegree,", binomialsInEachSample = ",toString binomialsInEachSample,", sampleSize = ",toString sampleSize);
expos = []; 
ideals ={};
S = QQ[x_0..x_(numVars-1)];
scan(sampleSize,i-> (
	bins = new Array from randomBinomials(S,maxDegree,binomialsInEachSample,Homogeneous=>homogeneous);
	assert(#bins == binomialsInEachSample); -- just to make sure I got the correct sample size; if wrong it'll print error on dialog so easy to spot!
	expos = expos| apply(bins,b->new Array from apply(exponents b,monexpo->new Array from monexpo));
	ideals = append(ideals,ideal toList bins);
	)
    )
print concatenate("writing to file ",filename);
writeData(expos, parameters, filename);
gbSizes = [];
-- TO DO: 
-- for generating larger samples, include a timer and throw out an ideal we can't compute in a few minutes. 
-- how to best handle this? 
scan(ideals, I-> gbSizes = gbSizes | [# flatten entries gens gb I]);
--gbSizes
assert(#gbSizes ==sampleSize) -- just to make sure we have same # of gbs as we do samples of binomial sets. 
concatenate(filename,".gbSizes.txt") << gbSizes << endl << close;




-- ******************************************************************
-- GENERATE binomials DATA 
-- run this section of code once 
-- to get JUST the random binomial ideals saved as a text file. 
-- See below for more! 
-- ******************************************************************
filename = concatenate("RandomBinomialDataSet.",toString numVars,"vars.deg",toString maxDegree,".sampleSize",toString sampleSize,".",toString currentTime());
parameters = concatenate("numVars = ",toString numVars,", maxDegree = ",toString maxDegree,", binomialsInEachSample = ",toString binomialsInEachSample,", sampleSize = ",toString sampleSize);
expos = []; 
S = QQ[x_0..x_(numVars-1)];
scan(sampleSize,i-> (
	bins = new Array from randomBinomials(S,maxDegree,binomialsInEachSample,Homogeneous=>homogeneous);
	assert(#bins == binomialsInEachSample); -- just to make sure I got the correct sample size; if wrong it'll print error on dialog so easy to spot!
	expos = expos| apply(bins,b->new Array from apply(exponents b,monexpo->new Array from monexpo));
	)
    )
print concatenate("writing to file ",filename);
writeData(expos, parameters, filename);

-- notes for future:
	--expos = append(expos, apply(bins,b->exponents b)); -- this makes everything with {}. We prefer [].
	-- At the moment, since we know we only have binomials: 
	expos = expos| apply(bins,b->[new Array from (exponents b)_0,new Array from (exponents b)_1]);
	-- but if we genereate a 0 then this is a problem; it assumes there is no 0 in the list! 
	-- if we want this to work for general polynomials, replace above line by the below one:
	--expos = expos| apply(bins,b->new Array from apply(exponents b,monexpo->new Array from monexpo));
-- ******************************************************************
 




-- ******************************************************************


-- TO DO 
-- [by Sonja]
-- I still want to write some code to use our Sample and Statistics part of the random package
-- because that way we can easliy get histograms of info on the data we are generating
-- rather than writing mroe code to get such data summaries. But it's not necessary
-- to get started, so I'll do it when I'm able. 

