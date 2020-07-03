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
numVars = 4  
maxDegree = 5 
homogeneous = true --  or false
binomialsInEachSample = 10 -- how many binomials in each sample 
sampleSize = 100  -- how many samples of the above many binomials each.

-- ******************************************************************
-- PREREQUISITES 
-- ******************************************************************
load"randomBinomialIdeals.m2"


-- ******************************************************************
-- GENERATE DATA -- run this sectino of code once to get 
-- ******************************************************************
filename = concatenate("RandomBinomialDataSet.",toString numVars,"vars.deg",toString maxDegree,".sampleSize",toString sampleSize);
parameters = concatenate("numVars = ",toString numVars,", maxDegree = ",toString maxDegree,", binomialsInEachSample = ",toString binomialsInEachSample,", sampleSize = ",toString sampleSize);
expos = {}; 
S = QQ[x_0..x_(numVars-1)];
scan(sampleSize,i-> (
	bins = randomBinomials(S,maxDegree,binomialsInEachSample,Homogeneous=>homogeneous);
	assert(#bins == binomialsInEachSample); -- just to make sure I got the correct sample size; if wrong it'll print error on dialog so easy to spot!
	expos = append(expos, apply(bins,b->exponents b));
	)
    )
writeData(expos, parameters, filename);
 




-- ******************************************************************


-- TO DO 
-- [by Sonja]
-- I still want to write some code to use our Sample and Statistics part of the random package
-- because that way we can easliy get histograms of info on the data we are generating
-- rather than writing mroe code to get such data summaries. But it's not necessary
-- to get started, so I'll do it when I'm able. 

