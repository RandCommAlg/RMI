-- ******************************************************************
-- README  
-- ******************************************************************

-- Random binomial data is written to a FILE. The filename for this data set is of
-- the format "RandomBinomialDataSet.3vars.deg4.sampleSize10". 
-- The format of the file is: 
-- 1st line contiains parameters for the data generated,
-- 2nd line is a list of lists of binomial exponents, printed as string from array,


-- ******************************************************************
-- PRESETS -- this is setup by hand for each data-generating run -- 
-- set up what kind of data you want 
-- ******************************************************************
numVars = 3
maxDegree = 2 
homogeneous = true --  or false
binomialsInEachSample = 5 -- how many binomials in each sample 
sampleSize = 3  -- how many samples of the above many binomials each.

-- ******************************************************************
-- PREREQUISITES 
-- ******************************************************************
load"randomBinomialIdeals.m2"

load"generateBinomialDataSets.m2" 


-- ******************************************************************
-- here is an example: 
bins = randomBinomials(QQ[a,b,c],2, 3,Homogeneous=>false)
-- for some reason this is broken at the moment though it's meant to showcase what the hell is going on. 
-- but the method (not this simple example) seems to be working so WHATEVER. 
expos = []; 
expos = expos| apply(bins,b->new Array from apply(exponents b,monexpo->new Array from monexpo));
--
-- END simple example. 
-- ******************************************************************
