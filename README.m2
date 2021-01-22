-- ******************************************************************
-- README  
-- ******************************************************************


-- ******************************************************************
-- ******************************************************************
-- option 1  for  generating  random toric ideals: kernel of  random monomial maps.
-- ******************************************************************
-- ******************************************************************

-- check the file  "randomBinomialIdeals.m2" starting at line 63!  


-- ******************************************************************
-- ******************************************************************
-- option 2  for  generating  random binomial ideals: differences of random monomials. 
-- ******************************************************************
-- ******************************************************************

-- Random binomial data is written to a FILE. The filename for this data set is of
-- the format "RandomBinomialDataSet.3vars.deg4.sampleSize10". 
-- The format of the file is: 
-- 1st line contiains parameters for the data generated,
-- 2nd line is a list of lists of binomial exponents, printed as string from array,


-- ******************************************************************
-- PREREQUISITES 
-- ******************************************************************
load"randomBinomialIdeals.m2"

load"generateBinomialDataSets.m2" 


-- ******************************************************************
-- PRESETS -- this is setup by hand for each data-generating run -- 
-- set up what kind of data you want 
-- ******************************************************************
numVars = 3
maxDegree = 2 
homogeneous = true --  do you prefer to have binomials all of maxDegree (true), or up to maxDegree (false)?
binomialsInEachSample = 5 -- how many binomials in each sample 
sampleSize = 3  -- how many samples of the above many binomials each.


-- ******************************************************************
-- DATA GENERATION CALLS 
-- ******************************************************************

-- probably the way we will be running this is : 
generateGBdata(numVars,maxDegree,binomialsInEachSample,sampleSize,Homogeneous=>homogeneous,OneAtATime=>true)

-- change homogeneous flat to check it's working: 
generateGBdata(numVars,maxDegree,binomialsInEachSample,sampleSize) -- ,Homogeneous=>false, OneAtATime=>true  are default anyway.

-- maybe you hate the "info line" at the start of the data file? you can remove it:
generateGBdata(numVars,maxDegree,binomialsInEachSample,sampleSize, InfoLine=>false)

-- not writing random binomials & gbs to files one at a time (not recommended except for small examples!) : 
generateGBdata(numVars,maxDegree,binomialsInEachSample,sampleSize,Homogeneous=>homogeneous, OneAtATime=>false) 
    