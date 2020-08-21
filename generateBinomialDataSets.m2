-- ******************************************************************
-- GENERATE GB training DATA 
-- run this section of code once to get BOTH
-- the random binomial ideals saved as a text file,
-- and the size of each minimal GB saved in a different text file. 
-- ******************************************************************
filename = "RandomBinomialDataSet."|toString numVars|"vars.deg"|toString maxDegree|".sampleSize"|toString sampleSize|"."|toString currentTime()|".txt";
  --for foolproof automation: --  filename=temporaryFileName()|".txt"; 
  --just add a "version number" at the end of the text file if one of the same name already exists. We do have a time stamp, so that should take care of most conflicts, but just in case, double-fool-proof?:)
  while  fileExists(filename)  do   filename = filename|"v"|toString random(10)|".txt";
-- store the preset parameters for this data set on the first line of the data file: 
parameters = "numVars = "|toString numVars|", maxDegree = "|toString maxDegree|", binomialsInEachSample = "|toString binomialsInEachSample|", sampleSize = "|toString sampleSize;
filename << parameters << endl << close;
   
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

-- writing data into a file one at a time
	 f := openOut filename;
	 line := "what to put into the file, assuming we want one entry per line";
	 f << line << endl;	 
	 close f;



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

