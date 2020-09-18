-- ******************************************************************
-- GENERATE GB training DATA 
-- run this section of code once to get BOTH
-- the random binomial ideals saved as a text file,
-- and the size of each minimal GB saved in a different text file. 
-- ******************************************************************


generateGBdata = method(TypicalValue => List, Options=>{Homogeneous=>false,OneAtATime=>true,InfoLine=>true})--IdealsOnly=>false,
generateGBdata(ZZ,ZZ,ZZ,ZZ) := List => o -> (numVars,maxDegree,binomialsInEachSample,sampleSize) -> ( 
    S = ZZ/32003[x_0..x_(numVars-1)];
    filename = "RandomBinomialDataSet."|toString numVars|"vars.deg"|toString maxDegree|".sampleSize"|toString sampleSize|"."|toString binomialsInEachSample|"binomialsEach."|toString currentTime()|".txt";
    --for foolproof automation: --  filename=temporaryFileName()|".txt"; 
    --just add a "version number" at the end of the text file if one of the same name already exists. We do have a time stamp, so that should take care of most conflicts, but just in case, double-fool-proof?: [adding random # at the end])
    while  fileExists(filename)  do   filename = "RandomBinomialDataSet."|toString numVars|"vars.deg"|toString maxDegree|".sampleSize"|toString sampleSize|"."|toString currentTime()|".v"|toString random(100)|".txt";
    
    -- store the preset parameters for this data set on the first line of the data file: 
    parameters = "numVars = "|toString numVars|", maxDegree = "|toString maxDegree|", binomialsInEachSample = "|toString binomialsInEachSample|", sampleSize = "|toString sampleSize;
    if o.InfoLine then (
	    f := openOut filename;
	    f << parameters << endl;
	    close f;
	    );
    
    if not o.OneAtATime then (
    	expos = []; 
    	ideals ={};
    	scan(sampleSize,i-> (
		-- generate a new random binomial set: 
		bins = new Array from randomBinomials(S,maxDegree,binomialsInEachSample,Homogeneous=>o.Homogeneous);
		assert(#bins == binomialsInEachSample); -- just to make sure I got the correct sample size; if wrong it'll print error on dialog so easy to spot!
		-- save exponents to the list for learning: 
		expos = expos|apply(bins,b->new Array from apply(exponents b,monexpo->new Array from monexpo));
		-- save ideals so we can compute GB down below: 
		ideals = append(ideals,ideal toList bins);
		)
    	    );
	f = openOutAppend filename;
    	f << replace(" ","",expos) << endl; 
	close f;
    	gbSizes = [];
    	-- TO DO: 
    	-- for generating larger samples, include a timer and throw out an ideal we can't compute in a few minutes. 
    	-- how to best handle this? try "try" command :D see Slack!  
    	scan(ideals, I-> gbSizes = gbSizes | [# flatten entries gens gb I]);
    	--gbSizes
    	assert(#gbSizes ==sampleSize); -- just to make sure we have same # of gbs as we do samples of binomial sets. 
    	concatenate(filename,".gbSizes.txt") << gbSizes << endl << close;

    ) else ( 
    	-- generate a random binomial set, write it to file f, then compute gb write size gb, and only THEN continue to next random set.
	scan(sampleSize,i-> (
		-- generate a new random binomial set: 
		bins = new Array from randomBinomials(S,maxDegree,binomialsInEachSample,Homogeneous=>o.Homogeneous);
		assert(#bins == binomialsInEachSample); -- just to make sure I got the correct sample size; if wrong it'll print error on dialog so easy to spot!
		-- save exponents to the open file f: 
		expos =  apply(bins,b->new Array from apply(exponents b,monexpo->new Array from monexpo));	
    		f := openOutAppend filename;
		f << replace(" ","",toString(expos));
		close f;
		gbf := openOutAppend concatenate(filename,".gbSizes.txt");
		gbf <<  [# flatten entries gens gb ideal toList bins];
		close gbf;
		)
	    );
	-- Eventually we will time out the operations and if gb doesn't complete we can save a 0 or a -1 in the gb sizes file. 
	);
    print("Saved data to file starting with "|filename);
    )

   
   end;
    -- STOP READING THE FILE.
   --  WHAT'S BELOW NEEDS TO BE EDITED 
    -- it's the stupid option of getting a bunch of binomial sets  but NOT  computing their gb. deal with later. 

--getFilename = () -> (
--  filename := temporaryFileName();
--  while fileExists(filename) 
--    or fileExists(filename|"PHCinput") 
--    or fileExists(filename|"PHCoutput") do filename = temporaryFileName();
--  filename
--)



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
writeData(replace(" ","",expos), parameters, filename);

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

