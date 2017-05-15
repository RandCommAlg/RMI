-------------------------------------------------------------------------------------
-- DO NOT LOAD THIS ENTIRE FILE. IT MAY TAKE FOREVER. THIS ACTUALLY RUNS A BUNCH OF SIMULATIONS!!!
-------------------------------------------------------------------------------------




-------------------------------------------------------------------------------------
-- Random Monomial Ideals (M2 Side)
-- Sept. 2016
-- Dane Wilburne


--- SP EDITS:18.9.2016.
-- run the following in M2 directly:
restart
-------------------------------------------------------------------------------------
------------------------------ Load Packages
--  [ moved from other sript file b/c can't reload the same package twice get error ]
-------------------------------------------------------------------------------------

--loadPackage("depth") --just use depth.m2
-- ^^ doesn't work for Sonja, so instead grabbed the isCM function and saved to local file:  
--
load("depth.m2")
loadPackage("BoijSoederberg")



{*
n=4  -- set number of vars
D=10  -- set max degree
N=1000 -- set desired sample size
p=0.1 -- set desired probability for the model.


-- -- path = prepend("/Users/petrovic/Dropbox/Random Monomial Ideals simulations - Sonja local",path)
basefilename="" -- change this if you want to save the output files to another directory;
-- for example i am using this:
basefilename = "/Users/petrovic/Dropbox/Random Monomial Ideals simulations - Sonja local/SimulationsOutput/"
load("RMIscript.m2")
*}


{*  
-- simulations needed for the paper: 
for n=2,3,4,5:
    for D=10,20,30,40:
         for 20 p values evenly spaced from D^(-n) to D^(-1):
              compute C-M, regularity, proj. dim
*}	       	      
	      
-- ITERATIVE SIMULATION RUNS AND STORING IN SEPARATE FOLDERS :
--
-- This is optional but I want *all* of my outputs in the folder below (else main folder gets flooded w/ files):
basefilename = "/Users/petrovic/Dropbox/Random Monomial Ideals simulations - Sonja local/SimulationsOutput/temp/";
basefilename = "/Users/sxp61/Dropbox/Random Monomial Ideals simulations - Sonja local/SimulationsOutputHomeDesktop/temp/";
-- in fact the folder above gets *deleted* upon end of iteration (!) and outputs are moved to a different place "storing parent folder"
-- This next one is used below to store each iteration into a sub-folder: 
iterativeStoringParentFolder = "/Users/petrovic/Dropbox/Random Monomial Ideals simulations - Sonja local/SimulationsOutput/";
iterativeStoringParentFolder = "/Users/sxp61/Dropbox/Random Monomial Ideals simulations - Sonja local/SimulationsOutputHomeDesktop/";
-- function for removing of directories function from M2 doc: 
rm = d -> if isDirectory d then removeDirectory d else removeFile d;




-- fix D, n, N and vary p: 
apply(20,i-> (
    -- make a tempororary place to store all files
    makeDirectory basefilename; 
    -- run i-th simulation: 
    n=2;
    N=1000;
    D=20;
    --for 20 p values evenly spaced from D^(-n) to D^(-1)  -- [20 values including endpts]
    p= sub(D^(-n),RR)+i*sub((D^(-1)-D^(-n))/19,RR);
    load("RMIscript.m2"); 
    -- move results to a place you want to store: 
    copyDirectory(basefilename, iterativeStoringParentFolder|"p"|toString(p)|"D"|toString(D)|"/");
    --    removeDirectory basefilename:
    scan(reverse findFiles basefilename, rm)
    )
);


-- an attempt to make 4 rows of the following (1st) table:
{*  
-- simulations needed for the paper: 
for n=2,3,4,5:
    for D=10,20,30,40:
         for 20 p values evenly spaced from D^(-n) to D^(-1):
              compute C-M, regularity, proj. dim
*}	       	      
apply(4,i-> (
	-- first select the D because that's the one that's varying across the row
     D=10+10*i;
     -- then, fix D, n, N and vary p: 
apply(20,i-> (
    -- make a tempororary place to store all files
    makeDirectory basefilename; 
    -- run i-th simulation: 
    n=2; -- set number of variables 
    N=500; -- set sample size 
    --for 20 different values of p, evenly spaced from D^(-n) to D^(-1)  -- [20 values including endpts]
    p= sub(D^(-n),RR)+i*sub((D^(-1)-D^(-n))/19,RR);
    load("RMIscript.m2"); 
    -- move results to a place you want to store: 
    copyDirectory(basefilename, iterativeStoringParentFolder|"p"|toString(p)|"D"|toString(D)|"/");
    --    removeDirectory basefilename:
    scan(reverse findFiles basefilename, rm)
    )
);
)
)



-- AN ATTEMPT TO GET 1ST 2 ROWS OF THE 2ND TABLE, WITH n INCREASING (!) -- shortened simulations: 
{*  
-- simulations needed for the paper: 
for D=2,3,4,5 (at least the 1st two!): 
    for n=10,20,30,40:
         for 20 p values evenly spaced from D^(-n) to D^(-1):
              compute C-M, regularity, proj. dim
*}	       	      

-- first row; D=2, n = 10,20,30,40:       	      
apply(4,i-> (
     n=10+10*i;
     -- then, fix D, n, N and vary p: 
apply(20,i-> (
    -- make a tempororary place to store all files
    makeDirectory basefilename; 
    -- run i-th simulation: 
    D=3; -- fixed across the row 
    N=500; -- set sample size (lower than the usual 1000) 
    --for 20 different values of p, evenly spaced from D^(-n) to D^(-1)  -- [20 values including endpts]
    p= sub(D^(-n),RR)+i*sub((D^(-1)-D^(-n))/19,RR);
    load("RMIscriptForPaper.m2"); 
    -- move results to a place you want to store: 
    copyDirectory(basefilename, iterativeStoringParentFolder|"p"|toString(p)|"D"|toString(D)|"/");
    --    removeDirectory basefilename:
    scan(reverse findFiles basefilename, rm)
    )
);
)
)

