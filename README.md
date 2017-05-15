# RandomMonomialIdeals in Macaulay2
RMI summer 2017 IIT working group

RandomMonomialIdeals.m2 is our main package file.

RMIs.m2 contains other functions (raw!) that need to be ported over to the package. The two script files (RMIscript.m2 and script.m2) are used to run simulations on random monomial ideals and are NOT a part of the package, but parts of them may be come part of the main documentation node (the top-level one), and can also be used to see how these methods work in practice.

To try this out, run the following 2 commands (after you download RandomMonomialIdeals.m2) : 

load "RandomMonomialIdeals.m2"

randomGeneratingSets(2,3,0.2,10)      

See if you understand what the output is. See comments in the .m2 file. We will discuss this Friday 5/12. 

Then to install a package:

installPackage("RandomMonomialIdeals")

or even better (to update the documentation): installPackage("RandomMonomialIdeals", RerunExamples=>true)
