# CompatibilityTests
For comparing two datasets to determine if they are consistent with each other. 

MakeCombineDatacard.C will make the datacards to imput into CMS HiggsCombineTool to perform the likelihood tests.  

UsePPAsComparison.C will run a standalone test for both chi2 and KS.  It generates the toy distributions and can output the p-values. 

CompareCompatibility.C will take the output from HiggsCombine and the standalone test and compare the results (as well as making the plots).  

HiggsCombine instructions are here (which uses the 2dNLL approach): 
https://twiki.cern.ch/twiki/bin/view/CMS/HiggsCombineTwoDatasetCompatibility
