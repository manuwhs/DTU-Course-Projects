/* This example:
	- Loads a dataset from the library "stat2" which is just a folder containing
	datasets and code.
	- prints a subset of variables 
	
	- Loads a correlation matrix instead of real data
*/

/* 
 * A library in SAS is just a folder with code and databases.
 * We indicate the existance of the library with the "libname" statement and we tell where it is.
 * We can access the files using name_library.name_file.
 * We would usually load the libraries in the AutoExec file so that we do not need to add the line in the programs
 */
libname stat2 '/folders/myfolders/data/stat2data'; * where stat2's data resides;
proc print data=stat2.sundhed; * ells SAS to print the data  tells SAS which dataset to print ;
var alder vegt; *tells SAS which variables in the dataset to print
run; *tells SAS that no more inputs are coming and to execute the PROC statementproc print tells SAS to print the data;

/*
 * SAS allows to load sufficient statistics of the data instead of the data itself.
 * If we assume Gaussianity, then the sufficient data is: Number of samples, Mean, STDs and Corr matrix
 */
data cement(type=corr);
infile datalines missover;
input _type_ $ _name_ $ BLAINE XSO3 LSR STRGTH3 STRGTH7 STRGTH28 ;
datalines;
N nobs 198 198 198 198 198 198
MEAN . 3095.71717 1.94126 40.29268 263.89394 382.28788 515.27778
STD . 234.22485 0.25448 8.55410 36.88332 36.70451 32.67124
CORR BLAINE 1
CORR XSO3 0.49077 1
CORR LSR -0.19604 -0.52069 1
CORR STRGTH3 0.80042 0.47398 -0.17769 1
CORR STRGTH7 0.73643 0.37344 -0.07175 0.89557 1
CORR STRGTH28 0.51413 0.27931 -0.08501 0.60826 0.74561 1
;
Title 'Cement Data';
proc print data=cement;
run;

/*
 * Then we can apply any other function that only needs the sufficient statistics.
 */
proc princomp data=cement cov; 