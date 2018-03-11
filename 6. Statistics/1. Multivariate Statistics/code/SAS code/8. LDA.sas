
libname lab '/folders/myfolders/data';
data lab2;
set lab.labm1m2;
if _N_ < 36 then pr=0;
else pr=1;
run;

/*
 * Plot the loaded data
 */
/*
proc print data=lab2;
run;
*/

/*
 * Perform LDA:
 * 	- The class is indicated by the variable PR
 * 	- The variables of the distribution are indicated by var 
 *  - pool = yes means that we assume the same covariance matrix for the classes.
 *  - Priors indicate the needed priors.
*/

proc discrim data=lab2 pool=yes;
class pr;
var LM2 aM2 bM2;
priors '0'=0.666666  '1'=0.33333333;
run;
/*
 * Perform LDA:
 * 	- Performs Crossvalidation 
*/

proc discrim data=lab2 pool=yes crossvalidate;
class pr;
var LM2 aM2 bM2;
run;