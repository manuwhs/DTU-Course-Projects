libname train '/folders/myfolders/data';
data train;
set train.train;
run;

/* 
 * We print the first 25 random observations.
 * Only the variables MJ1 and from MR1, MR2, ..., MR6
 */
Title 'First 25 out of 100 randomly chosen observations';
proc print data=train (obs=25);
var MJ1 MR1-MR6;
run;

/* 
 * We train the full GLM model
 */


Title 'Full model';
proc reg data=train;
model MJ1 = MR1-MR6;
run;


/* 
 * We train the full GLM model perform selection of variables by RSquare.
 * We can indidicate different elimination methods:
 * 		- selection = stepwise
 * 		- selection = Rsquare
 *  	- selection = Forward
 *  	- selection = Backward
 */

Title 'Model selected by Backward Eliminiation';
proc reg data=train;
model MJ1 = MR1-MR6/selection=Backward;
run;


Title 'Full model with intercept';
proc reg data=train;
model MJ1 = MR1-MR6/covb;
run;


/* 
 * We use intercept and the second variable
 */

Title 'Model Band2';
proc reg data=train;
model MJ1 = MR2/covb;
run;