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
 * We train a simple model where we only use the MR1.
 *  - We will analyze:  
 *   TODO
 */


Title 'Simple Model';
proc reg
data=train
plots(label)=(CooksD RStudentByLeverage DFFITS DFBETAS);
model MJ1 = MR1/r clb cli clm influence;
run;
