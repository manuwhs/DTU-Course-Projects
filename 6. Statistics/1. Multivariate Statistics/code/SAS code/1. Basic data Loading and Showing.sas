/* This example:
	- creates a mini dataset and 
	- prints a subset of variables 
*/

data xample;
input x1 x2 x3; * two variables;
datalines;
45.9 98.4 52
3 56.3 41
45.3 -42  89
;
run;
proc print data=xample;
var x1 x2;
run; 