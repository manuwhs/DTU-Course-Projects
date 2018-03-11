/* This example:
	- We see an example of a 2-sided ANOVA, that is that we just have 2 Factors !!
	- The data is given in the format:
	    - (Combinations of Factors) + List of Y with that combination.
	- TODO: I have to learn more about the interaction plots !
*/


data cement;
input mixer crusher @;
do i=1 to 4;
input strength @;
output;
end;
datalines;
1 1 28 52 -24 80
1 2 -66 -60 2 120
1 3 -84 18 32 -40
2 1 -58 28 58 -10
2 2 34 -12 -4 120
2 3 -82 -20 -40 -52
3 1 36 116 68 50
3 2 72 -24 62 56
3 3 -54 -7 -32 60
;
proc print data=cement;
run;
ods graphics on;
proc glm plot=meanplot(cl);
class mixer crusher;
model strength=mixer crusher
mixer*crusher;
lsmeans mixer crusher;
run;
ods graphics off;
