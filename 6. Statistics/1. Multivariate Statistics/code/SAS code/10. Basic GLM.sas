/* This example:
	- We will create first 2 different modelling of the variables for the database.
		- 1st Model: 
		  we give a Factor and 2 variables, the vector of data is [F,X1,X2]
		  We will just perform a GLM in which we predict the yield based on the separate parameters
		  and in the interaction X1(F).
		  Guess there is no interaction between X1 and X2 since they 
		  uncorrelated by definition, but could there be a linear relationship between it an X2 that is 
		  not captured and we need to also specify ? 
		    · I tried it and yes, the error decreases althought the parameters are not statistically significant.
	    
	    - 2nd Model: 
	   	  Now we create separate variables for X1 depending on the values of F. That is:
	   	  F is split into different vectors (F1, F2,...) according to its cardinality.
	   	  And then the variable X1 is also split, but not X2. 
	   	  This way we will be able to analyze the individual constribution of the separate classes of F ?
	   	  In the specification of the model now we do not need to indicate any interaction,
	   	  the independent fake variables F1, F2, F1X1, F2X1 acount for the interaction FX1
	   	  
	  THESE 2 MODELS are equivalent in that they yield the same total result, but in the second one
	  we can break down the contributions even further knowing how much they contribute.
	  	- Their parmeters are the same.
	  	- But we can compute the type I and type III error in the second one with more detail ?
	  	
	  	
	  Overall statements:
	     - The more parameters we have -> The least SSres. Names Error in the Tables. 
	     	Probably the mean error could decrease if a parameter is non informative, and then we would be dividing by more.
	     	
*/

/* 1st Model */

Data penicillin;
input yield sugar $ conc
sqconc;
cards;
0.606 lac -3 9
0.660 lac -1 1
0.984 lac 1 1
0.908 lac 3 9
0.761 can -3 9
0.933 can -1 1
1.072 can 1 1
0.979 can 3 9
;
run;
proc glm;
class sugar;
model yield = sugar
conc(sugar) sqconc/noint
solution;
run;

/* 2nd Model */

data penicillin2;
input yield intlac conclac
intcan conccan sqconc;
cards;
0.606 1 -3 0 0 9
0.660 1 -1 0 0 1
0.984 1 1 0 0 1
0.908 1 3 0 0 9
0.761 0 0 1 -3 9
0.933 0 0 1 -1 1
1.072 0 0 1 1 1
0.979 0 0 1 3 9
;
run;
proc glm;
model yield=intlac conclac
intcan conccan
sqconc/noint solution;
run;

