/* This example:
	- creates a mini dataset and 
	- prints a subset of variables 
*/

proc univariate data=stat2.sundhed plot freq normal;
var vegt;
/* This example will output:
	- The main momentums of the variable (mean, variance, skewness, kurtosis...)
	- Basic statistical measures (mean, median, range...)
	- Statistical test to test if:
		- The mean is significantly bigger than 0
		- The sign is random. (Same number of + variables and - variables)
		- The signed rank which I dont know what it is.
	- Some statistical tests for Normality to check if the samples are normally distributed
	- Quantiles, Extreme Observations and frequency counts.
	- Q-Q plot for normality
*/

proc plot data=stat2.sundhed; 
plot vegt*alder;
/* This program makes a quick and dirt scatter plot of the variables 
*/

proc sgplot data=stat2.sundhed;
scatter x=vegt y=alder;

/* This program makes a better scatter plot of the variables 
*/

proc corr data=stat2.sundhed cov;
/*
 * Compute the correlation and Covariance of the dataset.
 * For the Correlation it computes the statistical test that they are individually 0.
 * If the p-values is too low we reject the null hypothesis, and they are statistically significant
 */


proc princomp data=stat2.sundhed cov;
/*
 * Compute the eigenvectors and eigenvalues
 */
