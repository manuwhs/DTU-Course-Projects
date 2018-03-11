* Exam 2008 problem 1;
* The data are part of a larger study by Dr. Rick Linthurst while he was 
doing his PhD at Noth Carolina State University. 
He considered the aerial biomass of a certain marsh grass Spartina Alterniflora in the
Cape Fear Estuary of North Carolina.
There are corresponding observations of aerial biomass, and soil characteristics: salinity, acidity, potassium,
sodium, zinc. Furthermore, there are recordings of locality and type of vegetation.;

Data spartina;
input Obs Location $ Type $ Biomass Salinity pH K Na Zn;
datalines;
1 OI DVEG 676 33 5.00 1441.67 35184.5 16.4524
2 OI DVEG 516 35 4.75 1299.19 28170.4 13.9852
3 OI DVEG 1052 32 4.20 1154.27 26455.0 15.3276
4 OI DVEG 868 30 4.40 1045.15 25072.9 17.3128
5 OI DVEG 1008 33 5.55 521.62 31664.2 22.3312
6 OI SHRT 436 33 5.05 1273.02 25491.7 12.2778
7 OI SHRT 544 36 4.25 1346.35 20877.3 17.8225
8 OI SHRT 680 30 4.45 1253.88 25621.3 14.3516
9 OI SHRT 640 38 4.75 1242.65 27587.3 13.6826
10 OI SHRT 492 30 4.60 1281.95 26511.7 11.7566
11 OI TALL 984 30 4.10 553.69 7886.5 9.8820
12 OI TALL 1400 37 3.45 494.74 14596.0 16.6752
13 OI TALL 1276 33 3.45 525.97 9826.8 12.3730
14 OI TALL 1736 36 4.10 571.14 11978.4 9.4058
15 OI TALL 1004 30 3.50 408.64 10368.6 14.9302
16 SI DVEG 396 30 3.25 646.65 17307.4 31.2865
17 SI DVEG 352 27 3.35 514.03 12822.0 30.1652
18 SI DVEG 328 29 3.20 350.73 8582.6 28.5901
19 SI DVEG 392 34 3.35 496.29 12369.5 19.8795
20 SI DVEG 236 36 3.30 580.92 14731.9 18.5056
21 SI SHRT 392 30 3.25 535.82 15060.6 22.1344
22 SI SHRT 268 28 3.25 490.34 11056.3 28.6101
23 SI SHRT 252 31 3.20 552.39 8118.9 23.1908
24 SI SHRT 236 31 3.20 661.32 13009.5 24.6917
25 SI SHRT 340 35 3.35 672.15 15003.7 22.6758
26 SI TALL 2436 29 7.10 528.65 10225.0 0.3729
27 SI TALL 2216 35 7.35 563.13 8024.2 0.2703
28 SI TALL 2096 35 7.45 497.96 10393.0 0.3205
29 SI TALL 1660 30 7.45 458.38 8711.6 0.2648
30 SI TALL 2272 30 7.40 498.25 10239.6 0.2105
31 SM DVEG 824 26 4.85 936.26 20436.0 18.9875
32 SM DVEG 1196 29 4.60 894.79 12519.9 20.9687
33 SM DVEG 1960 25 5.20 941.36 18979.0 23.9841
34 SM DVEG 2080 26 4.75 1038.79 22986.1 19.9727
35 SM DVEG 1764 26 5.20 898.05 11704.5 21.3864
36 SM SHRT 412 25 4.55 989.87 17721.0 23.7063
37 SM SHRT 416 26 3.95 951.28 16485.2 30.5589
38 SM SHRT 504 26 3.70 939.83 17101.3 26.8415
39 SM SHRT 492 27 3.75 925.42 17849.0 27.7292
40 SM SHRT 636 27 4.15 954.11 16949.6 21.5699
41 SM TALL 1756 24 5.60 720.72 11344.6 19.6531
42 SM TALL 1232 27 5.35 782.09 14752.4 20.3295
43 SM TALL 1400 26 5.50 773.30 13649.8 19.5880
44 SM TALL 1620 28 5.50 829.26 14533.0 20.1328
45 SM TALL 1560 28 5.40 856.96 16892.2 19.2420
;

/*
 * Q1.2: WE compute the conditional correlation between the variables and ph
 */


proc corr data=spartina;
var Biomass K;
partial ph;
run;
/*
 * Q1.3: WE compute the correlation between the variables and use them in the formula to compute
 * the partial correlation needed.
 * We first so it with the correlations and then with the covariance and we see that we get the same result
 */
proc corr data=spartina;
var Salinity pH Zn;
run;

* Q 1.4;
/* We are going to first compute the correlation matrix with proc corr.
Then we will read the results, extract the matrices and compute the multiple correlation value

*/
proc corr data=spartina outp=test;
run;

proc iml;  * General procedure;

/* Read the correlation matrix from the total data outputed by proc corr*/
use test;
   read all var _NUM_ into testData; 
close test;

print testData;
/* Obtain the submatrices to compute their determinant */
Sigma_i = testData[5:10,2:7];
sigma_ii = testData[5,2];
Sigma_22 = testData[6:10,3:7];
print Sigma_i sigma_ii Sigma_22;
multcorr = 1 - det(Sigma_i)/(sigma_ii * det(Sigma_22));
print multcorr;
run;



proc iml;
Sigma_i = {
	435703.71	-253.31	637.24	-40200.66	-1235995.43	-3412.59,
	-253.31	13.84	-0.24	-22.84	4154.12	-12.96,
	637.24	-0.24	1.55	7.14	-323.73	-7.46,
	-40200.66	-22.84	7.14	88567.15	1622390.93	181.38,
	-1235995.43	4154.12	-323.73	1622390.93	47367751.47	6669.93,
	-3412.59	-12.96	-7.46	181.38	6669.93	68.56
	};
sigma_ii = 435703.71;
Sigma_22 = {
	13.84	-0.24	-22.84	4154.12	-12.96,
	-0.24	1.55	7.14	-323.73	-7.46,
	-22.84	7.14	88567.15	1622390.93	181.38,
	4154.12	-323.73	1622390.93	47367751.47	6669.93,
	-12.96	-7.46	181.38	6669.93	68.56
};
print Sigma_i sigma_ii Sigma_22;
multcorr = 1 - det(Sigma_i)/(sigma_ii * det(Sigma_22));
print multcorr;
run;

/*
proc factor data=spartina nfactors=3 rotate=varimax;
var Salinity pH K Na Zn;
run;
*/
