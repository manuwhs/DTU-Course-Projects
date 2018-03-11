/*
 * In this document we implement tests for the difference in 2-sample ANOVA multivariate.
 * That is:
 *  - We have 1 Factor with 2 Classes-
 *  - We have Several correlated variables [X1, X2, X3]
 * - We want to see if the Multivatiate means of the classes is different !
 */


data heatclimate;
input sex $ height evap temp;
cards;
m 177 18.1 33.9
m 189 18.8 33.2
m 181 20.4 33.9
m 184 19.5 33.8
m 183 30.5 33.3
m 178 22.2 33.6
m 162 19.4 34.2
m 176 26.7 33.2
m 190 16.6 33.2
m 180 45.4 33.5
m 179 24.0 33.9
m 175 34.6 33.8
m 183 21.3 33.5
m 177 33.3 33.9
m 185 22.9 33.8
m 176 18.6 33.5
f 160 14.6 32.9
f 171 27.0 33.7
f 168 27.6 32.3
f 171 20.2 33.1
f 169 30.8 33.4
f 169 17.4 33.5
f 167 21.1 33.0
f 170 19.3 34.1
f 162 21.5 33.8
f 160 15.2 33.0
f 168 15.4 33.7
f 157 25.2 33.9
f 161 13.9 34.8
f 164 20.2 31.9
f 161 25.3 34.0
f 180 12.6 33.5
;
proc discrim pool=yes anova distance pcov wcov;
class sex;
var height evap temp;
run;
proc discrim pool=yes distance pcov;
class sex;
var evap temp;
run;