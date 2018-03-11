/*
 * In this document we implement tests for the difference in 2-sample ANOVA multivariate.
 * That is:
 *  - We have 1 Factor with 2 Classes-
 *  - We have Several correlated variables [X1, X2, X3]
 * - We want to see if the Multivatiate means of the classes is different !
 */


proc iml;
/*Reading the matrices corresponding to the estimation*/
invxx=
{1.55920 -0.16549 -0.47258 -0.05010 0.41826 -0.00235 -0.42289,
-0.16549 0.85139 -0.17981 -0.01327 0.63774 -0.01759 -0.69467,
-0.47258 -0.17981 1.77862 -0.10728 -0.29340 0.01164 -0.02184,
-0.05010 -0.01327 -0.10728 0.02253 0.12325 -0.00441 -0.17012,
0.41826 0.63774 -0.29340 0.12325 5.25546 -0.08437 -7.04885,
-0.00235 -0.01759 0.01164 -0.00441 -0.08437 0.00243 0.11182,
-0.42289 -0.69467 -0.02184 -0.17012 -7.04885 0.11182 10.11541};
theta=
{0.28400 0.42731,
0.79508 0.22230,
-0.02573 0.02607,
-0.01151 0.06290,
-0.14467 -0.16756,
0.00307 0.01103,
0.10614 0.03463};
sigma=
{0.85897 0.07870,
0.07870 0.29444};

/*Reading the matrices corresponding to testing theta(4,2)=0*/
A={0 0 0 1 0 0 0};
B={0 1};
C={0};
/*Compute U-test statistic with DF=(1, 1, 169) and F-test statistic with
DF=(1, 169), t=1*/
delta=A*theta*B`-C;
r=169*sigma;
e=B*r*B`;
h=delta`*inv(A*invxx*A`)*delta;
lambda=det(e)/det(e+h);
f=169*(1-lambda)/lambda;
print lambda f;
run;








proc iml;
/*Reading the matrices corresponding to the estimation*/
invxx=
{1.55920 -0.16549 -0.47258 -0.05010 0.41826 -0.00235 -0.42289,
-0.16549 0.85139 -0.17981 -0.01327 0.63774 -0.01759 -0.69467,
-0.47258 -0.17981 1.77862 -0.10728 -0.29340 0.01164 -0.02184,
-0.05010 -0.01327 -0.10728 0.02253 0.12325 -0.00441 -0.17012,
0.41826 0.63774 -0.29340 0.12325 5.25546 -0.08437 -7.04885,
-0.00235 -0.01759 0.01164 -0.00441 -0.08437 0.00243 0.11182,
-0.42289 -0.69467 -0.02184 -0.17012 -7.04885 0.11182 10.11541};
theta=
{0.28400 0.42731,
0.79508 0.22230,
-0.02573 0.02607,
-0.01151 0.06290,
-0.14467 -0.16756,
0.00307 0.01103,
0.10614 0.03463};
sigma=
{0.85897 0.07870,
0.07870 0.29444};

/*Reading the matrices corresponding to testing theta(4,1)=theta(4,2)=0*/
A={0 0 0 1 0 0 0};
B={1 0,
0 1};
C={0 0};
/*Compute U-test statistic with DF=(2, 1, 169) and F-test statistic with
DF=(2, 168), t=1*/
delta=A*theta*B`-C;
r=169*sigma;
e=B*r*B`;
h=delta`*inv(A*invxx*A`)*delta;
lambda=det(e)/det(e+h);
f=(168/2)*(1-lambda)/lambda;
print lambda f;
run;












proc iml;
/*Reading the matrices corresponding to the estimation*/
invxx=
{1.55920 -0.16549 -0.47258 -0.05010 0.41826 -0.00235 -0.42289,
-0.16549 0.85139 -0.17981 -0.01327 0.63774 -0.01759 -0.69467,
-0.47258 -0.17981 1.77862 -0.10728 -0.29340 0.01164 -0.02184,
-0.05010 -0.01327 -0.10728 0.02253 0.12325 -0.00441 -0.17012,
0.41826 0.63774 -0.29340 0.12325 5.25546 -0.08437 -7.04885,
-0.00235 -0.01759 0.01164 -0.00441 -0.08437 0.00243 0.11182,
-0.42289 -0.69467 -0.02184 -0.17012 -7.04885 0.11182 10.11541};
theta=
{0.28400 0.42731,
0.79508 0.22230,
-0.02573 0.02607,
-0.01151 0.06290,
-0.14467 -0.16756,
0.00307 0.01103,
0.10614 0.03463};
sigma=
{0.85897 0.07870,
0.07870 0.29444};


/*Reading the matrices corresponding to testing theta(5,1)=…=theta(7,2)=0*/
A={0 0 0 0 1 0 0,
0 0 0 0 0 1 0,
0 0 0 0 0 0 1};
B={1 0,
0 1};
C={0 0,
0 0,
0 0};
/*Compute U-test statistic with DF=(2,3,169) and F-test st. with DF=(6,336), t=2*/
delta=A*theta*B`-C;
r=169*sigma;
e=B*r*B`;
h=delta`*inv(A*invxx*A`)*delta;
lambda=det(e)/det(e+h);
f=(336/6)*(1-sqrt(lambda))/sqrt(lambda);
print lambda f;
run;








proc iml;
/*Reading the matrices corresponding to the estimation*/
invxx=
{1.55920 -0.16549 -0.47258 -0.05010 0.41826 -0.00235 -0.42289,
-0.16549 0.85139 -0.17981 -0.01327 0.63774 -0.01759 -0.69467,
-0.47258 -0.17981 1.77862 -0.10728 -0.29340 0.01164 -0.02184,
-0.05010 -0.01327 -0.10728 0.02253 0.12325 -0.00441 -0.17012,
0.41826 0.63774 -0.29340 0.12325 5.25546 -0.08437 -7.04885,
-0.00235 -0.01759 0.01164 -0.00441 -0.08437 0.00243 0.11182,
-0.42289 -0.69467 -0.02184 -0.17012 -7.04885 0.11182 10.11541};
theta=
{0.28400 0.42731,
0.79508 0.22230,
-0.02573 0.02607,
-0.01151 0.06290,
-0.14467 -0.16756,
0.00307 0.01103,
0.10614 0.03463};
sigma=
{0.85897 0.07870,
0.07870 0.29444};

/*Reading the matrices corresponding to testing all theta(i,j)=0*/
A=I(7);
B=I(2);
C=J(7,2,0);
/*Compute U-test statistic with DF=(2,7,169) and F-test st. with F=(14,336), t=2*/
delta=A*theta*B`-C;
r=169*sigma;
e=B*r*B`;
h=delta`*inv(A*invxx*A`)*delta;
lambda=det(e)/det(e+h);
f=(336/14)*(1-sqrt(lambda))/sqrt(lambda);
print lambda f;
run;
