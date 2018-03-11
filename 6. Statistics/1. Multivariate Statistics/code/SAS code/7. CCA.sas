data acam(type=corr);
infile cards missover;
_type_='corr';
input _type_ $ _Name_ $ Y1 Y2 X1 X2;
cards;
N . 25 25 25 25
Corr Y1 1
Corr Y2 0.5 1
Corr X1 0.1 0.5 1
Corr X2 -0.5 0.1 0.5 1
;
run;


proc print data=acam;
run;

/*
 * To perform CCA is enough to just appy this line.
 * It will perform the CCA Y1 + Y2 = X1 + X2
 */
proc cancorr data=acam outstat=outacam;
var Y1 Y2;
with X1 X2;
run;
