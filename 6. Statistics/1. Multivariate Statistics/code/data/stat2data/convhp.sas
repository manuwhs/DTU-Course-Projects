filename andban 'andban.fil';
filename andban2 'andban2.fil';
filename astma 'astma.fil';
filename iris 'iris.fil';
filename iriscali 'iriscali.fil';
filename iristest 'iristest.fil';
filename opgze5 'opg_e5.fil';
filename sundhed 'sundhed.fil';
filename xmpl3z5 'xmpl3_5.fil';
libname stat2 '.';

data stat2.andban;
infile  andban;
        input
          x1
          x2
          x3
          y
          x1x2;

data stat2.andban2;
infile  andban2;
        input
          x1
          x2
          x3
          y
          x1x1
          x1x2
          x1x3
          x2x2
          x2x3
          x3x3;

data stat2.astma;
infile  astma;
        input
          middel $
          respuge1
          respuge2
          respuge3
          respuge4;

data stat2.iris;
infile  iris;
        input
          sepallen
          sepalwid
          petallen
          petalwid
          species $;

data stat2.iriscali;
infile  iriscali;
        input
          sepallen
          sepalwid
          petallen
          petalwid
          species $;

data stat2.iristest;
infile  iristest;
        input
          sepallen
          sepalwid
          petallen
          petalwid
          species $;

data stat2.opg_e5;
infile  opgze5;
        input
          aar
          kvartal
          stoerrel
          selskab
          lagerbev;

data stat2.sundhed;
infile  sundhed;
        input
          alder
          vegt
          ilt
          loebetid
          hvilpuls
          loebpuls
          maxpuls;

data stat2.xmpl3_5;
infile  xmpl3z5;
        input
          bakterie
          buffer
          forsoeg
          udbytte;

proc contents data=stat2._all_ position;

run;
