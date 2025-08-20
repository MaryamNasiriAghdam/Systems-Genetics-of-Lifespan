/* ============================================================
   Title: Table S1. Lifespan data and analyses.
     (C) Analyses of variance of lifespan for the 37 parental DGRP lines. 
         All flies were reared at 25°C. Data are from Huang et al. 2020.
   Manuscript: Systems Genetics of Lifespan and Senescence in 
               Drosophila melanogaster (Nasiri Aghdam et al., 2025)
   Author: Maria E. Adonay (@amalgamaria)
   Date: 2025-08-20
   Description: Type III Analysis of Variance (ANOVA) using PROC MIXED - Reduced Model, By Sex
   SAS Version: SAS® Studio (Enterprise) 5.2 as part of SAS® Viya® 3.5
                (Copyright 2023, SAS Institute Inc., Cary, NC, USA),
                accessed through Clemson University’s licensed site.
   ============================================================ */

/* ----------------------------- */
/* Import data */

filename reffile filesrvc folderpath='/path/to/input'
  filename='filename.csv';

proc import datafile=reffile dbms=csv out=work.import;
  guessingrows=max;
  getnames=yes;
run;

proc contents data=work.import;
run;

/* ----------------------------- */
/* Sort data */

proc sort data=work.import out=work.import_sort;
  by sex;
run;

/* ----------------------------- */
/* Run model */
/* Save tables */

ods noproctitle;
ods graphics / imagemap=off;

ods output Type3=work.type3_table;
ods output CovParms=work.covparam_table;

proc mixed data=work.import_sort method=type3 alpha=0.05;
  class line sex rep;
  model lifespan= /;
  random line rep(line) /;
  by sex;
run;

/* ----------------------------- */
/* Add scientific notation to p-values */
/* Adjust precision of estimates */

data work.type3_table_scientific;
  set work.type3_table;
  ProbF_scientific = put(ProbF, BEST32. -L);
run;

data work.covparam_table_precise;
  set work.covparam_table;
  Value_precise = put(Estimate, BEST32. -L);
run;
