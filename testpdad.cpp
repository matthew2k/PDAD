/*
PROGRAM: testpdad.cpp, DOB=6/17/2013                                 LAST=1/20/2014

PURPOSE: Computing p-value by fitting a weighted sum distribution co max likelihood
*/

#include "C:\dqian\cpp\src\f1_head.h"
#include "C:\dqian\cpp\src\f_pdad.h"
#include "testpdad.h"

int main()
{
   int iPST; iPST=getTmInt();

   // (1) FOR TEST RUNS:
   pw_try(5,9,24,100,10,10.,50., iPST,"test0.out"); // time=00:00

   //pbbm(9,24, "1,2,3,4,5,11","100,200,500",3,10,10000, iPST,"test1.out"); // time=00:07

   //tsbm(9,24, "1,11","100,200",10000, iPST,"test2.out"); // time=00:28

   //nbm(2,3,1,1, 9,24,1000, iPST,"test3.out",1); // time=


   // (2) FOR MANUSCRIPT: Computing p-value by fitting weighted sum distributions
   //pbbm(9,24, "1,2,3,4,5,11","100,200,500,1000,2000",6,200,10000000, iPST,"r1_pbbm_130912.out"); // time=35:56 (was 36:31)

   //tsbm(9,24, "1,2,3,4,5,11","100,200,500",1000000, iPST,"r2_tsbm_130920.out"); // time=50:28:30

   //nbm_mrun(9,24,1000,1, iPST,"r3_nbm_1309xx.out"); // time=xx:xx:xx (was 45:11:58_hooke)

   // huge computation for sample size estimation, expected 100+ hours:
   //nbm_mrun(9,24,1000,5, iPST,"fit_nbm_100a.out"); // fig: time=117:03:07=4.9d

   return 0;
}
