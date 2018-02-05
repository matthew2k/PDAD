/*
PROGRAM: f1_head.h, 
PURPOSE: Include all std::xxx, function headers & header files for C++ and g++44 programs.
*/


#ifndef F1_HEAD_H
#define F1_HEAD_H

////////////////////////////////////////////////////////////////////////////////
// #1: include C++ headers and global parameters
////////////////////////////////////////////////////////////////////////////////

#define _CRT_SECURE_NO_DEPRECATE // disable 'safe' stringFunInMsvs2005 (strncpy_s strcpy_s)

#include <iostream>  // cout,cin,cerr,endl,flush ios,fixed,left,right,showpoint,ostream,istream
#include <iomanip>   // setw, setprecision, setiosflags, resetiosflags, setfill
#include <fstream>   // ifstream, ofstream
#include <cstring>   // strcpy, strncpy
#include <cstdlib>   // srand, rand, exit
#include <new>       // new
#include <cmath>     // pow
#include <ctime>     // time
#include <ctype.h>   // isdigit, isalpha
#include <stdio.h>   // sprintf_s
#include <string>    // string
#include <sstream>   // string manipulation: ostringstream

using namespace std;

#define PI 3.14159265358979
#define INVSQRT2PI 0.398942280401433
#define LS 50000


////////////////////////////////////////////////////////////////////////////////
// #2: include downloaded headers/packages
////////////////////////////////////////////////////////////////////////////////

// RNALIB: downloaded/modified from http://www.netlib.org/random/ranlib.c.tar.gz on 1/3/2012:
#include "ranlib.h"

// BOOST: downloaded from http://www.boost.org/users/history/version_1_52_0.html on 1/4/2013
// (1) For Windows, use boost_1_52_0.7z, VERY ACCURATE nonCtrChiSqDist, rv functions, etc.
// To add "boost" directory into a VisualStudioSolution: Project >> Properties >> C/C++
// >> Additional Include Directories >> enter "C:/.../boost_1_52_0" >> OK
// (2) For Linux,

// BOOST: include headers (works for Windows C++, not work for Linux g++44 yet: 5/30/2013)
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/weibull.hpp>

using namespace boost::math;


////////////////////////////////////////////////////////////////////////////////
// #3: include function prototypes BEFORE all internal header files
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// [3] RANDOM NUMBER FUNCTIONS of f3_rand.h
////////////////////////////////////////////////////////////

double Unfm();                            // uniform random variable (rv) U(0,1)
double UnfmAny(double, double);           // uniform rv U(a,b)

double Norm();                            // standard normal rv N(0,1)
double Norm(double mu, double sigma);     // normal rv N(mu, sigma)
double NormR(double, double, double, double);// restricted N(me,sd,mi,ma)
double NormDF(double z);                  // standard normal density function
double NormCDF(double z);                 // standard normal CDF (p from z)
double NormNlp2S(double z);               // -log10(p) of 2-sided normal pval (nlp from z)
double NormNlpR(double z);                // -log10(p) of rt tailed normal pval (nlp from z)
double NormInv(double p);                 // inverse standard normal CDF (z from p)

double rchisq(double df, double nc);       // chisq_1of4.RANLIB: rv of nc chisq dist w/ df+nc
double dchisq(double x, double, double);    // chisq_2of4.BOOST: density (PDF f from xQtl)
double pchisq(double x, double, double);    // chisq_3of4.BOOST: probability (CDF pLt from xQtl)
double qchisq(double p, double, double);    // chisq_4of4.BOOST: quantile (invCDF xQtl from pLt)

double rgamma(double a, double b);         // gamma_1of4.RANLIB: rv w/ a=shape, b=scale=1/rate
double dgamma(double x, double a, double b);// gamma_2of4.BOOST: density=PDF
double pgamma(double x, double a, double b);// gamma_3of4.BOOST: probability=CDF
double qgamma(double p, double a, double b);// gamma_4of4.BOOST: quantile=invCDF

double rlnorm(double mu, double sigma);    // lnorm_1of4: rv w/ mu+sigma
double dlnorm(double x, double, double);    // lnorm_2of4: density=PDF
double plnorm(double x, double, double);    // lnorm_3of4: probability=CDF
double qlnorm(double p, double, double);    // lnorm_4of4: quantile=invCDF

double rweibull(double a, double b);       // weibull_1of4: rv w/ a=shape, b=scale
double dweibull(double x, double, double);  // weibull_2of4.BOOST: density=PDF
double pweibull(double x, double, double);  // weibull_3of4.BOOST: probability=CDF
double qweibull(double p, double, double);  // weibull_4of4: quantile=invCDF

double rgumbel(double a, double b);        // gumbel_1of4: rv w/ a=location, b=scale
double dgumbel(double x, double, double);   // gumbel_2of4: density=PDF
double pgumbel(double x, double, double);   // gumbel_3of4: probability=CDF
double qgumbel(double p, double, double);   // gumbel_4of4: quantile=invCDF

double pgpareto(double x, double a, double b);// gpareto: probability=CDF


double Expn(double);                      // exponential rv co theta
double Gamma_int(int a, double b);         // (naive) gamma rv co int shape a & dbl scale b
double Gamma(double a, double b);         // gamma rv co dbl shape a>0 & dbl scale b>0
double Gamma_1mass(double p0, double x0, double a, double b); // gamma rv w/ 1 mass
double Beta(double, double);              // beta rv w/ a & b

int RdInt(int, int);                      // random integer b/ iMin and iMax
int RdIntEx1(int, int, int);                // random integer except 1 value
int RdIntEx2(int, int, int, int);            // random integer except 2 values
int RdIntPr(double *, int);               // random int 0:n-1 co pr[] w/ sum=1
int RdIntPrChk(double *, int);            // random int 0:n-1, check if sumPr is 1
int RdIntPrAny(double *, int);            // random int 0:n-1 co any pr[]
int RdIntPrSum(double *, double, int);      // random int 0:n-1 co pr[] & sumPr (faster)
void RdIntPrWR(double *, int, int *, int);   // get m int 0:n-1 w/ prob+replacement
void get2RdIntWR(int, int, int &, int &);    // get 2 random int w/ replacement
void get2RdIntNR(int, int, int &, int &);    // get 2 random int wo replacement
void get4RdIntNR(int, int, int &, int &, int &, int &);
void getMltRdIntWR(int, int, int *, int);    // get m random int w/ replacement
void bootstrap(double *v, int s, double *vb); // get bootstrapped sample vb[] from v[]
void getMltRdIntNR(int, int, int *, int);    // get m random int wo replacement
void asgnArray012(int *, int);             // assign int array w/ 0,1,...,s-1
void asgnArray012Perm(int *, int);         // assign int array w/ permutation of 0,1,...,s-1

int Pois(double);                         // poisson rv w/ mean lambda
int PoisT(double, int, int);              // truncated poisson rv
int PoisR(double, int, int);              // restricted poisson rv
int Binm(double, int, int);               // binomial w/ probability p1
int BinmAnyPr(double, double, int, int);     // binomial w/ any prob desity p1, p2
int Geom(double);                         // geometric random function
int GeomR(double, int, int);              // restricted geometric random function

double PoisDF(int, double);               // poisson density function
double PoisCDF(int, double);              // poisson cumulative distribution function
double PoisCDFT(int, double, double);       // poisson CDF trucated at pmax
void MltnmRdDF(int N, double *df);        // multinomial random density function


										  ////////////////////////////////////////////////////////////
										  // [4] UTILITY FUNCTIONS
										  ////////////////////////////////////////////////////////////

										  // [4.1] PRINT MESSAGE, STRING/MATRIX, TIME, SOFTWARE HEADER, ETC
void prtMsg_exit(const char *msg);              // prt msg & exit program
void prtMsg1i_exit(const char *msg, int i1);     // prt msg+1int & exit program

void errata1(int errCode);                      // err #1 code & exit program
void errata2(int errCode, char *vNm);           // err #2 code+varnm & exit program
void prt_err_msg(char *c);                      // err #3 msg & exit program
void errata4(char *c1, int i1, char *c2);        // err #4 msg+int+msg & exit program
void errata5(char *, int, char *, int, char *);     // err #5 msg+(int+msg)*2 & exit program
void erratai(char *vNm, int v);                 // err #i invalid vnm+vInt & exit program
void erratad(char *vNm, double v);              // err #d invalid vnm+vDbl & exit program

void errOpenIF_exit(char *fn);                  // prtErrMsgForCannotOpen
void errOpenIF_exit(const char *fn);
void errOpenOF_exit(char *fn);
void errOpenOF_exit(const char *fn);
void getCsz(char *fn, int &cSz, int &ncMi, int &ncMa); // get colSz+nDlmCol ofIptFile
void getRsz(char *fn, int &rSz, int &rEmpty);      // get rowSz+rowEpt ofIptFile
char getDlm(char *fn);                          // get dlm=(' ', ',', etc) ofIptFile
void prtFileInfo(ostream &out, char *fn);        // prt fileInfo=rSz+cSz+fmt ofIptFile

void setFixedShowPoint(ostream &, int);          // set fixed show point & precision
void prtDbl(ostream &, double);                  // prt double value in non-fixed format
void prtDblFW(ostream &, double, int);            // prt dbl in nonFixFmt+fixWidth
void prtDblFWLt(ostream &out, double d, int w);   // prt dbl in nonFixFmt+fixWidth+LtAlign

void prtCh(ostream &, char, int);                 // prt chars w/ specified length
void prtChNR(ostream &, char, int);

void prtHmsSecNR(ostream &out, int sec);         // hh:mm:ss from int sec w/o return
void prtHmsSec(ostream &out, int sec);           // hh:mm:ss from int sec w/ return
int getMilliSecUnitCoSystem();                  // get unit of millisecond co system
void prtHmmMsecNR(ostream &out, int ms);         // hh:mm:ss.sss from ms w/o return
void prtHmmMsec(ostream &out, int ms);           // hh:mm:ss.sss from ms w/ return

void prtHmsNR(ostream &out, int iBeg);           // hh:mm:ss from int iBeg w/o return
void prtHms(ostream &out, int iBeg);             // hh:mm:ss from int iBeg w/ return

void prtMdyHmsNR(ostream &out, time_t tIn);      // prtDtTm as mm/dd/yyyy hh:mm:ss woRtn
void prtMdyHms(ostream &out, time_t tIn);        // prtDtTm as mm/dd/yyyy hh:mm:ss wRtn
void prtMdyHmsNR(ostream &out);                 // prtCurDtTm as mm/dd/yyyy hh:mm:ss woRtn
void prtMdyHms(ostream &out);                   // prtCurDtTm as mm/dd/yyyy hh:mm:ss wRtn

void prtHmsDttmNR(ostream &out, int iBeg);       // prt hh:mm:ssFromBeg & curDtTm woRtn
void prtHmsDttm(ostream &out, int iBeg);         // prt hh:mm:ssFromBeg & curDtTm wRtn


												 // [4.2] SEARCH, COUNT, COPY/MODIFY/COMPARE & CONVERSION OF CHARACTER ARRAYS
int index(string s1, const char *s2);         // search s1 for 1st occurrence of s2
int index(char *s1, char *s2);                // search s1 for 1st occurrence of s2
int index(const char *s1, char *s2);
int index(char *s1, const char *s2);
int index(const char *s1, const char *s2);
int index(char *s1, const char *s2, int n);    // search s1 for nth occurrence of s2
int indexUC(char *s1, char *s2);             // search upcase & compressed s1 for 1st occu of s2
int indexWd(char *s1, const char *s2);        // search s1 for fstWd as whole s2
int indexAnyWdS2(char *s1, const char *s2);   // search s1 for fstAnyWd in s2
int indexc(char *s1, char *s2);               // search s1 for 1st occurrence of any char in s2
int indexc(char *s1, const char *s2);
int indexc(char *s1, const char *s2, int n);   // search s1 for nth occurrence of any char in s2
int indexci(char *s1, int i1);                // search s1 for 1st occurrence of i1
int indexi(int i1, int i2);                   // search integer i1 for 1st occurrence of i2

int strCount(char *s1, const char *s2);       // count occurrences of string s2 in string s1
int strCountAny(char *s1, char *s2);         // count occurrences of any s2 char in string s1
int dlmCount(char *s1, const char *s2);       // count s1 for any non-adjacent delimiters in s2
int wordCount(string sIn);                   // count no. of words in sIn
int wordCount(char *sIn);                    // count no. of words in sIn
int wordCount(const char *sIn);

void compressc(char *, char *, const char *);  // strcpy by removing any char in s1
void strcpy_rmbk(char *, char *);             // strcpy by removing beg+end+multi blks
void strcpy_rmbk_op(char *, char *);          // strcpy by remove b+e+multi+flanking_ops blks
void strcpy_up(char *sOut, char *sIn);       // strcpy with letters replaced by upper ones
void strcpy_up(char *sOut, const char *sIn);
void strcpy_upcp(char *sOut, char *sIn);     // strcpy w/ letters replaced by upper+compressd
void strcpy_lw(char *sOut, char *sIn);       // strcpy w/ letters replaced by lower ones
void strcpy_an(char *, char *, int);           // strcpy after column n (n=0 for strcpy)
void strcpy_an(char *, const char *, int);
void strcpy_bn(char *, char *, int);           // strcpy before col n
void strcpy_bn(char *, const char *, int);
void strcpy_bman(char *, char *, int, int);     // strcpy before col m & after col n (m <= n)
void strcpy_ambn(char *, char *, int, int);     // strcpy after col m & before col n (m+2 <= n)
void strcpy_fmtn(char *, char *, int, int);     // strcpy from col m to col n (m <= n)
void strcpy_adlm(char *, char *, const char *);// strcpy after any delimiter
void strcpy_eadlm(char *, char *, char *);     // strcpy eql+aft any delimiter
void strcpy_bdlm(char *, char *, const char *);// strcpy before any const delimiter
void strcpy_wbfd(char *, char *, const char *);// strcpy within & before any dlm
void strcpy_wbup(char *, char *, const char *);// strcpy within & before upcase dlm of s1
void strcpy_adbd(char *, char *, const char *, const char *);// strcpy aft/bfr dlm s1/s2
void strcpy_astr(char *, char *, char *);      // strcpy after string s1
void strcpy_astr(char *, char *, const char *);
void strcpy_asup(char *, char *, const char *);// strcpy after upcaseConstString s1
void strcpy_upas(char *, char *, const char *);// strcpy upcase after upcase string s1
void strcpy_bstr(char *, char *, const char *);// strcpy before string s1
void strcpy_bsup(char *, char *, const char *);// strcpy before upcsConstString s1
void strcpy_a1b2(char *, char *, const char *, const char *);//strcpy aft/bfr s1/s2
void strcpy_abup(char *, char *, const char *, const char *);//strcpy aft/bfr upcsConst s1/s2
void strcpy_baup(char *, char *, const char *, const char *);//strcpy bfr/aft upcsConst s1/s2
void strcpy_asbd(char *, char *, const char *, const char *);//strcpy aftStr s1 & bfrAnyDlm in s2

void strcpy_addC(char *, char *, char *);      // strcpy and add characters at the end
void strcpy_addC(char *, char *, const char *);
void strcpy_addC(char *, const char *, char *);
void strcpy_apd_sLk(string &, string, const char *);// string_app & add sLk in btw

void strcpy_newE(char *, char *, const char *); // strcpy and renew ext at the end
void strcpy_del4CAddC(char *, char *, char *);//strcpy, del last 4 chars, add new chars
void strcpy_newiExt(char *, char *, int);      // strcpy and renew 4-char ext from integer
void strcpy_addiExt(char *, char *, int);      // strcpy and add 4-char ext from interger
void iToC(int iIn, char *sOut);              // convert integer to character string
void strcpy_addI(char *, char *, int);         // strcpy and add integer at the end
void strcpy_addI_pad0(char *, char *, int, int);// strcpy and add int w/ padded zeros
void strcpy_addIC(char *, char *, int, char *); // strcpy and add int & chars at the end
void strcpy_addCI(char *, char *, char *, int); // strcpy and add chars & int at the end
void strcpy_addCI_pad0(char *, char *, char *, int, int);//strcpy, add chars & int pad 0
void strcpy_addI_keepE(char *, char *, int);   // strcpy,add int bfr ".",keep ext
void strcpy_addC_keepE(char *, char *, char *);// strcpy, add chars bfr ".", keep ext
void strcpy_addC_keepE(char *, const char *, char *);
void strcpy_addC_keepE(char *, char *, const char *);

void strcpy_addC_newE(char *, char *, char *, char *);//strcpy, add char bfr ".", renew ext
void strcpy_addC_newE(char *, const char *, char *, char *);
void strupd_BCNA(char *, char *, char *, int &);// strupd by backup,copy, (index of) new/app

int strcmp_extNCS(char *fn, const char *ext); // strcmp ext NonCaseSensitive

											  // STRING CONVERSION TO/FROM OTHER DATA TYPES:
string dToS(double x);                                      // cvt dbl to string
string dToS_sci(double x, int dec);                          // cvt dbl to string w/ scifmt
string dToS_sciFix(double x, double xmiS, int decS, int decF); // cvt dbl to str w/ sciOrFixFmt


															   // [4.3] INPUT VARIABLE, ARRAY, LINE
void ipt_avn(char *, const char *, int &);    // input after var name for 1 integer
void ipt_avn(char *, const char *, double &); // input after var name for 1 double
void ipt_avn_s1(char *, const char *, char *);// input after var name for 1st string
void ipt_avn_str(char *, const char *, char *);// input after var name for whole string
void ipt_avn_sa(char *, const char *, char **, int &);//iptAftVarNm for StrAyIn()By,
void ipt_att(char *sIn, char *vOut);         // input after "title" for 1 string

void ipt_bn(char *sIn, int n, int &vOut);      // input integer before column n
void ipt_an(char *sIn, int n, int &vOut);      // input integer after column n
void ipt_nth_int(char *s, int n, int &vOut); // input nth integer
void ipt_nth_int(const char *s, int n, int &vOut);
void ipt_nth_dbl(char *s, int n, double &vOut);// input nth double
void ipt_nth_dbl(const char *s, int n, double &vOut);

void ipt_nth_word(char *sIn, int n, char *vOut);
void ipt_nth_word(const char *sIn, int n, char *vOut);
void ipt_nth_word(string sIn, int n, char *vOut);// input nth word
void ipt_nth_wdup(char *sIn, int n, char *vOut);
void ipt_nth_wdup(string sIn, int n, char *vOut);// input nth word & upcase

int ipt_sz(char *sIn);                       // input anyIntDbl arraySz (w/o array)
void ipt_ia(char *sIn, int *v, int &vSz);      // input int array+sz
void ipt_ia(const char *sIn, int *v, int &vSz);
void ipt_da(char *sIn, double *v, int &vSz);   // input dbl array
void ipt_da(const char *sIn, double *v, int &vSz);
void ipt_ia_avn(char *, char *, int *v, int &vSz);//input int array after vNm
void ipt_da_avn(char *, char *, double *v, int &vSz);//input dbl array after vNm
void ipt_ia_ml(istream &, char *, int *, int);  //input int array in 1/multi lines
void ipt_da_ml(istream &, char *, double *, int);//input dbl array in 1/multi lines

void ipt_akw_szPli_chkSz(char **Pgm, int Npl, const char *kw, int szMa, int &vSz, int &Pli);//ipt akw intDblArSz
void ipt_akw_i1Pli(char **, int, const char *, int &v1, int &Pli);//ipt akw for 1int+pgmLnIdx
void ipt_akw_i1Pli_errNF(char **, int, const char *, int &v1, int &Pli);//ipt 1int+Pli & prtErrForNotFd
void ipt_akw_i1_Nu(char **, int, const char *, int &, int);        //ipt akw 1int w/ nuVal
void ipt_akw_iaSzPli(char **, int, const char *, int *, int &, int &);//ipt akw intArray+Sz+Pli
void ipt_akw_iaSz(char **, int, const char *, int *, int &);       //ipt akw intArray+Sz
void ipt_akw_ia_chkSz(char **, int, const char *, int *, int);     //ipt akw intArray w/ chkSz
void ipt_akw_ia_chkSzRg(char **, int, const char *, int *, int, int, int);//ipt akw intArray w/ chkSz+chkRg
void ipt_akw_i2(char **, int, const char *, int &, int &);         //ipt akw 2int
void ipt_akw_i4(char **, int, const char *, int &, int &, int &, int &);//ipt akw 4int

void ipt_akw_d1Pli(char **, int, const char *, double &, int &);   //ipt akw 1dbl+Pli
void ipt_akw_d1_Nu(char **, int, const char *, double &, double);  //ipt akw 1dbl w/ nuVal
void ipt_akw_d1_chkRg(char **, int, const char *, double &, double, double);//ipt akw 1dbl w/ chkRg
void ipt_akwPli_daSzPli(char **, int, const char *, int, double *, int &, int &);//ipt akw+pliP dblAr+Sz+Pli
void ipt_akw_daSzPli(char **, int, const char *, double *, int &, int &);//ipt akw dblArray+Sz+Pli
void ipt_akw_daSz(char **, int, const char *, double *, int &);    //ipt akw dblArray+Sz
void ipt_akw_daSz_chkSz(char **, int, const char *, double *, int &, int, int); //ipt akw da+sz w/ chkSz
void ipt_akw_d2(char **, int, const char *, double &, double &);   //ipt akw 2dbl
void ipt_akw_d4(char **, int, const char *, double &, double &, double &, double &); //ipt akw 4dbl

void ipt_akw_i1d1_Nu(char **Pgm, int Npl, const char *kw, int &i1, double &d1, int i1Nu, double d1Nu); //ipt 1i1d_nu
void ipt_akw_i1d2_Nu(char **, int, const char *, int &, double &, double &, int, double, double);//ipt 1i2d_nu
void ipt_akw_d1i1_Nu(char **Pgm, int Npl, const char *kw, double &d1, int &i1, double d1Nu, int i1Nu);//ipt_1d1i_nu

void prt_akw_err(char **Pgm, const char *kw, int Pli);//prt akw for kwNotFd or errKw
void prt_akw_err_msg(char **Pgm, const char *kw, int Pli, const char *msg);//prt akw for errKw+msg
void prt_akw_errSz(char **Pgm, const char *kw, int Pli, int szData, int szTrue);//prt akw for kwNotFd/errKw/szWrg
int chk_kw(char **Pgm, int Npl, const char *kw); //chk kw from pgmLn
int chk_akw_sa(char **Pgm, int Npl, const char *kw, const char *saI);//chk akw forEachStrInStrArr saI
int chk_akw_saAny(char **Pgm, int Npl, const char *kw, const char *saI);//chk akw forAnyStrInStrArr saI

void ipt_akw_s1(char **Pgm, int Npl, const char *kw, char *s1, int &Pli); //ipt akw for 1str
void ipt_akw_s2(char **Pgm, int Npl, const char *kw, char *s1, char *s2, int &Pli);
void ipt_akw_s3(char **Pgm, int Npl, const char *kw, char *s1, char *s2, char *s3, int &Pli);
void ipt_akw_s4(char **Pgm, int Npl, const char *kw, char *s1, char *s2, char *s3, char *s4, int &Pli);

int getDsz(char *fn, int vhead);                          // get dataSz=nObs
void ipt_var_i1(char *fn, int vhead, int vn, int *v);      // iptVarArFmDfFor1IntVar
void ipt_var_d1(char *fn, int vhead, int vn, double *v);   // iptVarArFmDfFor1DblVar


														   // [4.4] OPEN+MANAGE FILE
void openOF(ofstream &out, const char *fn);
void openOF(ofstream &out, char *fn);
void openOFNewExt(ofstream &out, const char *fn, char *sExt);
void openOFAddExt(ofstream &out, const char *fn, char *sExt);
void openOFNewExtInt(ofstream &out, const char *fn, int iExt);
void openOFAddExtInt(ofstream &out, const char *fn, int iExt);
void openOF_addSuf_keepE(ofstream &out, const char *fn, char *sSuf);
void apndOF(ofstream &out, const char *fn);
void apndOF(ofstream &out, char *fn);
void apndOFNewExt(ofstream &out, const char *fn, char *sExt);
void openApndOF(ofstream &out, const char *fn, int o1a2);
void openApndOF(ofstream &out, char *fn, int o1a2);
void openApndOF_addC_keepE(ofstream &out, const char *fn, char *sNew, int o1a2);

void apnd2Files(char *f1, char *f2, char *fout, int n1l, int n1u, int n2l, int n2u);
void apnd2Files(char *f1, char *f2, char *fout);


// [4.5] CALCULATION FUNCTIONS
// [4.5a] SIMPLE FUNCTIONS AND MATRIX COMPUTATION
int sign(double x);                    // 1=(x>0), 0=(x=0), -1=(x<0)
double sqrt(int i);
double logit(double p);                // logit(p)=(inverse of logistic)
double logistic(double x);             // logistic(x)=(inverse of logit)

void sumArrayConst(double *a, int s1, double csum); // update array for const sum
void MatSub(double **, int, int, int, double **);//complementary sub-matrix of a_ij
double MatDet(double **, int);         // get determinant value of a matrix
void MatInv(double **, double **, int);  // get the inverse matrix of a matrix
int sharedItv(int *h1, int *h2, int M, int m0); // sharedItvSurMkBtw2Hts

												// STATISTICAL TESTS: FISHER'S EXACT TEST, 2 GROUP COMPARISONS, ETC
double lfact(int n);                                  // get log-factorial of n
double pFisherEtTab(int a, int b, int c, int d);         // get probFisherExactTest of 2by2 table
double pFisherEt2Side(int a0, int b0, int c0, int d0);   // get 2-sided pvalue of FisherExactTest

														 // RANK, PERCENTAGE AND PERCENTILE
int RankLg(double *xs, int n, double x);                // get rank from largest
int RankSm(double *xs, int n, double x);                // get rank from smallest

double Percentage(double *xs, int n, double x);         // get pctg of x in xs[]
void pctgVec(double *xs, int n, double *pctg);         // get pctg[] of xs[]
double pctgCorr(double *xs, int n);                    // corr(xs[],pctg[])

double Percentile(double *xs, int n, double p);        // get pctl of p in xs[]


													   // SIMPLE ALGORITHMS
double newton_raphson(double f, double fd, double x0, double eps, int itm); // Newton-Raphson mtd


																			////////////////////////////////////////////////////////////////////////////////
																			// #4: include internal header files
																			////////////////////////////////////////////////////////////////////////////////

#include "f2_temp.h"
#include "f3_rand.h"
#include "f4_util.h"
#include "f5_pdad.h"


#endif
