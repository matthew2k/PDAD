/*
PROGRAM: f2_temp.h, 
PURPOSE: All FUNCTION TEMPLATES
*/


#ifndef F2_TEMP_H
#define F2_TEMP_H

// [2.1] MEMORY ALLOCATION OF ARRAYS (templates)
template<class T> void Array1(T *&a, int I);                // 1-d array
template<class T> void Array1(T *&a, int I, T aIni);        // 1-d array w/ iniVal
void Array1(char *&a, int I,const char *aIni);              // 1-d charAr w/ iniVal
template<class T> void A1epd(T *&a, int Ipre,int Inew);     // 1-d array expand size
template<class T> void Drray1(T *&a, int I);                // 1-d deallo array

template<class T> void Array2(T **&a, int J, int I);        // 2-d array
template<class T> void Array2(T **&a, int J, int I,T aIni); // 2-d array w/ iniVal
void Array2(char **&a, int J, int I,const char *aIni);      // 2-d charAr w/ iniVal
void Array2(string **&a, int J, int I,string aIni);    // 2-d strAr w/ iniVal
template<class T> void A2ep2(T **&a, int J, int Ipre,int Inew); // 2-d ar expand size2
template<class T> void Drray2(T **&a, int J, int I);        // 2-d deallo array

template<class T> void Array3(T ***&a, int K, int J, int I); // 3-d array
template<class T> void Array3(T ***&a, int K, int J, int I,const char *aIni);
template<class T> void Drray3(T ***&a, int K, int J, int I);

template<class T> void Array4(T ****&a, int L, int K, int J, int I); // 4-d array
template<class T> void Array4(T ****&a, int L, int K, int J, int I,const char *aIni);
template<class T> void Drray4(T ****&a, int L, int K, int J, int I);

template<class T> void Array5(T *****&a, int M, int L, int K, int J, int I); // 5-d array
template<class T> void Drray5(T *****&a, int M, int L, int K, int J, int I);


// [2.2] MANIPULATOR OF ARRAY/MATRIX/VALUE (templates)
template<class T> void Swap(T &e1, T &e2); // swap two values between variables
template<class T> void Swap(T *e1Ptr,T *e2Ptr); // swap two elements with pointers (RARELY used)
template<class T> int idxMinVal(T *a,int s);
template<class T> int idxMaxVal(T *a,int s);
template<class T> int idxMinValPrtl(T *a,int s1,int s2);
template<class T> int idxMaxValPrtl(T *a,int s1,int s2);
template<class T> void indexMin2Values(T *a,int s,int &iMin1,int &iMin2);
template<class T> void indexMax2Values(T *a,int s,int &iMax1,int &iMax2);

template<class T> void passArray(T *aIn,T *aOut, int s1); //pass 1-d array aIn[] to array aOut[]
template<class T> void passArrayClean(T *a1,T *a2, int s1,int s2); //pass a1[s1] & clean a2[~s2]
template<class T> void passArray_IP(T *aIn,T *aOut, int,int);//pass aIn[] to aOut[] w/ iptIdxOnPartialRg
template<class T> void passArray_OP(T *aIn,T *aOut, int,int);//pass aIn[] to aOut[] w/ outIdxOnPartialRg
template<class T> void passArray_IOP(T *,T *,int,int);//pass aIn[] to aOut[] w/ iptAndOutIdxOnPartialRg
template<class T> void passArray_IPOP(T *,T *,int,int,int);//pass a[] to b[] w/ iptParIdxToOutParIdx
template<class T> void passArray_rdSubset(T *a,int s1, T *b,int s2); // pass rdSubsetAr fm a[] to b[]
template<class T> void passArray_rdSubseg(T *a,int s1, T *b,int s2); // pass rdSubsegmentAr--faster
template<class T> void passArray_lgst(T *a,int s1, T *b,int s2); // pass lgstValues fm a[] to b[]

template<class T> void passArray(T **a,T **b,int s1,int s2); //pass 2-d array
template<class T> void passArrayClean(T **a,T **b,int s1,int s2,int s2Ma); //pass a[][s2]+clean b[][~s2Ma]
template<class T> void passArray(T **aIn,T **aOut,int s1,int *s2);
template<class T> void passArray_I2P(T **,T **,int s1,int s2Min,int s2Max);//pass 2-d w/ 2ndIptIdxPartial
template<class T> void passArray_O1P(T **,T **,int s1Min,int s1Max,int s2);//pass 2-d w/ 1stOutIdxPartial
template<class T> void passArray_O2P(T **,T **,int s1,int s2Min,int s2Max);//pass 2-d w/ 2ndOutIdxPartial
template<class T> void passArray_IO2P(T **,T **,int s1,int s2Min,int s2Max);//pass 2-d 2ndIptOutIdxPartial
template<class T> void passArray(T ***aIn,T ***aOut,int s1,int s2,int s3); //pass 3-d array
template<class T> void passArray_I3P(T ***aIn,T ***aOut,int s1,int s2,int s3Min,int s3Max);
template<class T> void passArray_O3P(T ***aIn,T ***aOut,int s1,int s2,int s3Min,int s3Max);
template<class T> void passArray(T ****aIn,T ****aOut,int s1,int s2,int s3,int s4); //pass 4-d array
template<class T> void passArray(T ****aIn,T ****aOut,int s1,int *s2,int s3,int s4);

template<class T> void asgnArrayConst(T *a, int s, T b); // assign array w/ a constant for all elements
template<class T> void asgnArrayConst(T **a,int s1,int s2,T b);
template<class T> void asgnArrayConst(T **a,int s1,int *s2,T b);
template<class T> void asgnArrayConst(T ***a,int s1,int s2,int s3,T b);

template<class T> void sortArrayAsc(T *a, int s); // sort an array into ascending order
template<class T> void sortArrayDes(T *a, int s); // sort an array into descending order
template<class T> void sortArrayAsc(T **a,int s1, int s2);
template<class T> void sortArrayAsc(T **a,int s1,int *s2);

template<class T> int getRandIndex_RangeArray(T *a, int s, T mi, T ma);//get rdIdxWiRgOfArray
template<class T> void permuteArray(T *a,T *aPerm, int s1); //permute array a[s1] to array aPerm[s1]
template<class T> void permuteArray(T **a,T **b,int s1,int *s2);//permute a[s1][s2[s1]] to b[s1][s2[s1]]
template<class T> int compareArray(T *a,T *b,int s); // compare 2 arrays & return 0=same, 1=different


// [2.3] NUMERICAL STATS (COUNT/SUM/MEAN/SD/COV/CORR/MIN/MAX/MEDIAN/PERCENTILE) (templates)
// [2.3a] COUNT & SIZE (number of elements) OF AN ARRAY
template<class T> int Size(T **a,int s1,int *s2);
template<class T> int CountEq(T *a,int s1, T Val);
template<class T> int CountEq(T **a,int s1,int s2, T Val);
template<class T> int CountEq(T ***a,int s1,int s2,int s3, T Val);
template<class T> int CountMinMax(T *a,int s,T mi,T ma);
template<class T> int CountParMinexMax(T *a,int sF,int sL,T miEx,T ma);
template<class T> int CountGe(T *a,int s1, T Val);
template<class T> int CountGe(T **a,int s1,int s2, T Val);
template<class T> int CountLe(T *a,int s1, T Val);
template<class T> int CountLe(T **a,int s1,int s2, T Val);
template<class T> int CountEqMin(T *a, int s);
template<class T> int CountEqMax(T *a, int s);
template<class T> int CountFix2(T **a,int s1,int jFix, T Val);
template<class T> int CountMinMax(T **a,int s1,int s2, T mi,T ma);

// [2.3b] SUMMATION
template<class T> T Sum(T *a, int s);
template<class T> T SumPar(T *a, int sF, int sL);
template<class T> T Sum_MinMax(T *a, int s, T mi, T ma);
template<class T> double SumSqrt(T *a,int s); // sum of sqrt
template<class T> double SumSqrtPar(T *a,int sF,int sL); // sumOfSqrtBtwIdx SF to sL
template<class T1, class T2> T1 SumWt(T1 *a, T2 *fr, int s); // weighted sum
template<class T> T Sum(T **a,int s1, int s2);
template<class T> T Sum(T **a,int s1,int *s2);
template<class T> T Sum2Fix(T **a,int s1, int j);
template<class T> T Sum1Par2Fix(T **a, int sF, int sL, int j);
template<class T> T Sum1Par(T **a, int sF, int sL, int s2);
template<class T> T Sum1Par(T **a, int sF, int sL,int *s2);
template<class T> T Sum(T ***a,int s1,int s2,int s3);
template<class T> T Sum(T ****a,int s1,int s2,int s3,int s4);
template<class T> void Sum123(T ****a,int s1,int s2,int s3,int s4,T *sum);

// [2.3c] MEAN, VAR, SD, COV & CORR, SKEW, KURT
template<class T> double Mean(T a, T b);
template<class T> double Mean(T a, T b, T c);
template<class T> double Mean(T a, T b, T c, T d);

// 2.3c.1: Mean of a[], a[][], etc
template<class T> double Mean(T *a,int s);
template<class T> double Mean(T **a,int s1, int s2);
template<class T> double Mean(T ***a,int s1,int s2,int s3);
template<class T> double Mean(T ****a,int s1,int s2,int s3,int s4);

// 2.3c.2: Var/SD of a[], a[][], etc
template<class T> double Var(T *a,int s);
template<class T> double SD(T *a,int s);
template<class T> double SD(T **a,int s1,int s2);
template<class T> double SD(T ***a,int s1,int s2,int s3);
template<class T> double SD(T ****a,int s1,int s2,int s3,int s4);

// 2.3c.3: MeanPar of a[], a[][], etc
template<class T> double MeanPar(T *a,int sF,int sL);
template<class T> double MeanPar(T **a,int s1F,int s1L,int s2F,int s2L);
template<class T> double MeanPar(T ***a,int s1F,int s1L,int s2F,int s2L,int s3F,int s3L);

// 2.3c.4: SDPar of a[], a[][], etc
template<class T> double SDPar(T *a,int sF,int sL);
template<class T> double SDPar(T **a,int s1F,int s1L,int s2F,int s2L);
template<class T> double SDPar(T ***a,int s1F,int s1L,int s2F,int s2L,int s3F,int s3L);

// 2.3c.5: MedianPar of a[], a[][], etc
template<class T> double MedianPar(T *a,int s1F,int s1L);
template<class T> double MedianPar(T **a,int s1F,int s1L,int s2F,int s2L);
template<class T> double MedianPar(T ***a,int s1F,int s1L,int s2F,int s2L,int s3F,int s3L);


template<class T> double Mean_Min(T *a, int s,T aMi);
template<class T> double SS(T *a, int s);
template<class T> double SSE(T *a, int s);
template<class T> double VarSk(T *a, int s); // variance of skewDist
template<class T> double SD_Min(T *a, int s,T aMi);

template<class T> double Mean_Wt(T *a,T *wt, int s);
template<class T> double SD_Wt(T *a,T *wt, int s);

template<class T1,class T2> double MeanWt(T1 *a, T2 *fr, int s);

template<class T> double Mean_Min(T **a,int s1,int s2,T aMi);
template<class T> double SD_Min(T **a,int s1,int s2,T aMi);

template<class T> double Mean2Fix(T **a,int s1, int j);
template<class T> double SD2Fix(T **a,int s1, int j);
template<class T> double Mean1Partial2Fix(T **a, int sF, int sL, int j);
template<class T> double SD1Partial2Fix(T **a, int sF, int sL, int j);
template<class T> void MeanSDEtc_UpRtNoDiag(T **a, int s, int &n, double &mean,double &sd, T &min, T &max);
template<class T> void MeanSDEtc_Diag1Up(T **a, int s, int &n, double &mean,double &sd, T &min, T &max);
template<class T> double MeanDiag1Up(T **a, int s);
template<class T> double SDDiag1Up(T **a, int s);


template<class T> double Cov(T *a,T *b, int s);
template<class T> double Corr(T *a,T *b, int s);

template<class T> double Skew(T *a, int s);
template<class T> double Skew_Min(T *a,int s,T aMi);
template<class T> double Kurt(T *a, int s);

template<class T> double Cmt1(T *a,int n); // get k1=cumulant1
template<class T> double Cmt2(T *a,int n); // get k2=cumulant2
template<class T> double Cmt3(T *a,int n); // get k3=cumulant3
template<class T> double Cmt4(T *a,int n); // get k4=cumulant4
void getkSdk1234_dfNc(double df,double nc, double *k,double *sdk);


// [2.3d] MIN, MAX, MEDIAN AND PERCENTILES
// minimums, maximums
template<class T> T Min(T a, T b);
template<class T> T Max(T a, T b);
template<class T> T Min(T a, T b, T c);
template<class T> T Max(T a, T b, T c);
template<class T> T Min(T a, T b, T c, T d);
template<class T> T Max(T a, T b, T c, T d);
template<class T> T Min(T a,T b,T c,T d,T e);
template<class T> T Max(T a,T b,T c,T d,T e);
template<class T> T Min(T a,T b,T c,T d,T e,T f);
template<class T> T Max(T a,T b,T c,T d,T e,T f);
template<class T> T Min(T *a, int s);
template<class T> T Max(T *a, int s);
template<class T> T MinPartial(T *a, int sMin, int sMax);
template<class T> T MaxPartial(T *a, int sMin, int sMax);
template<class T> T MinNth(T *xs,int n, int nth); // get nth min in xs[n]
template<class T> T MaxNth(T *xs,int n, int nth); // get nth max in xs[n]

template<class T> T Min(T **a,int s1, int s2);
template<class T> T Max(T **a,int s1, int s2);
template<class T> T Min(T **a,int s1,int *s2);
template<class T> T Max(T **a,int s1,int *s2);
template<class T> T Min2Fix(T **a,int s1, int j);
template<class T> T Max2Fix(T **a,int s1, int j);
template<class T> T Max1Sum2(T **a,int s1, int s2);

template<class T> T Min(T ***a,int s1,int s2,int s3);
template<class T> T Max(T ***a,int s1,int s2,int s3);
template<class T> T Min(T ****a,int s1,int s2,int s3,int s4);
template<class T> T Max(T ****a,int s1,int s2,int s3,int s4);
template<class T> T Min(T ****a,int s1,int *s2,int s3,int s4);
template<class T> T Max(T ****a,int s1,int *s2,int s3,int s4);

// medians, percentage+rank, percentile
template<class T> double Median(T *a, int s);
template<class T> double Q1(T *a, int s);
template<class T> double Q3(T *a, int s);



// [2.4] CALCULATIONS AND UTILITIES (templates)
// [2.4a] DIGIT/DICIMAL
template<class T> int digit(T v); // get number of digits of an integer or double
template<class T> int digit(T *v,int s1);
template<class T> int digit(T **v,int s1,int s2);
template<class T> int digit(T **v,int s1,int *s2);

int decimal(double d); // get number of decimals of any double
int decimal(double *d,int s1);
int decimal(double **d,int s1,int s2);
int decimal(double **d,int s1,int *s2);

int width(double d); // get width of any double

template<class T> int decnnz(T d,int nnz); // get number of decimals allowing n- non-zeros
template<class T> int decnnz(T *d,int s1,int nnz);
template<class T> int decnnz(T **d,int s1,int s2,int nnz);

template<class T> T Floor(T a,double unit);// get floor+ceil+round of int/dbl
template<class T> int Floor(T a);
template<class T> T Ceil(T a,double unit);
template<class T> int Ceil(T a);
template<class T> T Round(T a,double unit);
template<class T> int Round(T a);

// [2.4b] CHECK RANGE (INTEGER OR DOUBLE)
template<class T> void chkVal(T v, char *vNm, T val); // check a value (int/dbl)
template<class T> void chkV2(T v, char *vNm, T v1,T v2);// check 2 values (int/dbl)
template<class T> void chkRng(T v, char *vNm, T vmin, T vmax);// chk range of int/dbl
template<class T> void chkRng(T *v,int vsz, char *vNm, T vmin, T vmax);// check range of an array
template<class T> void chkMin(T v, char *vNm, T vmin);// check min of a value
template<class T> void chkMin(T *v,int vsz, char *vNm, T vmin);// check min of an array
template<class T> void chkMax(T v, char *vNm, T vmax);// check max of a value

// [2.4c] PRINT ARRAY
template<class T> void prtAr1(ostream &out,T *a,int s1);             // prt 1-d array w/ rtn

template<class T> void getWSetDec_Nnz(ostream &out,T *a,int s1,int nnz, int &w); // getW+setDec co Nnz
template<class T> void getWSetDec_Nnz(ostream &out,T **a,int s1,int s2,int nnz, int &w);
template<class T> void prtArray1(ostream &out,T *a,int s1,char *rt); // prt 1-d array w/ rtn & fix2nz

void prtArray1(ostream &out,char **a,int s1,char *rt);               // prt 1-d char array w/ or wo rtn
void prtArray2(ostream &out,char ***a,int s1,int s2);                // prt 2-d char array w/ rtn
void prtArray3(ostream &out,char ****a,int s1,int s2,int s3);        // prt 3-d char array w/ rtn
void prtArray3(ostream &out,char ****a,int s1,int *s2,int s3);       // prt 3-d char array w/ rtn

template<class T> void prtArray1_P(ostream &out,T *a,int sF,int sL,char *rt);// prt1da parIdx sF-sL wRtn
template<class T> void prtArray1BRNR(ostream &out,T *a,int s1, int nPR);// prt1da byRowWoRtn
template<class T> void prtArray1BR(ostream &out,T *a,int s1, int nPR);  // prt1da byRowWRtn
template<class T> void prtArray1FWNR(ostream &out,T *a,int s1, int w);  // prt1da w/ fixWidthNR
template<class T> void prtArray1FW(ostream &out,T *a,int s1, int w);    // prt1da w/ fixWidth
template<class T> void prtArray1FDNR(ostream &out,T *a,int s1, int dec);// prt1da w/ fixDecNR
template<class T> void prtArray1FD(ostream &out,T *a,int s1, int dec);  // prt1da w/ fixDec
template<class T> void prtArray1FWD(ostream &out,T *a,int s1, int w,int dec);
template<class T> void prtArray1FWBRNR(ostream &out,T *a,int s1, int w,int nPR);
template<class T> void prtArray1FWBR(ostream &out,T *a,int s1, int w,int nPR);
template<class T> void prtArray1FWDBRNR(ostream &out,T *a,int s1, int w,int dec,int nPR);
template<class T> void prtArray1FWDBR(ostream &out,T *a,int s1, int w,int dec,int nPR);
template<class T> void prtArray1DlmNR(ostream &out,T *a,int s1, char *c);//prt 1-d WcharDlmWoRtn
template<class T> void prtArray1Dlm(ostream &out,T *a,int s1, char *c);  //prt 1-d WcharDlmWRtn
template<class T> void prtArray1DlmBRNR(ostream &out,T *a,int s1, char *c,int nPR);
template<class T> void prtArray1DlmBR(ostream &out,T *a,int s1, char *c,int nPR);

template<class T> void prtArray1BRTF(T *a,int s1,int nPR, char *fn); // prt 1-d array to file

template<class T> void prtArray2(ostream &out,T **a,int s1,int s2, int nnz);// prt 2-d array w/ nnz
template<class T> void prtArray2(ostream &out,T **a,int s1,int s2);
template<class T> void prtArray2FWD(ostream &out,T **a,int s1,int s2, int w,int dec);
template<class T> void prtMtx_UpRtNoDiag(ostream &out,T **a,int s1); // prt mtx up rt w/o diag
template<class T> void prtArray2_2P(ostream &out,T **a,int s1,int s2Min,int s2Max);

template<class T> void prtArray2_idCsv(char *fn,int *idx,T **a,int,int,int nnz); // prt idx+2dAr inCsvFmt
template<class T> void prtArray2_idCsv(char *fn,int *idx,T **a,int s1,int s2);

template<class T> void prtArray3(ostream &out,T ***a,int s1,int s2,int s3);
template<class T> void prtArray3_1P(ostream &out,T ***a,int s1Min,int s1Max,int s2,int s3);
template<class T> void prtArray3_3P(ostream &out,T ***a,int s1,int s2,int s3Min,int s3Max);
template<class T> void prtArray4(ostream &out,T ****a,int s1,int s2,int s3,int s4);
template<class T> void prtArray4(ostream &out,T ****a,int s1,int *s2,int s3,int s4);
template<class T> void prtArray5(ostream &out,T *****a,int s1,int s2,int s3,int s4,int s5);

// [2.4d] PRINT OR GET STATS nGroup+Freq+MeanSd
template<class T> int nGroup(T *a, int s);// get number of groups in 1-d array
template<class T> int nGroup(T *a,int s, T mi,T ma);// get nOfGps wi-a-rg in 1-d array
template<class T> int nGroup(T **a,int s1, int s2);// get number of groups in 2-d array
template<class T> void Freq(T *a, int s,T *gVal, int *gCnt);// get freq (output=gVal[], gCnt[])
template<class T> void Freq(T **a,int s1,int s2,T *gVal, int *gCnt);// get freq in 2-d array
template<class T> void prtFreq(ostream &out,T *a, int s);// print frequency in 1-d array
template<class T> void prtFreqHet(ostream &out,T *a, int s);
template<class T> void prtFreq(ostream &out,T **a,int s1,int s2);// print frequency in 2-d array
template<class T> void MeanSDByGroup(T *a,T *aGroup, int s,T *gVal,
   int *gCnt, double *gMean, double *gSD,T *gMin,T *gMax);// get mean and sd by group in 1-d array
template<class T> void MeanSDByGroup(T **a,T **aGroup,int s1,int s2,T *gVal,
   int *gCnt, double *gMean, double *gSD);// get mean and sd by group in 2-d array
template<class T> void prtMeanSDByGroup(ostream &out,T *a,T *aGroup, int s);// print mean+sd by group
template<class T> void prtMeanSDByGroup(ostream &out,T **,T **,int,int);//prtMean+sd by gp in 2-d array


//--------------------------------------------------------------------
// [2.1] MEMORY ALLOCATION OF ARRAYS (templates)
//--------------------------------------------------------------------

// 1-d array
template<class T> void Array1(T *&a, int I)
{
   int i; T aIni=(T)(-1);
   a=new T[I];
   for (i=0; i<I; i++) a[i]=aIni;
}
template<class T> void Array1(T *&a, int I, T aIni)
{
   int i;
   a=new T[I];
   for (i=0;i<I;i++) a[i]=aIni;
}
void Array1(char *&a, int I,const char *aIni)
{
   a=new char[I];
   strcpy(a,aIni);
}
void Array1(string *&a, int I,const char *aIni)
{
   int i;
   a=new string[I];
   for (i=0;i<I;i++) a[i]=aIni;
}

// 1-d array expand from size Ipre to Inew
template<class T> void A1epd(T *&a, int Ipre,int Inew)
{
   int i; T *aKp,aIni=(T)(-1);

   if (Inew<=Ipre) { cerr<<"\nERROR: The expanded array size ("<<Inew
      <<") is not > the previous one ("<<Ipre<<").\n\n"; exit(1); }
   Array1(aKp,Ipre); passArray(a,aKp,Ipre);
   a=new T[Inew];
   for (i=0;i<Inew;i++) { if (i<Ipre) a[i]=aKp[i]; else a[i]=aIni; }

   Drray1(aKp,Ipre);
}

template<class T> void Drray1(T *&a, int I)
{
   if (a!=0 && I>=1) { delete [] a; }
}

// 2-d array
template<class T> void Array2(T **&a, int J, int I)
{
   int j, i; T aIni=(T)(-1);
   a=new T*[J]; a[0]=new T[J*I];
   for (j=0; j<J; j++)
   {
      if (j>=1) a[j]=a[0]+(j*I);
      for (i=0; i<I; i++) a[j][i]=aIni;
   }
}
template<class T> void Array2(T **&a, int J, int I,T aIni)
{
   int j, i;
   a=new T*[J]; a[0]=new T[J*I];
   for (j=0; j<J; j++)
   {
      if (j>=1) a[j]=a[0]+(j*I);
      for (i=0; i<I; i++) a[j][i]=aIni;
   }
}
void Array2(char **&a, int J, int I,const char *aIni)
{
   int j;
   a=new char*[J]; a[0]=new char[J*I];
   for (j=0; j<J; j++)
   {
      if (j>=1) a[j]=a[0]+(j*I);
      strcpy(a[j],aIni);
   }
}
void Array2(string **&a,int J,int I,string aIni)
{
   int j, i;
   a=new string*[J]; a[0]=new string[J*I];
   for (j=0; j<J; j++)
   {
      if (j>=1) a[j]=a[0]+(j*I);
      for (i=0; i<I; i++) a[j][i]=aIni;
   }
}

// 2-d array expand of size2 from Ipre to Inew
template<class T> void A2ep2(T **&a, int J, int Ipre,int Inew)
{
   int j,i; T **aKp,aIni=(T)(-1);

   if (Inew<=Ipre) { cerr<<"\nERROR: The expanded array size ("<<Inew
      <<") is not > the previous one ("<<Ipre<<").\n\n"; exit(1); }
   Array2(aKp,J,Ipre); passArray(a,aKp,J,Ipre);
   a=new T*[J]; a[0]=new T[J*Inew];
   for (j=0;j<J;j++)
   {
      if (j>=1) a[j]=a[0]+(j*Inew);
      for (i=0; i<Inew; i++) { if (i<Ipre) a[j][i]=aKp[j][i]; else a[j][i]=aIni; }
   }
}

template<class T> void Drray2(T **&a, int J, int I)
{
   if (a!=0 && J>=1 && I>=1)
   {
      delete [] a[0];
      delete [] a;
   }
}

// 3-d array
template<class T> void Array3(T ***&a, int K, int J, int I)
{
   int k, j, i; T aIni=-1;
   a=new T**[K]; a[0]=new T*[K*J]; a[0][0]=new T[K*J*I];
   for (k=0; k<K; k++)
   {
      if (k>=1) a[k]=a[0]+(k*J);
      for (j=0; j<J; j++)
      {
         if (k>=1||j>=1) a[k][j]=a[0][0]+(k*J*I+j*I);
         for (i=0; i<I; i++) a[k][j][i]=aIni;
      }
   }
}
template<class T> void Array3(T ***&a, int K, int J, int I,const char *aIni)
{
   int k, j;
   a=new T**[K]; a[0]=new T*[K*J]; a[0][0]=new T[K*J*I];
   for (k=0; k<K; k++)
   {
      if (k>=1) a[k]=a[0]+(k*J);
      for (j=0; j<J; j++)
      {
         if (k>=1||j>=1) a[k][j]=a[0][0]+(k*J*I+j*I);
         strcpy(a[k][j],aIni);
      }
   }
}

template<class T> void Drray3(T ***&a, int K, int J, int I)
{
   if (a!=0 && K>=1 && J>=1 && I>=1) { delete [] a[0][0]; delete [] a[0]; delete [] a; }
}

// 4-d array
template<class T> void Array4(T ****&a, int L, int K, int J, int I)
{
   int l,k,j,i; T aIni=-1;
   a=new T***[L]; a[0]=new T**[L*K]; a[0][0]=new T*[L*K*J]; a[0][0][0]=new T[L*K*J*I];
   for (l=0;l<L;l++)
   {
      if (l>=1) a[l]=a[0]+(l*K);
      for (k=0; k<K; k++)
      {
         if (l>=1||k>=1) a[l][k]=a[0][0]+(l*K*J+k*J);
         for (j=0; j<J; j++)
         {
            if (l>=1||k>=1||j>=1) a[l][k][j]=a[0][0][0]+(l*K*J*I+k*J*I+j*I);
            for (i=0; i<I; i++) a[l][k][j][i]=aIni;
         }
      }
   }
}
template<class T> void Array4(T ****&a, int L, int K, int J, int I,const char *aIni)
{
   int l,k,j;
   a=new T***[L]; a[0]=new T**[L*K]; a[0][0]=new T*[L*K*J]; a[0][0][0]=new T[L*K*J*I];
   for (l=0;l<L;l++)
   {
      if (l>=1) a[l]=a[0]+(l*K);
      for (k=0; k<K; k++)
      {
         if (l>=1||k>=1) a[l][k]=a[0][0]+(l*K*J+k*J);
         for (j=0; j<J; j++)
         {
            if (l>=1||k>=1||j>=1) a[l][k][j]=a[0][0][0]+(l*K*J*I+k*J*I+j*I);
            strcpy(a[l][k][j],aIni);
         }
      }
   }
}

template<class T> void Drray4(T ****&a, int L, int K, int J, int I)
{
   if (a!=0 && L>=1 && K>=1 && J>=1 && I>=1)
   { delete [] a[0][0][0]; delete [] a[0][0]; delete [] a[0]; delete [] a; }
}

// 5-d array
template<class T> void Array5(T *****&a, int M, int L, int K, int J, int I)
{
   int m, l, k, j, i; T aIni=-1;
   a=new T****[M]; a[0]=new T***[M*L]; a[0][0]=new T**[M*L*K]; a[0][0][0]=new T*[M*L*K*J];
   a[0][0][0][0]=new T[M*L*K*J*I];
   for (m=0; m<M; m++)
   {
      if (m>=1) a[m]=a[0]+(m*L);
      for (l=0;l<L;l++)
      {
         if (m>=1||l>=1) a[m][l]=a[0][0]+(m*L*K+l*K);
         for (k=0; k<K; k++)
         {
            if (m>=1||l>=1||k>=1) a[m][l][k]=a[0][0][0]+(m*L*K*J+l*K*J+k*J);
            for (j=0; j<J; j++)
            {
               if (m>=1||l>=1||k>=1||j>=1) a[m][l][k][j]=a[0][0][0][0]+(m*L*K*J*I+l*K*J*I+k*J*I+j*I);
               for (i=0; i<I; i++) a[m][l][k][j][i]=aIni;
            }
         }
      }
   }
}
template<class T> void Drray5(T *****&a, int M, int L, int K, int J, int I)
{
   if (a!=0 && M>=1 && L>=1 && K>=1 && J>=1 && I>=1)
   { delete[] a[0][0][0][0]; delete[] a[0][0][0]; delete[] a[0][0]; delete[] a[0]; delete[] a; }
}


//--------------------------------------------------------------------
// [2.2] MANIPULATOR OF ARRAY/MATRIX/VALUE (templates)
//--------------------------------------------------------------------

// swap two values between variables
template<class T> void Swap(T &e1, T &e2)
{
   T hold = e1; e1 = e2; e2 = hold;
}

// swap two elements with pointers (RARELY used)
template<class T> void Swap(T *e1Ptr,T *e2Ptr)
{
   T hold = *e1Ptr; *e1Ptr = *e2Ptr; *e2Ptr = hold;
}

template<class T> int idxMinVal(T *a,int s)
{
   T min; int i,nMin,iMinRtn;
   T *aMin;   Array1(aMin,s);
   int *iMin; Array1(iMin,s);

   min=Min(a,s); nMin=0;
   for (i=0;i<s;i++) { if (a[i]==min) { aMin[nMin]=a[i]; iMin[nMin]=i; nMin++; } }
   i=RdInt(0,nMin-1); iMinRtn=iMin[i];
   Drray1(aMin,s);
   Drray1(iMin,s);
   return iMinRtn;
}

template<class T> int idxMaxVal(T *a,int s)
{
   T max; int i,nMax,iMaxRtn;
   T *aMax;   Array1(aMax,s);
   int *iMax; Array1(iMax,s);

   max=Max(a,s); nMax=0;
   for (i=0;i<s;i++) { if (a[i]==max) { aMax[nMax]=a[i]; iMax[nMax]=i; nMax++; } }
   i=RdInt(0,nMax-1); iMaxRtn=iMax[i];
   Drray1(aMax,s);
   Drray1(iMax,s);
   return iMaxRtn;
}

template<class T> int idxMinValPrtl(T *a,int s1,int s2)
{
   T min; int i,nMin,iMinRtn;
   T *aMin;   Array1(aMin,s2+1);
   int *iMin; Array1(iMin,s2+1);

   min=MinPartial(a,s1,s2); nMin=0;
   for (i=s1;i<=s2;i++) { if (a[i]==min) { aMin[nMin]=a[i]; iMin[nMin]=i; nMin++; } }
   i=RdInt(0,nMin-1); iMinRtn=iMin[i];
   Drray1(aMin,s2+1);
   Drray1(iMin,s2+1);
   return iMinRtn;
}

template<class T> int idxMaxValPrtl(T *a,int s1,int s2)
{
   T max; int i,nMax,iMaxRtn;
   T *aMax;   Array1(aMax,s2+1);
   int *iMax; Array1(iMax,s2+1);

   max=MaxPartial(a,s1,s2); nMax=0;
   for (i=s1;i<=s2;i++) { if (a[i]==max) { aMax[nMax]=a[i]; iMax[nMax]=i; nMax++; } }
   i=RdInt(0,nMax-1); iMaxRtn=iMax[i];
   Drray1(aMax,s2+1);
   Drray1(iMax,s2+1);
   return iMaxRtn;
}

template<class T> void indexMin2Values(T *a,int s,int &iMin1,int &iMin2)
{
   T min; int i,i1,i2,nMin;
   T *aMin;   Array1(aMin,s);
   int *iMin; Array1(iMin,s);

   min=Min(a,s); nMin=0;
   for (i=0;i<s;i++) { if (a[i]==min) { aMin[nMin]=a[i]; iMin[nMin]=i; nMin++; } }
   if (nMin>=2)
   {
      get2RdIntNR(0,nMin-1, i1,i2); iMin1=iMin[i1]; iMin2=iMin[i2];
   }
   else if (nMin==1)
   {
      iMin1=iMin[0];
      for (i=0;i<s;i++) { if (i<iMin1) aMin[i]=a[i]; else if (i>iMin1) aMin[i-1]=a[i]; }
      iMin2=idxMinVal(aMin,s-1); if (iMin2>=iMin1) iMin2++;
   }
}

template<class T> void indexMax2Values(T *a,int s,int &iMax1,int &iMax2)
{
   T max; int i,i1,i2,nMax;
   T *aMax;   Array1(aMax,s);
   int *iMax; Array1(iMax,s);

   max=Max(a,s); nMax=0;
   for (i=0;i<s;i++) { if (a[i]==max) { aMax[nMax]=a[i]; iMax[nMax]=i; nMax++; } }
   if (nMax>=2)
   {
      get2RdIntNR(0,nMax-1, i1,i2); iMax1=iMax[i1]; iMax2=iMax[i2];
   }
   else if (nMax==1)
   {
      iMax1=iMax[0];
      for (i=0;i<s;i++) { if (i<iMax1) aMax[i]=a[i]; else if (i>iMax1) aMax[i-1]=a[i]; }
      iMax2=idxMaxVal(aMax,s-1); if (iMax2>=iMax1) iMax2++;
   }
}

// pass array aIn[s1] to array aOut[s1]
template<class T> void passArray(T *aIn,T *aOut, int s1)
{
   for (int i=0; i<s1; i++) aOut[i]=aIn[i];
}
// pass array a1[s1] to array a2[s1], and clean a2[]:i=s1~s2-1
template<class T> void passArrayClean(T *a1,T *a2, int s1,int s2)
{
   int i,sMi,sMa; sMi=Min(s1,s2); sMa=Max(s1,s2);
   for (i=0;i<sMa;i++) a2[i]=(i<sMi?a1[i]:-1);
}

// pass array aIn[n] to array aOut[n] w/ Input index on Partial range and output from [0]
template<class T> void passArray_IP(T *aIn,T *aOut, int sMin,int sMax)
{
   for (int i=sMin;i<=sMax;i++) aOut[i-sMin]=aIn[i];
}
// pass array aIn[n] to array aOut[n] w/ Output index on Partial range and input from [0]
template<class T> void passArray_OP(T *aIn,T *aOut, int sMin,int sMax)
{
   for (int i=sMin;i<=sMax;i++) aOut[i]=aIn[i-sMin];
}
// pass array aIn[n] to array aOut[n] w/ both Input & Output indexes on Partial ranges
template<class T> void passArray_IOP(T *aIn,T *aOut, int sMin,int sMax)
{
   for (int i=sMin;i<=sMax;i++) aOut[i]=aIn[i];
}
// pass a[] to b[] w/ both input idx partial to output idx partial
template<class T> void passArray_IPOP(T *a,T *b,int saMi,int saMa,int sbMi)
{
   for (int i=saMi;i<=saMa;i++) b[i-saMi+sbMi]=a[i];
}

// pass random-subset-array from a[s1] to b[s2] co s2<=s1
template<class T> void passArray_rdSubset(T *a,int s1, T *b,int s2)
{
   int i,*iPm; Array1(iPm,s1);

   if (s2>s1) { cerr<<"\nERROR: passArray_rdSubset(), invalid array sizes.\n\n"; exit(1); }
   asgnArray012Perm(iPm,s1);
   for (i=0;i<s2;i++) b[i]=a[iPm[i]];

   Drray1(iPm,s1);
}
// pass random-subsegment-array from a[s1] to b[s2] co s2<=s1
template<class T> void passArray_rdSubseg(T *a,int s1, T *b,int s2)
{
   int i,iRd;

   if (s2>s1) { cerr<<"\nERROR: passArray_rdSubseg(), invalid array sizes.\n\n"; exit(1); }
   iRd=RdInt(0,s1-s2); for (i=0;i<s2;i++) b[i]=a[iRd+i];
}

// pass s2 largest values from a[s1] to b[s2] co s2<=s1
template<class T> void passArray_lgst(T *a,int s1, T *b,int s2)
{
   int i,j,j1; T ami,a1;

   ami=Min(a,s1)-1; for (j=0;j<s2;j++) b[j]=ami;
   for (i=0;i<s1;i++)
   {
      a1=a[i];
      if (a1>b[0])
      {
         for (j=0;j<s2;j++)
         {
            if ((j<=s2-2 && a1>b[j] && a1<=b[j+1]) || (j==s2-1 && a1>b[s2-1]))
            {
               for (j1=0;j1<j;j1++) b[j1]=b[j1+1]; b[j]=a1; break;
            }
         }
      }
   }
}


//pass 2-d array a[s1][s2] to b[][]
template<class T> void passArray(T **a,T **b,int s1,int s2)
{
   int i;
   for (i=0;i<s1;i++) for (int j=0;j<s2;j++) b[i][j]=a[i][j];
}
// pass 2-d a[s1][s2] to b[][] & clean b[s1][~s2Ma]
template<class T> void passArrayClean(T **a,T **b,int s1,int s2,int s2Ma)
{
   int i;
   for (i=0;i<s1;i++) for (int j=0;j<s2Ma;j++) b[i][j]=(j<s2?a[i][j]:-1);
}

template<class T> void passArray(T **aIn,T **aOut,int s1,int *s2)
{
   for (int i=0;i<s1;i++) for (int j=0; j<s2[i]; j++) aOut[i][j]=aIn[i][j];
}

template<class T> void passArray_I2P(T **aIn,T **aOut,int s1,int s2Min,int s2Max)
{
   int i,j;
   for (i=0;i<s1;i++) for (j=s2Min;j<=s2Max;j++) aOut[i][j-s2Min]=aIn[i][j];
}

template<class T> void passArray_O1P(T **aIn,T **aOut,int s1Min,int s1Max,int s2)
{
   int i,j;
   for (i=s1Min;i<=s1Max;i++) for (j=0;j<s2;j++) aOut[i][j]=aIn[i-s1Min][j];
}
template<class T> void passArray_O2P(T **aIn,T **aOut,int s1,int s2Min,int s2Max)
{
   int i,j;
   for (i=0;i<s1;i++) for (j=s2Min;j<=s2Max;j++) aOut[i][j]=aIn[i][j-s2Min];
}
template<class T> void passArray_IO2P(T **aIn,T **aOut,int s1,int s2Min,int s2Max)
{
   int i,j;
   for (i=0;i<s1;i++) for (j=s2Min;j<=s2Max;j++) aOut[i][j]=aIn[i][j];
}

template<class T> void passArray(T ***aIn,T ***aOut,int s1,int s2,int s3)
{
   int i,j,k;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++) aOut[i][j][k]=aIn[i][j][k];
}
template<class T> void passArray_I3P(T ***aIn,T ***aOut,int s1,int s2,int s3Min,int s3Max)
{
   int i,j,k;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=s3Min;k<=s3Max;k++) aOut[i][j][k-s3Min]=aIn[i][j][k];
}
template<class T> void passArray_O3P(T ***aIn,T ***aOut,int s1,int s2,int s3Min,int s3Max)
{
   int i,j,k;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=s3Min;k<=s3Max;k++) aOut[i][j][k]=aIn[i][j][k-s3Min];
}

template<class T> void passArray(T ****aIn,T ****aOut,int s1,int s2,int s3,int s4)
{
   int i,j,k,l;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++) for(l=0;l<s4;l++)
      aOut[i][j][k][l]=aIn[i][j][k][l];
}

template<class T> void passArray(T ****aIn,T ****aOut,int s1,int *s2,int s3,int s4)
{
   int i,j,k,l;
   for (i=0;i<s1;i++) for (j=0;j<s2[i];j++) for (k=0;k<s3;k++) for(l=0;l<s4;l++)
      aOut[i][j][k][l]=aIn[i][j][k][l];
}

// assign array w/ a constant for all elements
template<class T> void asgnArrayConst(T *a, int s, T b)
{
   for (int i=0;i<s;i++) a[i]=b;
}
template<class T> void asgnArrayConst(T **a,int s1,int s2,T b)
{
   for (int i=0;i<s1;i++) for (int j=0;j<s2;j++) a[i][j]=b;
}
template<class T> void asgnArrayConst(T **a,int s1,int *s2,T b)
{
   for (int i=0;i<s1;i++) for (int j=0;j<s2[i];j++) a[i][j]=b;
}
template<class T> void asgnArrayConst(T ***a,int s1,int s2,int s3,T b)
{
   for (int i=0;i<s1;i++) for (int j=0;j<s2;j++) for (int k=0;k<s3;k++) a[i][j][k]=b;
}

// sort an array into ascending order
template<class T> void sortArrayAsc(T *a, int s)
{
   for (int pass=1; pass<s; pass++) for (int j=0; j<s-1; j++)
      if (a[j]>a[j+1]) Swap( &a[j], &a[j+1] );
}

// sort an array into descending order
template<class T> void sortArrayDes(T *a, int s)
{
   for (int pass=1; pass<s; pass++) for (int j=0; j<s-1; j++)
      if (a[j]<a[j+1]) Swap( &a[j], &a[j+1] );
}

template<class T> void sortArrayAsc(T **a,int s1, int s2)
{
   int i, j, ij, total=s1*s2;
   T *aJk; Array1(aJk, total);

   ij=0;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) aJk[ij++]=a[i][j];
   sortArrayAsc(aJk, total);

   ij=0;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) a[i][j]=aJk[ij++];

   Drray1(aJk, total);
}

template<class T> void sortArrayAsc(T **a,int s1,int *s2)
{
   int i, j, ij, total; total=Sum(s2, s1);
   T *aJk; Array1(aJk, total);

   ij=0;
   for (i=0;i<s1;i++) for (j=0; j<s2[i]; j++) aJk[ij++]=a[i][j];
   sortArrayAsc(aJk, total);

   ij=0;
   for (i=0;i<s1;i++) for (j=0; j<s2[i]; j++) a[i][j]=aJk[ij++];

   Drray1(aJk, total);
}


// get random index within range values of an array
template<class T> int getRandIndex_RangeArray(T *a, int s, T mi, T ma)
{
   int i,nX,iRX,iRXJk,iRd;

   nX=CountMinMax(a,s,mi,ma);
   iRX=RdInt(0,nX-1); iRXJk=-1; iRd=-1;
   for (i=0;i<s;i++)
   {
      if (a[i]>=mi && a[i]<=ma) iRXJk++;
      if (iRXJk==iRX) { iRd=i; break; }
      if (i==s-1 && iRd==-1)
      {
         cerr<<"\nERROR: No index found within the specified range of array values.\n\n";
         exit(1);
      }
   }
   return iRd;
}

// permute array a[s1] to array aPerm[s1]
template<class T> void permuteArray(T *a,T *aPerm, int s1)
{
   int i,*iPm; Array1(iPm,s1);

   asgnArray012Perm(iPm, s1);
   for (i=0;i<s1;i++) aPerm[i]=a[iPm[i]];

   Drray1(iPm,s1);
}

// permute array a[s1][s2[s1]] to array b[s1][s2[s1]]
template<class T> void permuteArray(T **a,T **b,int s1,int *s2)
{
   int i, j, iJk, ii, jj, iiJk, totSize=Sum(s2, s1);
   int *perm; Array1(perm, totSize);

   asgnArray012Perm(perm, totSize);

   iJk=0;
   for (i=0;i<s1;i++) for (j=0; j<s2[i]; j++)
   {
      iiJk=0;
      for (ii=0; ii<s1; ii++) for (jj=0; jj<s2[ii]; jj++)
      {
         if (iiJk==perm[iJk])
         {
            b[i][j]=a[ii][jj];
            goto next_bij;
         }
         else iiJk++;
      }

next_bij:
      iJk++;
   }

   Drray1(perm, totSize);
}


// compare 2 arrays & return 0=same, 1=different
template<class T> int compareArray(T *a,T *b,int s)
{
   int dif,i; dif=0;
   for (i=0;i<s;i++) { if (a[i]!=b[i]) { dif=1; break; } }
   return dif;
}


///--------------------------------------------------------------------
// [2.3] NUMERICAL STATS (COUNT/SUM/MEAN/SD/COV/CORR/MIN/MAX/MEDIAN/PERCENTILE) (templates)
//--------------------------------------------------------------------

// [2.3a] COUNT & SIZE (number of elements) OF AN ARRAY
template<class T> int Size(T **a,int s1,int *s2)
{
   int size = 0;
   for (int i=0;i<s1;i++) for (int j=0; j<s2[i]; j++) size++;
   return size;
}

template<class T> int CountEq(T *a,int s1, T Val)
{
   int i, n=0;
   for (i=0;i<s1;i++) { if (a[i]==Val) n++; }
   return n;
}
template<class T> int CountEq(T **a,int s1,int s2, T Val)
{
   int i, j, n = 0;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) { if (a[i][j]==Val) n++; }
   return n;
}
template<class T> int CountEq(T ***a,int s1,int s2,int s3, T Val)
{
   int i,j,k,n=0;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++) { if (a[i][j][k]==Val) n++; }
   return n;
}

template<class T> int CountMinMax(T *a, int s, T mi, T ma)
{
   int i, n = 0;
   for (i=0; i<s; i++) { if (a[i]>=mi && a[i]<=ma) n++; }
   return n;
}

// get count of a[]=(miEx, maEq] with partial idx range sL idx1Fix & idx2Par
template<class T> int CountParMinexMax(T *a,int sF,int sL,T miEx,T ma)
{
   int i, n = 0;
   for (i=sF;i<=sL;i++) { if (a[i]>miEx && a[i]<=ma) n++; }
   return n;
}

template<class T> int CountGe(T *a,int s1, T Val)
{
   int i,n;
   n=0; for (i=0;i<s1;i++) { if (a[i]>=Val) n++; }
   return n;
}
template<class T> int CountGe(T **a,int s1,int s2, T Val)
{
   int i,j,n;
   n=0; for (i=0;i<s1;i++) for (j=0;j<s2;j++) { if (a[i][j]>=Val) n++; }
   return n;
}

template<class T> int CountLe(T *a,int s1, T Val)
{
   int i,n;
   n=0; for (i=0;i<s1; i++) { if (a[i]<=Val) n++; }
   return n;
}
template<class T> int CountLe(T **a,int s1,int s2, T Val)
{
   int i,j,n;
   n=0; for (i=0;i<s1; i++) for (j=0;j<s2;j++) { if (a[i][j]<=Val) n++; }
   return n;
}

template<class T> int CountEqMin(T *a, int s)
{
   T min; int i,n;
   min = a[0]; for (i=1;i<s;i++) { if (min>a[i]) min=a[i]; }
   n=0; for (i=0;i<s;i++) { if (a[i]==min) n++; }
   return n;
}

template<class T> int CountEqMax(T *a, int s)
{
   T max; int i,n;
   max = a[0]; for (i=1;i<s;i++) { if (max<a[i]) max=a[i]; }
   n=0; for (i=0;i<s;i++) { if (a[i]==max) n++; }
   return n;
}

template<class T> int CountFix2(T **a,int s1,int jFix, T Val)
{
   int i, n = 0;
   for (i=0;i<s1;i++) { if (a[i][jFix]==Val) n++; }
   return n;
}

template<class T> int CountMinMax(T **a,int s1,int s2, T mi,T ma)
{
   int i, j, n = 0;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++)
   {
      if (a[i][j]>=mi && a[i][j]<=ma) n++;
   }
   return n;
}


// [2.3b] SUMMATION
template<class T> T Sum(T *a, int s)
{
   int i; T sum=0;
   for (i=0;i<s;i++) sum += a[i];
   return sum;
}

template<class T> T SumPar(T *a, int sF, int sL)
{
   int i; T sum=0;
   for (i=sF;i<=sL;i++) sum += a[i];
   return sum;
}

template<class T> T Sum_MinMax(T *a, int s, T mi, T ma)
{
   int i; T sum=0;
   for (i=0;i<s;i++) { if (a[i]>=mi && a[i]<=ma) sum += a[i]; }
   return sum;
}

template<class T> double SumSqrt(T *a, int s) // sum of sqrt
{
   int i; double sum=0;
   for (i=0;i<s;i++) sum += sqrt(a[i]);
   return sum;
}

template<class T> double SumSqrtPar(T *a,int sF,int sL) // sumOfSqrtBtwIdx SF to sL
{
   int i; double sum=0;
   for (i=sF;i<=sL;i++) sum += sqrt(a[i]);
   return sum;
}

template<class T1, class T2> T1 SumWt(T1 *a, T2 *fr, int s) // weighted sum
{
   int i; T1 sum=0;
   for (i=0;i<s;i++) sum += a[i]*fr[i];
   return sum;
}

template<class T> T Sum(T **a,int s1, int s2)
{
   int i,j; T sum=0;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) sum += a[i][j];
   return sum;
}

template<class T> T Sum(T **a,int s1,int *s2)
{
   int i,j; T sum=0;
   for (i=0;i<s1;i++) for (j=0;j<s2[i];j++) sum += a[i][j];
   return sum;
}

template<class T> T Sum2Fix(T **a,int s1, int j)
{
   int i; T sum=0;
   for (i=0;i<s1;i++) sum += a[i][j];
   return sum;
}

template<class T> T Sum1Par2Fix(T **a, int sF, int sL, int j)
{
   int i; T sum=0;
   for (i=sF;i<=sL;i++) sum += a[i][j];
   return sum;
}

template<class T> T Sum1Par(T **a, int sF, int sL, int s2)
{
   int i,j; T sum=0;
   for (i=sF;i<=sL;i++) for (j=0;j<s2;j++) sum += a[i][j];
   return sum;
}

template<class T> T Sum1Par(T **a, int sF, int sL,int *s2)
{
   int i,j; T sum=0;
   for (i=sF;i<=sL;i++) for (j=0;j<s2[i];j++) sum += a[i][j];
   return sum;
}

template<class T> T Sum(T ***a,int s1,int s2,int s3)
{
   int i,j,k; T sum=0;
   for (i=0;i<s1;i++) for (int j=0;j<s2;j++) for (k=0;k<s3;k++) sum += a[i][j][k];
   return sum;
}

template<class T> T Sum(T ****a,int s1,int s2,int s3,int s4)
{
   int i,j,k,l; T sum=0;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++) for (l=0;l<s4;l++)
      sum += a[i][j][k][l];
   return sum;
}

template<class T> void Sum123(T ****a,int s1,int s2,int s3,int s4,T *sum)
{
   for (int l=0;l<s4;l++)
   {
      sum[l] = 0;
      for (int i=0;i<s1;i++) for (int j=0;j<s2;j++) for (int k=0;k<s3;k++)
         sum[l] += a[i][j][k][l];
   }
}


// [2.3c] MEAN, VAR, SD, SKEW, COV & CORR
template<class T> double Mean(T a, T b) { return (double)(a+b)/2; }
template<class T> double Mean(T a, T b, T c) { return (double)(a+b+c)/3; }
template<class T> double Mean(T a, T b, T c, T d) {return (double)(a+b+c+d)/4;}


// 2.3c.1: Mean of a[], a[][], etc
template<class T> double Mean(T *a,int s)
{
   int i; double sum,me;

   sum=0; for (i=0;i<s;i++) sum+=a[i]; me=sum/s;

   return me;
}

template<class T> double Mean(T **a,int s1,int s2)
{
   int i,j; double sum,me;

   sum=0; for (i=0;i<s1;i++) for (j=0;j<s2;j++) sum+=a[i][j]; me=sum/(s1*s2);

   return me;
}

template<class T> double Mean(T ***a,int s1,int s2,int s3)
{
   int i,j,k; double sum,me;

   sum=0; for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++) sum+=a[i][j][k];
   me=sum/(s1*s2*s3);

   return me;
}

template<class T> double Mean(T ****a,int s1,int s2,int s3,int s4)
{
   int i,j,k,l; double sum,me;

   sum=0; for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++) for (l=0;l<s4;l++)
      sum+=a[i][j][k][l]; me=sum/(s1*s2*s3*s4);

   return me;
}


// 2.3c.2: Var/SD of a[], a[][], etc
template<class T> double Var(T *a,int s)
{
   int i; double me,sse,a1,var;

   me=Mean(a,s); sse=0; for (i=0;i<s;i++) { a1=a[i]-me; sse+=a1*a1; }
   if (s>=2) var=sse/(s-1); else var=0;

   return var;
}

template<class T> double SD(T *a,int s)
{
   int i; double me,sse,a1,sd;

   me=Mean(a,s);
   sse=0; for (i=0;i<s;i++) { a1=a[i]-me; sse+=a1*a1; }
   if (s>=2) sd=sqrt(sse/(s-1)); else sd=0;

   return sd;
}

template<class T> double SD(T **a,int s1,int s2)
{
   int s,i,j; double me,sse,a1,sd;

   me=Mean(a,s1,s2); s=s1*s2;
   sse=0; for (i=0;i<s1;i++) for (j=0;j<s2;j++) { a1=a[i][j]-me; sse+=a1*a1; }
   if (s>=2) sd=sqrt(sse/(s-1)); else sd=0;

   return sd;
}

template<class T> double SD(T ***a,int s1,int s2,int s3)
{
   int s,i,j,k; double me,sse,a1,sd;

   me=Mean(a,s1,s2,s3); s=s1*s2*s3;
   sse=0; for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++) {
      a1=a[i][j][k]-me; sse+=a1*a1; }
   if (s>=2) sd=sqrt(sse/(s-1)); else sd=0;

   return sd;
}

template<class T> double SD(T ****a,int s1,int s2,int s3,int s4)
{
   int s,i,j,k,l; double me,sse,a1,sd;

   me=Mean(a,s1,s2,s3,s4); s=s1*s2*s3*s4;
   sse=0; for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++) for (l=0;l<s4;l++) {
      a1=a[i][j][k][l]-me; sse+=a1*a1; }
   if (s>=2) sd=sqrt(sse/(s-1)); else sd=0;

   return sd;
}


// 2.3c.3: MeanPar of a[], a[][], etc
template<class T> double MeanPar(T *a, int sF, int sL)
{
   int i; double sum,me;

   sum=0; for (i=sF;i<=sL;i++) sum+=a[i]; me=sum/(sL-sF+1);

   return me;
}

template<class T> double MeanPar(T **a,int s1F,int s1L,int s2F,int s2L)
{
   int i,j; double sum,me;

   sum=0; for (i=s1F;i<=s1L;i++) for (j=s2F;j<=s2L;j++) sum+=a[i][j];
   me=sum/((s1L-s1F+1)*(s2L-s2F+1));

   return me;
}

template<class T> double MeanPar(T ***a,int s1F,int s1L,int s2F,int s2L,int s3F,int s3L)
{
   int i,j,k; double sum,me;

   sum=0; for (i=s1F;i<=s1L;i++) for (j=s2F;j<=s2L;j++) for (k=s3F;k<=s3L;k++) sum+=a[i][j][k];
   me=sum/((s1L-s1F+1)*(s2L-s2F+1)*(s3L-s3F+1));

   return me;
}


// 2.3c.4: SDPar of a[], a[][], etc
template<class T> double SDPar(T *a,int sF,int sL)
{
   int s,i; double me,sse,a1,sd;

   me=MeanPar(a,sF,sL); s=sL-sF+1;
   sse=0; for (i=sF;i<=sL;i++) { a1=a[i]-me; sse+=a1*a1; }
   if (s>=2) sd=sqrt(sse/(s-1)); else sd=0;

   return sd;
}

template<class T> double SDPar(T **a,int s1F,int s1L,int s2F,int s2L)
{
   int s,i,j; double me,sse,a1,sd;

   me=MeanPar(a,s1F,s1L,s2F,s2L); s=(s1L-s1F+1)*(s2L-s2F+1);
   sse=0; for (i=s1F;i<=s1L;i++) for (j=s2F;j<=s2L;j++) { a1=a[i][j]-me; sse+=a1*a1; }
   if (s>=2) sd=sqrt(sse/(s-1)); else sd=0;

   return sd;
}

template<class T> double SDPar(T ***a,int s1F,int s1L,int s2F,int s2L,int s3F,int s3L)
{
   int s,i,j,k; double me,sse,a1,sd;

   me=MeanPar(a,s1F,s1L,s2F,s2L,s3F,s3L); s=(s1L-s1F+1)*(s2L-s2F+1)*(s3L-s3F+1);
   sse=0; for (i=s1F;i<=s1L;i++) for (j=s2F;j<=s2L;j++) for (k=s3F;k<=s3L;k++) {
      a1=a[i][j][k]-me; sse+=a1*a1; }
   if (s>=2) sd=sqrt(sse/(s-1)); else sd=0;

   return sd;
}

// 2.3c.5: MedianPar of a[], a[][], etc
template<class T> double MedianPar(T *a,int s1F,int s1L)
{
   T *a1; int s,i1,i; double med; s=s1L-s1F+1; Array1(a1,s);

   i1=0; for (i=s1F;i<=s1L;i++) { a1[i1]=a[i]; i1++; }
   med=Median(a1,s); Drray1(a1,s); return med;
}
template<class T> double MedianPar(T **a,int s1F,int s1L,int s2F,int s2L)
{
   T *a1; int s,i1,i,j; double med; s=(s1L-s1F+1)*(s2L-s2F+1); Array1(a1,s);

   i1=0; for (i=s1F;i<=s1L;i++) for (j=s2F;j<=s2L;j++) { a1[i1]=a[i][j]; i1++; }
   med=Median(a1,s); Drray1(a1,s); return med;
}
template<class T> double MedianPar(T ***a,int s1F,int s1L,int s2F,int s2L,int s3F,int s3L)
{
   T *a1; int s,i1,i,j,k; double med; s=(s1L-s1F+1)*(s2L-s2F+1)*(s3L-s3F+1); Array1(a1,s);

   i1=0; for (i=s1F;i<=s1L;i++) for (j=s2F;j<=s2L;j++) for (k=s3F;k<=s3L;k++) {
      a1[i1]=a[i][j][k]; i1++; }
   med=Median(a1,s); Drray1(a1,s); return med;
}



//////////////////////////////////////////////////////////////////
template<class T> double Mean_Min(T *a, int s,T aMi)
{
   int i,n; T sum; double me;

   n=0; sum=0; for (i=0;i<s;i++) { if (a[i]>=aMi) { n++; sum+=a[i]; } }
   if (n>=1) me=(double)sum/n;
   return me;
}

template<class T> double SS(T *a, int s)
{
   int i; double ss=0;
   for (i=0;i<s;i++) ss += (double)a[i]*a[i];
   return ss;
}
template<class T> double SSE(T *a, int s)
{
   double sum,ss, sse;
   sum=(double)Sum(a,s); ss=(double)SS(a,s); sse=ss-sum*sum/s;
   return Max(sse,0.);
}
template<class T> double VarSk(T *a, int s) // get sample variance of skewedDist
{
   double k2,k4,kur,var;

   if      (s>=2)
   {
      k2=Cmt2(a,s); k4=Cmt4(a,s); kur=k4/(k2*k2);
      var=SSE(a,s)/(s-1.5-0.25*kur);
   }
   else if (s==1) var=0;
   else if (s<1) { cerr<<"\nERROR: Cannot compute variance for n<1.\n\n"; exit(1); }
   return var;
}

template<class T> double SD_Min(T *a,int s,T aMi)
{
   int n,i; double sum,ss,sd;

   n=0; sum=0.; ss=0.;
   for (i=0;i<s;i++)
   {
      if (a[i]>=aMi) { n++; sum += (double)a[i]; ss += (double)a[i]*a[i]; }
   }
   if (n>=2) sd=pow((ss-sum*sum/n)/(n-1),0.5); else sd=0.;

   return sd;
}

template<class T> double Mean_Wt(T *a,T *wt, int s)
{
   int i; double sum,me;

   sum=0.; for (i=0;i<s;i++) sum += a[i]*wt[i]; me=sum/s;
   return me;
}
template<class T> double SD_Wt(T *a,T *wt, int s)
{
   int i; double sum,ss,d1, sd;

   sum=0.; ss=0.;
   for (i=0;i<s;i++) { d1=a[i]*wt[i]; sum+=d1; ss+=d1*d1; }
   sd=pow((ss-sum*sum/s)/(s-1),0.5);

   return sd;
}


template<class T1,class T2> double MeanWt(T1 *a, T2 *fr, int s)
{
   return (double)SumWt(a,fr,s)/Sum(fr,s);
}




template<class T> double Mean_Min(T **a,int s1,int s2,T aMi)
{
   int i,j,n; T sum; double me;

   n=0; sum=0;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++)  { if (a[i][j]>=aMi) { n++; sum+=a[i][j]; } }
   if (n>=1) me=(double)sum/n; else me=-9.;
   return me;
}

template<class T> double SD_Min(T **a,int s1,int s2,T aMi)
{
   int n,i,j; double sum,ss,d1,sd;

   n=0; sum=0.; ss=0.;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++)
   {
      if (a[i][j]>=aMi) { n++; d1=(double)a[i][j]; sum+=d1; ss+=d1*d1; }
   }
   if (n>=2) sd=pow((ss-sum*sum/n)/(n-1),0.5); else sd=0.;

   return sd;
}

template<class T> double Mean2Fix(T **a,int s1, int j)
{
   return (double)Sum2Fix(a, s1, j)/s1;
}

template<class T> double SD2Fix(T **a,int s1, int j)
{
   int i; double sd=0, sum;
   for (i=0;i<s1;i++) sd += a[i][j]*a[i][j];
   sum = Sum2Fix(a, s1, j);
   sd = pow( (sd-sum*sum/s1)/(s1-1), 0.5 );
   return sd;
}

template<class T> double Mean1Partial2Fix(T **a, int sF, int sL, int j)
{
   return (double)Sum1Par2Fix(a, sF, sL, j)/(sL-sF+1);
}

template<class T> double SD1Partial2Fix(T **a, int sF, int sL, int j)
{
   int i; double sd=0, sum;
   for (i=sF; i<=sL; i++) sd += a[i][j]*a[i][j];
   sum = Sum1Par2Fix(a, sF, sL, j);
   sd = pow( (sd-sum*sum/(sL-sF+1))/(sL-sF), 0.5 );
   return sd;
}

template<class T> void MeanSDEtc_UpRtNoDiag(T **a, int s, int &n, double &mean,
                                             double &sd, T &min, T &max)
{
   double sum;
   n=s*(s-1)/2; sum=0; sd=0; min=a[0][1]; max=a[0][1];
   for (int i=0; i<s-1; i++) for (int j=i+1; j<s; j++)
   {
      sum += a[i][j];
      sd+=a[i][j]*a[i][j];
      if (min > a[i][j]) min=a[i][j];
      if (max < a[i][j]) max=a[i][j];
   }
   mean=sum/n; sd=pow((sd-sum*sum/n)/(n-1), 0.5);
}

template<class T> void MeanSDEtc_Diag1Up(T **a, int s, int &n, double &mean,
                                          double &sd, T &min, T &max)
{
   double sum;
   n=s-1; sum=0; sd=0; min=a[0][1]; max=a[0][1];
   for (int i=0; i<s-1; i++)
   {
      sum += a[i][i+1]; sd+=a[i][i+1]*a[i][i+1];
      if (min > a[i][i+1]) min=a[i][i+1];
      if (max < a[i][i+1]) max=a[i][i+1];
   }
   mean=sum/n; sd=pow((sd-sum*sum/n)/(n-1), 0.5);
}

template<class T> double MeanDiag1Up(T **a, int s)
{
   double sum=0;
   for (int i=0; i<s-1; i++) sum += a[i][i+1];
   return (double)sum/(s-1);
}

template<class T> double SDDiag1Up(T **a, int s)
{
   double sum=0, sd=0;
   for (int i=0; i<s-1; i++) { sum+=a[i][i+1]; sd+=a[i][i+1]*a[i][i+1]; }
   sd=sd-sum*sum/(s-1); sd=pow(sd/(s-2), 0.5);
   return sd;
}




template<class T> double Cov(T *a,T *b, int s)
{
   int i; double cov=0;
   for (i=0; i<s; i++) cov += a[i]*b[i];
   cov = cov-Sum(a,s)*Sum(b,s)/s;
   return cov;
}

template<class T> double Corr(T *a,T *b, int s)
{
   double cov,se1,se2,cc;
   cov=Cov(a,b,s); se1=SSE(a,s); se2=SSE(b,s);
   if (se1==0. || se2==0.) cc=0.; else cc=cov/pow(se1*se2,0.5);
   return cc;
}

// get skewness sk=m3/sg^3
template<class T> double Skew(T *a, int s)
{
   int i; double sum,ss,me,v,sk;

   if (s>=2)
   {
      sum=0; ss=0; for (i=0;i<s;i++) { sum+=(double)a[i]; ss+=(double)a[i]*a[i]; }
      me=(double)sum/s; v=(ss-sum*sum/s)/s;
      sk=0.; for (i=0;i<s;i++) sk += pow((double)a[i]-me,3);
      sk=(sk/s)/pow(v,1.5); //sk*=sqrt(s*(s-1))/(s-2);
   }
   else sk=0.;

   return sk;
}
template<class T> double Skew_Min(T *a,int s,T aMi)
{
   int n,i; double sum,ss,me,sd,sk;

   n=0; sum=0.; ss=0.;
   for (i=0;i<s;i++)
   {
      if (a[i]>=aMi) { n++; sum+=(double)a[i]; ss+=(double)a[i]*a[i]; }
   }
   if (n>=2)
   {
      me=(double)sum/n; sd=pow((ss-sum*sum/n)/(n-1),0.5);
      sk=0.; for (i=0;i<s;i++) { if (a[i]>=aMi) { sk += pow(((double)a[i]-me)/sd,3); } }
      sk=(sk/n)*(sqrt(n*(n-1))/(n-2));
   }
   else sk=0.;
   return sk;
}

// get kurtosis ku=m4/sg^4
template<class T> double Kurt(T *a, int s)
{
   int i; double sum,ss,me,v,ku;

   if (s>=2)
   {
      sum=0; ss=0; for (i=0;i<s;i++) { sum+=(double)a[i]; ss+=(double)a[i]*a[i]; }
      me=(double)sum/s; v=(ss-sum*sum/s)/s;
      ku=0.; for (i=0;i<s;i++) ku += pow((double)a[i]-me,4);
      ku=(ku/s)/(v*v);
   }
   else ku=0.;

   return ku;
}


// get k1=Cumulant_1 [= nc.m1=non-central-moment1=Mean(a,n)] = sum_a/n
template<class T> double Cmt1(T *a,int n)
{
   int i; double k1;

   k1=0; for (i=0;i<n;i++) k1+=a[i]; k1=k1/n;
   return k1;
}

// get k2=Cmt2 [= m2=centralMmt2=Var(a,n)] = sum_a^2/(n-1)
template<class T> double Cmt2(T *a,int n)
{
   int i; double k1,k2,d1;

   k1=Cmt1(a,n); k2=0; for (i=0;i<n;i++) { d1=a[i]-k1; k2+=d1*d1; }
   k2=k2/(double)(n-1);
   return k2;
}

// get k3=Cmt3 [~ m3=centralMmt3] = sum_a^3*n/((n-1)*(n-2))
template<class T> double Cmt3(T *a,int n)
{
   int i; double k1,m3,k3,d1;

   k1=Cmt1(a,n); m3=0; for (i=0;i<n;i++) { d1=a[i]-k1; m3+=d1*d1*d1; }
   m3/=n; k3=m3*((double)n/(n-1))*((double)n/(n-2));
   return k3;
}

// get k4=Cmt4 [~ m4-3*m2*m2] = (n/(n-2))*(n/(n-3))*(m4*(n+1)/(n-1)-3*m2*m2)
template<class T> double Cmt4(T *a,int n)
{
   int i; double k1,m2,m4,d1,k4; // m2/m4=2nd/4th_centralMoments

   k1=Cmt1(a,n); m2=0; m4=0; for (i=0;i<n;i++) { d1=a[i]-k1; m2+=d1*d1; m4+=pow(d1,4.); }
   m2/=n; m4/=n; k4=((double)n/(n-2))*((double)n/(n-3))*(m4*(double)(n+1)/(n-1)-3*m2*m2);
   return k4;
}

// get k1-k4=cumulants 1-4 and sd(k1),...,sd(k4)
void getkSdk1234_dfNc(double df,double nc, double *k,double *sdk)
{
   int i; double d1,k1,k2,k3,k4,k5,k6,k8,*da; Array1(da,8);

   for (i=0;i<8;i++) { if (i==0) d1=1; else d1=d1*2*i; da[i]=d1*(df+(i+1)*nc); }
   for (i=0;i<4;i++) k[i]=da[i];
   k1=da[0]; k2=da[1]; k3=da[2]; k4=da[3]; k5=da[4]; k6=da[5]; k8=da[7];
   sdk[0]=sqrt(k2);
   sdk[1]=sqrt(k4+2*k2*k2);
   sdk[2]=sqrt(k6+9*k2*k4+9*k3*k3+6*k2*k2*k2);
   sdk[3]=sqrt(k8+16*k2*k6+48*k3*k5+34*k4*k4 +72*k2*k2*k4+144*k2*k3*k3+24*pow(k2,4.));

   Drray1(da,8);
}



// [2.3d] MIN, MAX, MEDIAN AND PERCENTILES

// minimums, maximums
template<class T> T Min(T a, T b) { return (a<b) ? a : b; }
template<class T> T Max(T a, T b) { return (a>b) ? a : b; }

template<class T> T Min(T a, T b, T c) { return Min(Min(a, b), c); }
template<class T> T Max(T a, T b, T c) { return Max(Max(a, b), c); }

template<class T> T Min(T a, T b, T c, T d) { return Min(Min(a,b),Min(c,d)); }
template<class T> T Max(T a, T b, T c, T d) { return Max(Max(a,b),Max(c,d)); }

template<class T> T Min(T a,T b,T c,T d,T e) { return Min(Min(a,b,c),Min(d,e)); }
template<class T> T Max(T a,T b,T c,T d,T e) { return Max(Max(a,b,c),Max(d,e)); }

template<class T> T Min(T a,T b,T c,T d,T e,T f) { return Min(Min(a,b,c),Min(d,e,f)); }
template<class T> T Max(T a,T b,T c,T d,T e,T f) { return Max(Max(a,b,c),Max(d,e,f)); }

template<class T> T Min(T *a, int s)
{
   int i; T min = a[0];
   for (i=0;i<s;i++) { if (a[i]<min) min=a[i]; }
   return min;
}
template<class T> T Max(T *a, int s)
{
   int i; T max = a[0];
   for (i=0;i<s;i++) { if (a[i]>max) max=a[i]; }
   return max;
}

template<class T> T MinPartial(T *a, int sMin, int sMax)
{
   int i; T min=a[sMin];
   for (i=sMin;i<=sMax;i++) { if (a[i]<min) min=a[i]; }
   return min;
}
template<class T> T MaxPartial(T *a, int sMin, int sMax)
{
   int i; T max=a[sMin];
   for (i=sMin;i<=sMax;i++) { if (a[i]>max) max=a[i]; }
   return max;
}

// get nth minimal value in array xs[n]
template<class T> T MinNth(T *xs,int n, int nth)
{
   int i,j,j1; T d1,xNth,*xt; Array1(xt,nth); // xt[]=nMinXs[]TailValues

   d1=Max(xs,n); for (j=0;j<nth;j++) xt[j]=d1;
   for (i=0;i<n;i++)
   {
      d1=xs[i];
      if (d1<=xt[0]) // xt[0]=nth min
      {
         for (j=0;j<nth;j++)
         {
            if ((j<=nth-2 && d1<=xt[j] && d1>xt[j+1]) || (j==nth-1 && d1<=xt[nth-1]))
            {
               for (j1=0;j1<j;j1++) xt[j1]=xt[j1+1]; xt[j]=d1; break;
            }
         }
      }
   }
   xNth=xt[0]; Drray1(xt,nth);

   return xNth;
}

// get nth maximal value in array xs[n]
template<class T> T MaxNth(T *xs,int n, int nth)
{
   int i,j,j1; T d1,xNth,*xt; Array1(xt,nth); // xt[]=nMaxXs[]TailValues

   d1=Min(xs,n); for (j=0;j<nth;j++) xt[j]=d1;
   for (i=0;i<n;i++)
   {
      d1=xs[i];
      if (d1>=xt[0]) // xt[0]=nth max
      {
         for (j=0;j<nth;j++)
         {
            if ((j<=nth-2 && d1>=xt[j] && d1<xt[j+1]) || (j==nth-1 && d1>=xt[nth-1]))
            {
               for (j1=0;j1<j;j1++) xt[j1]=xt[j1+1]; xt[j]=d1; break;
            }
         }
      }
   }
   xNth=xt[0]; Drray1(xt,nth);

   return xNth;
}


template<class T> T Min(T **a,int s1, int s2)
{
   int i, j;
   T min = a[0][0];
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) { if (min > a[i][j]) min=a[i][j]; }
   return min;
}
template<class T> T Max(T **a,int s1, int s2)
{
   int i, j;
   T max = a[0][0];
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) { if (max < a[i][j]) max=a[i][j]; }
   return max;
}

template<class T> T Min(T **a,int s1,int *s2)
{
   int i, j;
   T min = a[0][0];
   for (i=0;i<s1;i++) for (j=0; j<s2[i]; j++) { if (min > a[i][j]) min=a[i][j]; }
   return min;
}
template<class T> T Max(T **a,int s1,int *s2)
{
   int i, j;
   T max = a[0][0];
   for (i=0;i<s1;i++) for (j=0; j<s2[i]; j++) { if (max < a[i][j]) max=a[i][j]; }
   return max;
}

template<class T> T Min2Fix(T **a,int s1, int j)
{
   T min = a[0][j];
   for (int i=0;i<s1;i++) { if (min > a[i][j]) min=a[i][j]; }
   return min;
}
template<class T> T Max2Fix(T **a,int s1, int j)
{
   T max = a[0][j];
   for (int i=0;i<s1;i++) { if (max < a[i][j]) max=a[i][j]; }
   return max;
}

template<class T> T Max1Sum2(T **a,int s1, int s2)
{
   T max; int i;

   max=Sum(a[0],s2); for (i=1;i<s1;i++) { max=Max(max,Sum(a[i],s2)); }
   return max;
}


template<class T> T Min(T ***a,int s1,int s2,int s3)
{
   int i, j, k;
   T min = a[0][0][0];
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++)
   {
      if (min > a[i][j][k]) min=a[i][j][k];
   }
   return min;
}
template<class T> T Max(T ***a,int s1,int s2,int s3)
{
   int i, j, k;
   T max = a[0][0][0];
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++)
   {
      if (max < a[i][j][k]) max=a[i][j][k];
   }
   return max;
}

template<class T> T Min(T ****a,int s1,int s2,int s3,int s4)
{
   int i, j, k, l;
   T min = a[0][0][0][0];
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++) for (l=0;l<s4;l++)
   {
      if (min > a[i][j][k][l]) min=a[i][j][k][l];
   }
   return min;
}
template<class T> T Max(T ****a,int s1,int s2,int s3,int s4)
{
   int i, j, k, l;
   T max = a[0][0][0][0];
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++) for (l=0;l<s4;l++)
   {
      if (max < a[i][j][k][l]) max=a[i][j][k][l];
   }
   return max;
}

template<class T> T Min(T ****a,int s1,int *s2,int s3,int s4)
{
   int i, j, k, l;
   T min = a[0][0][0][0];
   for (i=0;i<s1;i++) for (j=0;j<s2[i];j++) for (k=0;k<s3;k++) for (l=0;l<s4;l++)
   {
      if (min > a[i][j][k][l]) min=a[i][j][k][l];
   }
   return min;
}
template<class T> T Max(T ****a,int s1,int *s2,int s3,int s4)
{
   int i, j, k, l;
   T max = a[0][0][0][0];
   for (i=0;i<s1;i++) for (j=0;j<s2[i];j++) for (k=0;k<s3;k++) for (l=0;l<s4;l++)
   {
      if (max < a[i][j][k][l]) max=a[i][j][k][l];
   }
   return max;
}

// medians
template<class T> double Median(T *a, int s)
{
   int i; double median;
   T *a_jk; Array1(a_jk, s);

   for (i=0; i<s; i++) a_jk[i]=a[i];
   sortArrayAsc(a_jk, s);
   median=0.5*(a_jk[(s-1)/2]+a_jk[s/2]);

   Drray1(a_jk, s);
   return median;
}

template<class T> double Q1(T *a, int s)
{
   int i; double q1;
   T *a_jk; Array1(a_jk, s);

   for (i=0; i<s; i++) a_jk[i]=a[i];
   sortArrayAsc(a_jk, s);
   q1=0.5*(a_jk[(s-1)/4]+a_jk[s/4]);

   Drray1(a_jk, s);
   return q1;
}

template<class T> double Q3(T *a, int s)
{
   int i; double q3;
   T *a_jk; Array1(a_jk, s);

   for (i=0; i<s; i++) a_jk[i]=a[i];
   sortArrayAsc(a_jk, s);
   q3=0.5*(a_jk[(3*s-1)/4]+a_jk[3*s/4]);

   Drray1(a_jk, s);
   return q3;
}



//--------------------------------------------------------------------
// [2.4] CALCULATIONS AND UTILITIES (templates)
//--------------------------------------------------------------------

// [2.4a] DIGIT/DICIMAL

// get number of digits of an integer or double
template<class T> int digit(T v)
{
   int dgt;
   if (v>=1 || v<=-1) dgt=(int)log10(fabs((double)v))+1;
   else dgt=1;
   return dgt;
}
template<class T> int digit(T *v,int s1)
{
   int dgt,i;
   dgt=0; for (i=0;i<s1;i++) dgt=Max(dgt, digit(v[i]));
   return dgt;
}
template<class T> int digit(T **v,int s1,int s2)
{
   int dgt,i,j;
   dgt=0; for (i=0;i<s1;i++) for (j=0;j<s2;j++) dgt=Max(dgt, digit(v[i][j]));
   return dgt;
}
template<class T> int digit(T **v,int s1,int *s2)
{
   int dgt,i,j;
   dgt=0; for (i=0;i<s1;i++) for (j=0;j<s2[i];j++) dgt=Max(dgt, digit(v[i][j]));
   return dgt;
}


// get number of decimals of any double
int decimal(double d)
{
   char c[32],c1[32]; int sz,iExp,iNeg,iPos,iDot,i,dec;
   // reset to the default C++ print format (seems unnecessary)

//   sprintf_s(c,32,"%g",d); // "sprintf_s" is security enhanced "sprintf", but invalid in Linux
//   sz=(int)strlen(c);
   sz=sprintf(c,"%g",(double)d); // write formatted data to string, valid in both Windows & Linux

   iExp=0; for (i=0;i<sz;i++) { if (c[i]=='e') { iExp=i+1; break; } }
   iNeg=0; for (i=0;i<sz;i++) { if (c[i]=='-') { iNeg=i+1; break; } }
   iPos=0; for (i=0;i<sz;i++) { if (c[i]=='+') { iPos=i+1; break; } }
   iDot=0; for (i=0;i<sz;i++) { if (c[i]=='.') { iDot=i+1; break; } }
   if (iExp>=1 && iNeg==iExp+1)
   {
      strncpy(c1,c+iNeg,sz-iNeg); c1[sz-iNeg]='\0';
      dec=atoi(c1);
   }
   else if (iExp>=1 && iPos==iExp+1)
   {
      strncpy(c1,c+iPos,sz-iPos); c1[sz-iPos]='\0';
      dec=Max((iExp-3)-(atoi(c1)-1),0);
   }
   else dec=(iDot>=1?(int)strlen(c)-iDot:0);
   return dec;
}

// get maximum decimal of double array
int decimal(double *d,int s1)
{
   int i,dec; dec=0; for (i=0;i<s1;i++) dec=Max(dec, decimal(d[i]));
   return dec;
}
int decimal(double **d,int s1,int s2)
{
   int i,j,dec;
   dec=0; for (i=0;i<s1;i++) for (j=0;j<s2;j++) dec=Max(dec, decimal(d[i][j]));
   return dec;
}
int decimal(double **d,int s1,int *s2)
{
   int i,j,dec;
   dec=0; for (i=0;i<s1;i++) for (j=0;j<s2[i];j++) dec=Max(dec, decimal(d[i][j]));
   return dec;
}


// get width of any double value
int width(double d)
{
   int dec,w;

   dec=decimal(d);
   w=digit(d)+(dec>=1?dec+1:0);
   return w;
}


// get number of decimals allowing nnz=(n- non-zeros)
// ex: (dec|nnz=1)=3 for 0.00345, (dec|nnz=2)=4 for 0.00345
template<class T> int decnnz(T d,int nnz)
{
   int i,nZ,nNZ; double d1,d2;

   if (d<0) d=-d; d1=((double)d-(int)d)+1e-12;
   nZ=0; for (i=1;i<=8;i++) { if (d1<pow(0.1,(double)i)) nZ=i; else break; } if (nZ==8) nZ=0;
   nNZ=0; for (i=0;i<nnz;i++) { d2=d1*pow(10.,nZ+i); if ((d2-(int)d2)>0.001) nNZ++; }
   return nZ+nNZ;
}

// get maximum decimal of double array allowing nnz=(n- non-zeros)
template<class T> int decnnz(T *d,int s1,int nnz)
{
   int dec,i; dec=0; for (i=0;i<s1;i++) { dec=Max(dec, decnnz(d[i],nnz)); }
   return dec;
}
template<class T> int decnnz(T **d,int s1,int s2,int nnz)
{
   int dec,i,j; dec=0; for (i=0;i<s1;i++) for (j=0;j<s2;j++) { dec=Max(dec, decnnz(d[i][j],nnz)); }
   return dec;
}

// get floor, ceil and round of a value (integer or double)
template<class T> T Floor(T a,double unit)
{
   T aRtn; int n;
   if (unit<1.)
   {
      n=(int)pow(10.,(int)(-log10(unit)+0.000000001));
      aRtn=(T)((double)((int)(a*n))/n);
   }
   else if (unit>=1.)
   {
      n=(int)pow(10.,(int)(log10(unit)+0.000000001));
      aRtn=(T)((double)((int)(a/n))*n);
   }
   return aRtn;
}
template<class T> int Floor(T a)
{
   return (int)Floor(a,1.);
}

template<class T> T Ceil(T a,double unit)
{
   T aRtn; int n;
   if (unit<1.)
   {
      n=(int)pow(10.,(int)(-log10(unit)+0.000000001));
      aRtn=(T)((double)((int)((double)a*n+0.99999999))/n);
   }
   else if (unit>=1.)
   {
      n=(int)pow(10.,(int)(log10(unit)+0.000000001));
      aRtn=(T)((double)((int)((double)a/n+0.99999999))*n);
   }
   return aRtn;
}
template<class T> int Ceil(T a)
{
   return (int)Ceil(a,1.);
}

template<class T> T Round(T a,double unit)
{
   T aRtn; int n;
   if (unit<1.)
   {
      n=(int)pow(10.,(int)(-log10(unit)+0.0000000000001));
      if ((double)a>=0.) aRtn=(T)((double)((int)((double)a*n+0.5))/n);
      else               aRtn=(T)((double)((int)((double)a*n-0.5))/n);
   }
   else if (unit>=1.)
   {
      n=(int)pow(10.,(int)(log10(unit)+0.0000000000001));
      if ((double)a>=0.) aRtn=(T)((double)((int)((double)a/n+0.5))*n);
      else               aRtn=(T)((double)((int)((double)a/n-0.5))*n);
   }
   return aRtn;
}
template<class T> int Round(T a)
{
   return (int)Round(a,1.);
}



// [2.4b] CHECK RANGE (INTEGER OR DOUBLE)

// check a value (integer or double)
template<class T> void chkVal(T v, char *vNm, T val)
{
   if (v!=val)
   {
      cerr<<"\nERROR: Variable "<<vNm<<" = "<<v<<" is not equal to "<<val<<".\n\n"; exit(1);
   }
}

// check 2 values (integer or double)
template<class T> void chkV2(T v, char *vNm, T v1,T v2)
{
   if (v!=v1 && v!=v2)
   {
      cerr<<"\nERROR: Variable "<<vNm<<" = "<<v<<" is not equal to "<<v1<<" or "<<v2<<".\n\n"; exit(1);
   }
}

// check range of a value (integer or double)
template<class T> void chkRng(T v, char *vNm, T vmin, T vmax)
{
   if (v<vmin||v>vmax)
   {
      cerr<<"\nERROR: Variable "<<vNm<<" = "<<v<<" is out of range ("<<vmin<<", "<<vmax<<").\n\n"; exit(1);
   }
}

// check range of an array
template<class T> void chkRng(T *v,int vsz, char *vNm, T vmin, T vmax)
{
   for (int i=0;i<vsz;i++)
   {
      if (v[i]<vmin||v[i]>vmax)
      {
         cerr<<"\nERROR: Variable " << vNm << "[" << i << "] = " << v[i] << " is out of range ("
            <<vmin<<", "<<vmax<<").\n\n"; exit(1);
      }
   }
}

// check min of a value
template<class T> void chkMin(T v, char *vNm, T vmin)
{
   if (v<vmin)
   {
      cerr<<"\nERROR: Variable "<<vNm<<" = "<<v<<" is less than "<<vmin<<".\n\n"; exit(1);
   }
}

// check min of an array
template<class T> void chkMin(T *v,int vsz, char *vNm, T vmin)
{
   for (int i=0;i<vsz;i++)
   {
      if (v[i]<vmin)
      {
         cerr<<"\nERROR: Variable " << vNm << "[" << i << "] = " << v[i] << " is less than "
            <<vmin<<".\n\n"; exit(1);
      }
   }
}

// check max of a value
template<class T> void chkMax(T v, char *vNm, T vmax)
{
   if (v>vmax)
   {
      cerr<<"\nERROR: Variable "<<vNm<<" = "<<v<<" is greater than "<<vmax<<".\n\n"; exit(1);
   }
}



// [2.4c] PRINT ARRAY

// simple print arrays
template<class T> void prtAr1(ostream &out,T *a,int s1)
{
   int i;
   for (i=0;i<s1;i++) out<<(i==0?"":" ")<<a[i]; out<<endl;
}

// get w+dec+out co array & nnz=nNonZero decimals
template<class T> void getWSetDec_Nnz(ostream &out,T *a,int s1, int nnz, int &w)
{
   T ma,mi; int dec;
   ma=Max(a,s1); mi=Min(a,s1); w=Max((ma>0?digit(ma):0),(mi<0?digit(mi)+1:0)); dec=decnnz(a,s1,nnz);
   if (dec>=1) { w=w+dec+1; out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(dec); }
}
template<class T> void getWSetDec_Nnz(ostream &out,T **a,int s1,int s2, int nnz,int &w)
{
   T ma,mi; int dec;
   ma=Max(a,s1,s2); mi=Min(a,s1,s2); w=Max((ma>0?digit(ma):0),(mi<0?digit(mi)+1:0)); 
   dec=decnnz(a,s1,s2,nnz);
   if (dec>=1) { w=w+dec+1; out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(dec); }
}

// print 1-dim array w/ or w/o return
template<class T> void prtArray1(ostream &out,T *a,int s1,char rt[]="\n")
{
   int w,i; getWSetDec_Nnz(out,a,s1,2, w);
   for (i=0;i<s1;i++) out<<(i==0?"":" ")<<a[i]<<rt;
}

// print 1/2/3-dim char array
void prtArray1(ostream &out,char **a,int s1,char rt[]="\n")
{
   int i; for (i=0;i<s1;i++) out<<(i==0?"":" ")<<a[i]<<rt;
}
void prtArray2(ostream &out,char ***a,int s1,int s2)
{
   int i,j;
   for (i=0;i<s1;i++) { for (j=0;j<s2;j++) out<<(j==0?"":" ")<<a[i][j]; out<<endl; }
}
void prtArray3(ostream &out,char ****a,int s1,int s2,int s3)
{
   int i,j,k;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++)
   {
      for (k=0;k<s3;k++) out<<(k==0?"":" ")<<a[i][j][k]; out<<endl;
   }
}
void prtArray3(ostream &out,char ****a,int s1,int *s2,int s3)
{
   int i,j,k;
   for (i=0;i<s1;i++) for (j=0;j<s2[i];j++)
   {
      for (k=0;k<s3;k++) out<<(k==0?"":" ")<<a[i][j][k]; out<<endl;
   }
}


// prt 1-d array parIdxFm sF-sL
template<class T> void prtArray1_P(ostream &out,T *a,int sF, int sL,char rt[]="\n")
{
   int i; for (i=sF;i<=sL;i++) out<<(i==sF?"":" ")<<a[i]; out<<rt;
}


// print 1-dim array by row w/o return
template<class T> void prtArray1BRNR(ostream &out,T *a,int s1, int nPR)
{
   int w,i; getWSetDec_Nnz(out,a,s1, 2,w);

   for (i=0;i<s1;i++)
   {
      out<<(i%nPR==0?"":" ")<<setw(w)<<a[i]; if (i%nPR==(nPR-1) && i<s1-1) out<<endl;
   }
}
template<class T> void prtArray1BR(ostream &out,T *a,int s1, int nPR)
{
   prtArray1BRNR(out,a,s1,nPR); out<<endl;
}

template<class T> void prtArray1FWNR(ostream &out,T *a,int s1, int w)
{
   for (int i=0;i<s1;i++) out<<setw(w)<<a[i];
}
template<class T> void prtArray1FW(ostream &out,T *a,int s1, int w)
{
   prtArray1FWNR(out,a,s1,w); out<<endl;
}

template<class T> void prtArray1FDNR(ostream &out,T *a,int s1, int dec)
{
   int i;
   out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(dec);
   for (i=0;i<s1;i++) out<<(i==0?"":" ")<<(double)a[i];
}
template<class T> void prtArray1FD(ostream &out,T *a,int s1, int dec)
{
   prtArray1FDNR(out,a,s1,dec); out<<endl;
}

template<class T> void prtArray1FWD(ostream &out,T *a,int s1, int w,int dec)
{
   int i;
   out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(dec);
   for (i=0;i<s1;i++) out<<setw(w)<<(double)a[i]; out<<endl;
}

template<class T> void prtArray1FWBRNR(ostream &out,T *a,int s1, int w,int nPR)
{
   int i;
   for (i=0;i<s1;i++)
   {
      out<<((i%nPR==0 && i>=1)?"\n":"")<<setw(w)<<a[i];
   }
}
template<class T> void prtArray1FWBR(ostream &out,T *a,int s1, int w,int nPR)
{
   prtArray1FWBRNR(out,a, s1,w,nPR); out<<endl;
}

template<class T> void prtArray1FWDBRNR(ostream &out,T *a,int s1, int w,int dec,int nPR)
{
   int i,w0,w1,i1; w0=digit(Max(a,s1)); w1=w-(dec>=1?dec+1:0);
   out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(dec);
   for (i=0;i<s1;i++)
   {
      i1=digit(a[i]);
      out<<((i%nPR==0 && i>=1)?"\n":"")<<setw(i%nPR==0?w0-i1:w1-i1)<<""<<a[i];
   }
}
template<class T> void prtArray1FWDBR(ostream &out,T *a,int s1, int w,int dec,int nPR)
{
   prtArray1FWDBRNR(out,a,s1, w,dec,nPR); out<<endl;
}

// print 1-dim array w/ delimiter "c"
template<class T> void prtArray1DlmNR(ostream &out,T *a,int s1, char *c)
{
   for (int i=0;i<s1;i++) out<<(i==0?"":c)<<a[i];
}
template<class T> void prtArray1Dlm(ostream &out,T *a,int s1, char *c)
{
   prtArray1DlmNR(out,a, s1,c); out<<endl;
}
template<class T> void prtArray1DlmBRNR(ostream &out,T *a,int s1, char *c,int nPR)
{
   for (int i=0;i<s1;i++)
      out<<((i%nPR==0 && i>=1)?"\n":"")<<(i==0?"":c)<<a[i];
}
template<class T> void prtArray1DlmBR(ostream &out,T *a,int s1, char *c,int nPR)
{
   prtArray1DlmBRNR(out,a, s1,c,nPR); out<<endl;
}

// print 1-d array by row to file
template<class T> void prtArray1BRTF(T *a,int s1,int nPR, char *fn)
{
   ofstream out; openOF(out,fn);
   prtArray1BR(out,a,s1,10);
   out.close();
}


// print 2-dim array w/ nnz=(n- non-zeros)
template<class T> void prtArray2(ostream &out,T **a,int s1,int s2, int nnz)
{
   int w,i,j; getWSetDec_Nnz(out,a,s1,s2, nnz,w);
   for (i=0;i<s1;i++) { for (j=0;j<s2;j++) out<<(j==0?"":" ")<<setw(w)<<a[i][j]; out<<endl; }
}
template<class T> void prtArray2(ostream &out,T **a,int s1,int s2) { prtArray2(out,a,s1,s2,2); }

template<class T> void prtMtx_UpRtNoDiag(ostream &out,T **a,int s1)
{
   int w,i,j; getWSetDec_Nnz(out,a,s1,s1,2, w);

   for (i=0;i<s1-1;i++)
   {
      for (j=1;j<s1;j++) { out<<(j==1?"":" "); if (j>i) out<<setw(w)<<a[i][j]; else out<<setw(w)<<""; }
      out<<endl;
   }
}

template<class T> void prtArray2FWD(ostream &out,T **a,int s1,int s2, int w,int dec)
{
   int i,j;
   out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(dec);
   for (i=0;i<s1;i++) { for (j=0;j<s2;j++) out<<setw(w)<<a[i][j]; out<<endl; }
}

template<class T> void prtArray2_2P(ostream &out,T **a,int s1,int s2Min,int s2Max)
{
   int w,i,j; getWSetDec_Nnz(out,a,s1,s2Max,2, w);

   for (i=0;i<s1;i++)
   {
      for (j=s2Min;j<=s2Max;j++) out<<(j==s2Min?"":" ")<<setw(w)<<a[i][j]; out<<endl;
   }
}

// print index-by-row & 2-d array in csv format
template<class T> void prtArray2_idCsv(char *fn,int *idx,T **a,int s1,int s2, int nnz)
{
   int w,i,j;

   ofstream out; openOF(out,fn);
   getWSetDec_Nnz(out,a,s1,s2, nnz,w);
   for (i=0;i<s1;i++)
   {
      out<<idx[i]; for (j=0;j<s2;j++) out<<","<<setw(w)<<a[i][j]; out<<endl;
   }
}
template<class T> void prtArray2_idCsv(char *fn,int *idx,T **a,int s1,int s2)
{
   prtArray2_idCsv(fn,idx,a,s1,s2,2);
}


template<class T> void prtArray3(ostream &out,T ***a,int s1,int s2,int s3)
{
   int iDgt=digit(s1),jDgt=digit(s2), i,j,k;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++)
   {
      out<<setw(iDgt)<<i<<setw(jDgt+1)<<j<<":";
      for (k=0;k<s3;k++) out<<" "<<a[i][j][k]; out<<endl;
   }
}

template<class T> void prtArray3_1P(ostream &out,T ***a,int s1Min,int s1Max,int s2,int s3)
{
   int iDgt=digit(s1Max),jDgt=digit(s2),aDgt=digit(Max(a,s1Max+1,s2,s3)), i,j,k;
   for (i=s1Min;i<=s1Max;i++) for (j=0;j<s2;j++)
   {
      out<<setw(iDgt)<<i<<setw(jDgt+1)<<j<<":";
      for (k=0;k<s3;k++) out<<setw(aDgt+1)<<a[i][j][k]; out<<endl;
   }
}

template<class T> void prtArray3_3P(ostream &out,T ***a,int s1,int s2,int s3Min,int s3Max)
{
   int iDgt=digit(s1),jDgt=digit(s2), i,j,k;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++)
   {
      out<<setw(iDgt)<<i<<setw(jDgt+1)<<j<<":";
      for (k=s3Min;k<=s3Max;k++) out<<" "<<a[i][j][k]; out<<endl;
   }
}

template<class T> void prtArray4(ostream &out,T ****a,int s1,int s2,int s3,int s4)
{
   int iDgt=digit(s1),jDgt=digit(s2),kDgt=digit(s3), i,j,k,l;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++)
   {
      out<<setw(iDgt)<<i<<setw(jDgt+1)<<j<<setw(kDgt+1)<<k<<": ";
      for (l=0;l<s4;l++) out<<" "<<a[i][j][k][l]; out<<endl;
   }
}

template<class T> void prtArray4(ostream &out,T ****a,int s1,int *s2,int s3,int s4)
{
   int iDgt=digit(s1), jDgt=digit(Max(s2,s1)), kDgt=digit(s3), i,j,k,l;
   for (i=0;i<s1;i++) for (j=0; j<s2[i]; j++) for (k=0;k<s3;k++)
   {
      out<<setw(iDgt)<<i<<setw(jDgt+1)<<j<<setw(kDgt+1)<<k<<": ";
      for (l=0;l<s4;l++) out<<" "<<a[i][j][k][l]; out<<endl;
   }
}

template<class T> void prtArray5(ostream &out,T *****a,int s1,int s2,int s3,int s4,int s5)
{
   int iDgt=digit(s1),jDgt=digit(s2),kDgt=digit(s3),lDgt=digit(s4), i,j,k,l,m;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++) for (k=0;k<s3;k++) for (l=0;l<s4;l++)
   {
      out<<setw(iDgt)<<i<<setw(jDgt+1)<<j<<setw(kDgt+1)<<k<<setw(lDgt+1)<<l<<": ";
      for (m=0;m<s5;m++) out<<" "<<a[i][j][k][l][m]; out<<endl;
   }
}


// [2.4d] PRINT OR GET STATS nGroup+Freq+MeanSd

// get number of groups in 1-d array
template<class T> int nGroup(T *a, int s)
{
   int i, g, occurred, nGp;
   T *gVal; Array1(gVal, s);

   nGp=0;
   for (i=0; i<s; i++)
   {
      if (i==0) { gVal[nGp]=a[i]; nGp++; }
      else
      {
         occurred=0;
         for (g=0; g<nGp; g++)
         {
            if (abs(a[i]-gVal[g])<1e-9) { occurred=1; break; }
         }
         if (occurred==0) { gVal[nGp]=a[i]; nGp++; }
      }
   }

   Drray1(gVal, s);
   return nGp;
}
// get number of groups wi-a-rg in 1-d array
template<class T> int nGroup(T *a,int s, T mi,T ma)
{
   int i,s2,nGp; T *b; Array1(b,s);

   s2=0; for (i=0;i<s;i++) { if (a[i]>=mi && a[i]<=ma) { b[s2]=a[i]; s2++; } }
   nGp=nGroup(b,s2); return nGp;
}


// get number of groups in 2-d array
template<class T> int nGroup(T **a,int s1, int s2)
{
   int i, j, g, occurred, nGp;
   T *gVal; Array1(gVal, s1*s2);

   nGp=0;
   for (i=0;i<s1;i++) for (j=0;j<s2;j++)
   {
      if (i==0 && j==0) { gVal[nGp]=a[i][j]; nGp++; }
      else
      {
         occurred=0;
         for (g=0; g<nGp; g++) if (a[i][j]==gVal[g]) { occurred=1; break; }
         if (occurred==0) { gVal[nGp]=a[i][j]; nGp++; }
      }
   }

   Drray1(gVal, s1*s2);
   return nGp;
}


// get frequency in 1-d array (input=a[], s, output=gVal[], gCnt[])
template<class T> void Freq(T *a, int s,T *gVal, int *gCnt)
{
   int g, i, j, occurred, gSize; gSize=nGroup(a, s);
   T aMin, aMax;
   T *gVal2;   Array1(gVal2, gSize);
   int *gCnt2; Array1(gCnt2, gSize);

   aMin=Min(a, s); aMax=Max(a, s);
   for (g=0; g<gSize; g++)
   {
      if (g==0) gVal[g]=aMin;
      else if (g>=1)
      {
         gVal[g]=aMax;
         for (i=0; i<s; i++) if (a[i]>gVal[g-1] && gVal[g]>a[i]) gVal[g]=a[i];
      }
   }
   for (g=0; g<gSize; g++)
   {
      gCnt[g]=0;
      for (i=0;i<s;i++) if (a[i]==gVal[g]||fabs((double)(a[i]-gVal[g]))<0.0001) ++gCnt[g];
   }

   // sort freq by gCnt:
   for (g=0; g<gSize; g++)
   {
      gVal2[g]=aMin-1; gCnt2[g]=0;
      for (i=0; i<gSize; i++)
      {
         if (g==0 && (gCnt[i]>gCnt2[g]||(gCnt[i]==gCnt2[g] && gVal[i]<gVal2[g])))
         { gCnt2[g]=gCnt[i]; gVal2[g]=gVal[i]; }
         else if (g>=1 && gCnt[i]>=gCnt2[g] && gCnt[i]<=gCnt2[g-1])
         {
            occurred=0;
            for (j=0; j<g; j++) if (gVal[i]==gVal2[j]) occurred=1;
            if (occurred==0) { gCnt2[g]=gCnt[i]; gVal2[g]=gVal[i]; }
         }
      }
   }
   for (g=0; g<gSize; g++) { gCnt[g]=gCnt2[g]; gVal[g]=gVal2[g]; }

   Drray1(gVal2, gSize);
   Drray1(gCnt2, gSize);
}


// get frequencies in 2-d array
template<class T> void Freq(T **a,int s1,int s2,T *gVal, int *gCnt)
{
   int g, i, j, gSize;
   
   gSize=nGroup(a, s1, s2);
   for (g=0; g<gSize; g++) { gVal[g]=0; gCnt[g]=0; }
   for (g=0; g<gSize; g++)
   {
      if (g==0) gVal[g]=Min(a, s1, s2);
      else if (g>=1)
      {
         for (i=0;i<s1;i++) for (j=0;j<s2;j++)
         {
            if (a[i][j]>gVal[g-1] && (gVal[g]<=0||gVal[g]>a[i][j]))
               gVal[g]=a[i][j];
         }
      }
   }

   for (g=0; g<gSize; g++) for (i=0;i<s1;i++) for (j=0;j<s2;j++)
   {
      if (a[i][j]==gVal[g] || fabs(a[i][j]-gVal[g])<0.0001) ++gCnt[g];
   }
}


// print frequency in 1-d array
template<class T> void prtFreq(ostream &out,T *a, int s)
{
   int g, gSize=nGroup(a, s);
   T *gVal;   Array1(gVal, gSize);
   int *gCnt; Array1(gCnt, gSize);

   Freq(a, s, gVal, gCnt);

   out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6)
      <<"   Value   Count     Percent\n"
      <<"----------------------------\n";
   for (g=0; g<gSize; g++)
   {
      out<<setw(8)<<gVal[g]<<setw(8)<<gCnt[g]
         <<setw(12)<<(double)gCnt[g]/s<<endl;
   }
   out<<"----------------------------\n"
<< "   Total"<<setw(8)<<s<<setw(12)<<(double)Sum(gCnt,gSize)/s
<< endl;

   Drray1(gVal, gSize);
   Drray1(gCnt, gSize);
}


template<class T> void prtFreqHet(ostream &out,T *a, int s)
{
   int g, gSize=nGroup(a, s);
   double het;
   T *gVal;      Array1(gVal, gSize);
   int *gCnt;    Array1(gCnt, gSize);
   int *gCntCum; Array1(gCntCum, gSize);

   Freq(a, s, gVal, gCnt);
   gCntCum[0]=gCnt[0];
   for (g=1; g<gSize; g++) gCntCum[g]=gCntCum[g-1]+gCnt[g];
   het=1.0-SS(gCnt, gSize)/(s*s);

   out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(6)
      <<"                              Cumulative\n"
      <<"   Value   Count     Percent     Percent\n"
      <<"----------------------------------------\n";
   for (g=0; g<gSize; g++)
   {
      out<<setw(8)<<gVal[g]<<setw(8)<<gCnt[g]
         <<setw(12)<<(double)gCnt[g]/s
         <<setw(12)<<(double)gCntCum[g]/s
         <<endl;
   }
   out<<"----------------------------------------"
<< "\n   Total"<<setw(8)<<s<<setw(12)<<(double)Sum(gCnt,gSize)/s
<< "\nHeterogeneity = 1 - SS(Percent) = "<<het<<endl;

   Drray1(gVal, gSize);
   Drray1(gCnt, gSize);
   Drray1(gCntCum, gSize);
}


// print frequency in 2-d array
template<class T> void prtFreq(ostream &out,T **a,int s1,int s2)
{
   int g, gSize=nGroup(a, s1, s2);
   T *gVal;   Array1(gVal, gSize);
   int *gCnt; Array1(gCnt, gSize);

   Freq(a, s1, s2, gVal, gCnt);

   out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4)
      <<"    Value       Count\n"
      <<"---------------------"<<endl;
   for (g=0; g<gSize; g++)
   {
      out<<setw(7)<<gVal[g]<<setw(12)<<gCnt[g]<<endl;
   }
   out<<"---------------------\n"
      <<"    Total"<<setw(12)<<s1*s2<<endl;

   Drray1(gVal, gSize);
   Drray1(gCnt, gSize);
}


// get mean and sd by group in 1-d array
template<class T> void MeanSDByGroup(T *a,T *aGroup, int s,T *gVal,
                   int *gCnt, double *gMean, double *gSD,T *gMin,T *gMax)
{
   int g, i, gSize=nGroup(aGroup, s);

   Freq(aGroup, s, gVal, gCnt);

   for (g=0; g<gSize; g++)
   {
      gMean[g] = 0;
      gSD[g]=0;
      gMin[g]=10000;
      gMax[g]=-10000;

      for (i=0; i<s; i++)
      {
         if (aGroup[i]==gVal[g] || fabs(aGroup[i]-gVal[g])<0.0001)
         {
            gMean[g] += a[i];
            gSD[g] += a[i]*a[i];
            if (a[i]<gMin[g]) gMin[g]=a[i];
            if (a[i]>gMax[g]) gMax[g]=a[i];
         }
      }

      gMean[g] /= gCnt[g];
      gSD[g] = pow( (gSD[g]-gMean[g]*gMean[g]*gCnt[g])/(gCnt[g]-1), 0.5 );
   }
}


// get mean and sd by group in 2-d array
template<class T> void MeanSDByGroup(T **a,T **aGroup,int s1,int s2,T *gVal,
                   int *gCnt, double *gMean, double *gSD)
{
   int g, i, j, gSize=nGroup(aGroup, s1, s2);

   Freq(aGroup, s1, s2, gVal, gCnt);

   for (g=0; g<gSize; g++)
   {
      gMean[g] = 0;
      gSD[g]=0;
   }

   for (g=0; g<gSize; g++)
   {
      for (i=0;i<s1;i++) for (j=0;j<s2;j++)
      {
         if (aGroup[i][j]==gVal[g] || fabs(aGroup[i][j]-gVal[g])<0.0001)
         {
            gMean[g] += a[i][j];
            gSD[g] += a[i][j]*a[i][j];
         }
      }

      gMean[g] /= gCnt[g];
      gSD[g] = pow( (gSD[g]-gMean[g]*gMean[g]*gCnt[g])/(gCnt[g]-1), 0.5 );
   }
}


// print mean and sd by group in 1-d array
template<class T> void prtMeanSDByGroup(ostream &out,T *a,T *aGroup, int s)
{
   int g, gSize=nGroup(aGroup, s);
   T *gVal;
   int *gCnt;
   double *gMean, *gSD;
   T *gMin, *gMax;
   Array1(gVal, gSize);
   Array1(gCnt, gSize);
   Array1(gMean, gSize);
   Array1(gSD, gSize);
   Array1(gMin, gSize);
   Array1(gMax, gSize);

   MeanSDByGroup(a, aGroup, s, gVal, gCnt, gMean, gSD, gMin, gMax);

   out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4);
   out<<"    Group     Count      Mean        SD       Min       Max\n"
      <<"-----------------------------------------------------------\n";
   for (g=0; g<gSize; g++)
   {
      out<<setw(9)<<gVal[g]<<setw(10)<<gCnt[g] 
          <<setw(10)<<gMean[g]<<setw(10)<<gSD[g]
          <<setw(10)<<gMin[g]<<setw(10)<<gMax[g]
          <<endl;
   }
   out<<"-----------------------------------------------------------\n"
      <<"    Total"<<setw(10)<<s
      <<setw(10)<<Mean(a, s)<<setw(10)<<SD(a, s)
      <<setw(10)<<Min(a, s)<<setw(10)<<Max(a, s)
      <<endl;

   Drray1(gVal, gSize);
   Drray1(gCnt, gSize);
   Drray1(gMean, gSize);
   Drray1(gSD, gSize);
   Drray1(gMin, gSize);
   Drray1(gMax, gSize);
}


// print mean and sd by group in 2-d array
template<class T> void prtMeanSDByGroup(ostream &out,T **a,T **aGroup,int s1,int s2)
{
   int g, gSize=nGroup(aGroup, s1, s2);
   T *gVal;
   int *gCnt;
   double *gMean, *gSD;
   Array1(gVal, gSize);
   Array1(gCnt, gSize);
   Array1(gMean, gSize);
   Array1(gSD, gSize);

   MeanSDByGroup(a, aGroup, s1, s2, gVal, gCnt, gMean, gSD);

   out<<"    Group       Count        Mean          SD\n"
      <<"---------------------------------------------" 
      <<endl;
   
   for (g=0; g<gSize; g++)
   {
      out<<setw(9)<<gVal[g]<<setw(12)<<gCnt[g] 
         <<setw(12)<<gMean[g]<<setw(12)<<gSD[g]
         <<endl;
   }
   out<<"---------------------------------------------\n" 
       <<"    Total"<<setw(12)<<s1*s2
       <<setw(12)<<Mean(a, s1, s2)
       <<setw(12)<<SD(a, s1, s2)
       <<endl;

   Drray1(gVal, gSize);
   Drray1(gCnt, gSize);
   Drray1(gMean, gSize);
   Drray1(gSD, gSize);
}


#endif

// Line(date): 2071(9/10/2009 for Windows & Linux), 2128(10/3),2182(12/28).
// 2201(8/18/2010),2210(9/4).
// 2217(4/8/2011),2252(6/30).
// 2478(7/13/2012 for g++44 version)
