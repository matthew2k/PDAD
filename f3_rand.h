/*
PROGRAM: f3_rand.h, 
PURPOSE: All RANDOM NUMBER and cumulative distribution functions
*/


#ifndef F3_RAND_H
#define F3_RAND_H


////////////////////////////////////////////////////////////
// [3] RANDOM NUMBER FUNCTIONS of f3_rand.h
////////////////////////////////////////////////////////////

// uniform random variable function U(0,1) w/ high accuracy
double Unfm()
{
   double u;
redo_u01:
   u=(double)((long long)rand()*RAND_MAX+rand())/(((long long)RAND_MAX+1)*RAND_MAX);
   if (u>0 && u<1) return u; else goto redo_u01;
}

// uniform random variable function U(a,b) w/ high accuracy
double UnfmAny(double a, double b)
{
   double u; if (a>b) { cerr<<"\nERROR: Double value out of range: UnfmAny().\n\n"; exit(1); }
redo_u01:
   u=( (double)((long long)rand()*RAND_MAX+rand())/(((long long)RAND_MAX+1)*RAND_MAX) )*(b-a)+a;
   if (u>=a && u<=b) return u; else goto redo_u01;
}

// standard normal rv function N(0,1) -- taken from DCT pgm Random4f.h 9/17/96
double Norm()
{
   static int ir = 0; static double an = 0; double e, w, v1, v2;

   int newr=-1; double newn=0;
   if (newr == -1)
   {
      if (ir == 0)
      {
         do { v1=2.0*Unfm()-1.0; v2=2.0*Unfm()-1.0; w=v1*v1+v2*v2; } while (w > 1);
         e = sqrt((-2.0*log(w))/w); an = v1*e; ir = 1;
         return v2*e;
      }
      else { ir = 0; return an; }
   }
   else { ir = newr; an = newn; return 0; }
}

// normal rv function N(mu, sigma)
double Norm(double mu, double sigma)
{
   return (sigma*Norm()+mu);
}

// restructed normal rv function N(me,sd) w/ return value=[min,max]
double NormR(double me, double sd, double min,double max)
{
   double d;
   do
   {
      d=Norm(me,sd);
   }
   while (d>max || d<min);
   return d;
}

// standard normal density function
double NormDF(double z)
{
   return INVSQRT2PI*exp(-z*z/2);
}

// standard normal cumulative distribution function, p = Pr(Z_StdNor < z)
double NormCDF(double z)
{
   double T,Z,ANS, p;

   Z=fabs(z)/sqrt(2.0); T=1.0/(1.0+0.5*Z);
   ANS=T*exp(-Z*Z-1.26551223+T*(1.00002368+T*(0.37409196+T*(0.09678418+
      T*(-0.18628806+T*(0.27886807+T*(-1.13520398+T*(1.48851587+
      T*(-0.82215223+T*0.17087277)))))))));
   p=(z>=0.0 ? 0.5*(2-ANS) : 1.-0.5*(2-ANS));
   return p;
}

// -log10(p) of 2-sided normal p-value, p = Pr(X_StdNorm>|z|)
double NormNlp2S(double z)
{
   double nlp, c=exp(z*z/2), logc=z*z/2;

   if (z>=-7. && z<=7.)
   {
      nlp = -log10(2*(1.-NormCDF(fabs(z))));
   }
   else
   {
      nlp = 2*INVSQRT2PI*exp(-z*z/2+logc)/sqrt(1.+z*z);
      nlp = -log10(nlp)+log10(c);
   }

   return nlp;
}

// get nlp=-log10(p) from z byTheRightTailNormal p = Pr(X_StdNorm>z)
double NormNlpR(double z)
{
   double z1,zMax,nlp,c,logc;

   zMax=37.6771; z1=Min(z,zMax); c=exp(z1*z1/2); logc=z1*z1/2;
   if (z1>=-7. && z1<=7.)
   {
      nlp = -log10(1.-NormCDF(z1));
   }
   else if (z1>7.)
   {
      nlp = INVSQRT2PI*exp(-z1*z1/2+logc)/sqrt(1.+z1*z1);
      nlp = -log10(nlp)+log10(c);
      if (z>z1) nlp=nlp+(z-z1)*(NormNlpR(zMax)-NormNlpR(zMax-1));
   }
   else
   {
      nlp=0.;
   }

   return nlp;
}

// inverse standard normal CDF, get z at p = Pr(Z_StdNor < z). remark: this function 
// was downloaded & edited by Dajun Qian from http://www.mathfinance.de/FF/cpplib.html
double NormInv(double p)
{
   static double a[4]={2.50662823884,
                     -18.61500062529,
                      41.39119773534,
                     -25.44106049637};
   static double b[4]={-8.47351093090,
                       23.08336743743,
                      -21.06224101826,
                        3.13082909833};
   static double c[9]={0.3374754822726147,
                       0.9761690190917186,
                       0.1607979714918209,
                       0.0276438810333863,
                       0.0038405729373609,
                       0.0003951896511919,
                       0.0000321767881768,
                       0.0000002888167364,
                       0.0000003960315187};
   double x,z;
   x=p-0.5;
   if (fabs(x)<0.42)
   {
      z=x*x;
      z=x*(((a[3]*z+a[2])*z+a[1])*z+a[0]) / ((((b[3]*z+b[2])*z+b[1])*z+b[0])*z+1.0);
      return(z);
   }
   z=p;
   if(x>0.0) z=1.0-p;
   z=log(-log(z));
   z=c[0]+z*(c[1]+z*(c[2]+z*(c[3]+z*(c[4]+z*(c[5]+z*(c[6]+z*(c[7]+z*c[8])))))));
   if(x<0.0) z=-z;
   return(z);
}


///////////////////////////////////////
// chisq_1of4.RANLIB: random variable (RV) of non-central chi-square (CHISQ) distribution
double rchisq(double df,double nc)
{
   if (df<=0. || nc<0.) { cerr<<"\nERROR: Invalid parameters in rchisq("
      <<df<<","<<nc<<").\n\n"; exit(1); }
   return gennch(df,nc); // from RANLIB: faster than BOOST quantile transformation
}

// chisq_2of4.BOOST: density of nc chisq dist (PDF from xQtl)
double dchisq(double x,double df,double nc)
{
   if (df<=0. || nc<0.) { cerr<<"\nERROR: Invalid parameters in dchisq("
      <<df<<","<<nc<<").\n\n"; exit(1); }
   return (x>=0?pdf(non_central_chi_squared(df,nc),x):1e-50);
}

// chisq_3of4.BOOST: probability of nc chisq dist (CDF pLt from xQtl)
double pchisq(double x,double df,double nc)
{
   if (x<0. || df<=0. || nc<0.) { cerr<<"\nERROR: Invalid parameters in pchisq("
      <<x<<","<<df<<","<<nc<<").\n\n"; exit(1); }
   return Min(cdf(non_central_chi_squared(df,nc), x), 1-1.1e-16);
}

// chisq_4of4.BOOST: quantile of nc chisq dist (inv CDF xQtl from pLt)
double qchisq(double p,double df,double nc)
{
   if (p<0. || df<=0. || nc<0.) { cerr<<"\nERROR: Invalid parameters in qchisq("
      <<p<<","<<df<<","<<nc<<").\n\n"; exit(1); }
   return quantile(non_central_chi_squared(df,nc), p);
}


///////////////////////////////////////
double rgamma(double a,double b) // gamma_1of4.RANLIB: rv w/ a=shape, b=scale=1/rate
{
   return gengam(1./b, a);
}

double dgamma(double x,double a,double b) // gamma_2of4.BOOST: density=PDF
{
   return (x>=0?gamma_p_derivative(a,x/b)/b:1e-50);
}

double pgamma(double x,double a,double b) // gamma_3of4.BOOST: probability=CDF
{
   gamma_distribution<> gamma_dist(a,b); return Min(cdf(gamma_dist,x), 1-1.1e-16);
}

double qgamma(double p,double a,double b) // gamma_4of4.BOOST: quantile=invCDF
{
   gamma_distribution<> gamma_dist(a,b); return quantile(gamma_dist,p);
}


///////////////////////////////////////
double rlnorm(double mu,double sigma) // lnorm_1of4: rv w/ mu+sigma=(mean+sd of normal dist)
{
   return exp(Norm()*sigma+mu);
}

double dlnorm(double x,double mu,double sigma) // lnorm_2of4: density=PDF
{
   double z;
   z=(log(x)-mu)/sigma; return Max(INVSQRT2PI*exp(-z*z/2)/(x*sigma), 1e-50);
}

double plnorm(double x,double mu,double sigma) // lnorm_3of4: probability=CDF
{
   return NormCDF((log(x)-mu)/sigma);
}

double qlnorm(double p,double mu,double sigma) // lnorm_4of4: quantile=invCDF
{
   return exp(NormInv(p)*sigma+mu);
}


///////////////////////////////////////
double rweibull(double a,double b) // weibull_1of4: rv w/ a=shape, b=scale
{
   return b*pow(-log(Unfm()),1./a);
}

double dweibull(double x,double a,double b) // weibull_2of4.BOOST: density=PDF
{
   return (x>=0?pdf(weibull(a,b),x):1e-50);
}

double pweibull(double x,double a,double b) // weibull_3of4.BOOST: probability=CDF
{
   weibull_distribution<> weibull_dist(a,b); return Min(cdf(weibull_dist,x), 1-1.1e-16);
}

double qweibull(double p,double a,double b) // weibull_4of4: quantile=invCDF
{
   //weibull_distribution<> weibull_dist(a,b); return quantile(weibull_dist,p);
   return b*pow(-log(1-p),1./a); // seems same accuracy as BOOST version above
}


///////////////////////////////////////
double rgumbel(double a,double b) // gumbel_1of4: rv w/ a=location=anyDbl, b=scale>0
{
   return a-b*log(-log(Unfm()));
}

double dgumbel(double x,double a,double b) // gumbel_2of4: density=PDF
{
   return exp(-(x-a)/b-exp(-(x-a)/b))/b;
}

double pgumbel(double x,double a,double b) // gumbel_3of4: probability=CDF
{
   return exp(-exp(-(x-a)/b));
}

double qgumbel(double p,double a,double b) // gumbel_4of4: quantile=invCDF
{
   return a-b*log(-log(p));
}


///////////////////////////////////////
double pgpareto(double x,double a,double b) // gpareto: probability=CDF
{
   double p;
   if      (a>1e-10)  p=((x>=0 && x<(b/a))?(1.-pow(1.-a*x/b,1/a)):1.-1e-50);
   else if (a<-1e-10) p=1.-pow(1.-a*x/b, 1/a);
   else               p=1.-exp(-x/b);
   return p;
}


// exponential rv function w/ mean theta
double Expn(double theta)
{
   return -log(Unfm())*theta;
}

// (naive) gamma rv w/ int shape a & scale b (me=a*b, sd=a*(b*b))
double Gamma_int(int a,double b) // modified from Xuexia subfunction7.h on 7/19/2011
{
   int j; double am,e,s,v1,v2,x,y;

   if(a<1) { cerr<<"\nERROR: Shape parameter must >= 1.\n\n"; exit(1); }
   else if(a<6) // assume gamma=sum_exponential
   {
      x=1.; for(j=0;j<a;j++) x=x*Unfm();
      x=-log(x);
   }
   else
   {
      do
      {
         do
         {
            do
            {
               v1=Unfm(); v2=2.0*Unfm()-1.0;
            }
            while(v1*v1+v2*v2>1.0);

            y=v2/v1; am=a-1; s=sqrt(2.0*am+1.0); x=s*y+am;
         }
         while(x<=0.0);

         e=(1.0+y*y)*exp(am*log(x/am)-s*y);
      }
      while(Unfm()>e);
   }

   return(x*b);
}

/*
This gamma function generates rv using an acceptence-rejection algorithm based on:
Marsaglia G and Tsang WW (2000) A Simple Method for generating gamma variables, ACM Transactions 
on Mathematical Software, 26(3):363-72.

Modified from double gsl_ran_gamma (const gsl_rng * r, const double a, const double b) in gsl-1.14
(downloaded from http://www.gnu.org/s/gsl/, saved to C:\z_soft\gsl\gsl-1.14\randist\gamma.c).

Example runs for a=0.184, b=0.320 (mean=a*b=0.058880, sd=sqrt(a)*b=0.137265):
      n,     mean,       sd,   median,      min,      max
   1000, 0.055129, 0.125657, 0.004161, 0.000000, 0.987630
  10000, 0.060440, 0.140935, 0.005466, 0.000000, 2.177013
 100000, 0.058795, 0.137968, 0.004821, 0.000000, 2.614496
*/
// gamma rv w/ double shape a & double scale b (me=a*b, sd=sqrt(a)*b)
double Gamma(double a, double b)
{
   double d,c,x,v,u;

   if (a<=0.||b<=0.) { cerr<<"\nERROR: rv function Gamma("<<a<<","<<b<<") is invalid.\n\n"; exit(1); }
   else if (a<1)
   {
      u=Unfm();
      return Gamma(a+1.,b)*pow(u,1./a);
   }
   else
   {
      d=a-1./3.; c=1./sqrt(9.*d);
      do
      {
         do
         {
            x=Norm(); v=1.+c*x;
         }
         while(v<=0.);
         v=v*v*v; u=Unfm();
      }
      while ( (u>1.-.0331*(x*x)*(x*x)) && (log(u)>0.5*x*x+d*(1.-v+ log(v))));
      return b*d*v;
   }
}

// mixed gamma rv w/ p0 for 1 mass at x0 and 1-p0 for for gamma(a,b)
double Gamma_1mass(double p0,double x0,double a,double b)
{
   double u,x;

   u=Unfm();
   if (u<p0) x=x0;
   else x=Gamma(a,b);
   return x;
}


// beta rv function w/ a & b (UNFINISHED)
double Beta(double a, double b)
{
   a=1.0; b=1.0;
   return 1.0;
}

// random integer b/ iMin and iMax
int RdInt(int iMin, int iMax)
{
   int iRg=iMax-iMin+1, iRtn=-1;
   if (iRg<=0)
   {
      cerr<<"\nERROR: Input values out of range, RdInt()\n\n"; exit(1);
   }
   else if (iRg==1) return iMin;
redo_rand:
   if      (iRg>1 && iRg<=1000) iRtn=iMin+rand()%iRg;
   else if (iRg>1000) iRtn=iMin+(int)((double)iRg*Unfm());
   if (iRtn<iMin||iRtn>iMax) goto redo_rand;
   return iRtn;

//   return iMin + (int)((double)(iMax-iMin+1)*rand()/(RAND_MAX+1.0));
}

// random integer except 1 value
int RdIntEx1(int iMin, int iMax, int iEx)
{
   int iRtn=RdInt(iMin, iMax-1);
   if (iRtn>=iEx) iRtn++;
   return iRtn;
}

// random integer except 2 values
int RdIntEx2(int iMin, int iMax, int iEx1, int iEx2)
{
   int iRtn;

   do
   {
      iRtn=RdInt(iMin, iMax);
   }
   while (iRtn==iEx1||iRtn==iEx2);
   return iRtn;
}

// random integer from 0,1,...,n-1 w/ pr[] sum of 1
int RdIntPr(double *pr, int n)
{
   int i, iRtn; double cumPr, u01;

   iRtn=0; cumPr=pr[0]; u01=Unfm();
   for (i=1; i<n; i++) 
   {
      if (cumPr < u01) iRtn=i; else break;
      cumPr = cumPr + pr[i];
   }
   return iRtn;
}

// random integer from 0,1,...,n-1, and check if pr[] sum is 1
int RdIntPrChk(double *pr, int n)
{
   int i, iRtn=-1;
   double u01, cumPr;

   for (i=0; i<n; i++) if (pr[i]<-.0 || pr[i]>1.0)
   { cerr<<"\nERROR: Double value out of range: RdIntPrChk().\n\n"; exit(1); }
   cumPr=Sum(pr, n); if (cumPr<0.999 || cumPr>1.001)
   { cerr<<"\nERROR: Double value out of range: RdIntPrChk().\n\n"; exit(1); }

   u01=Unfm(); cumPr=0;
   for (i=0; i<n; i++) 
   {
      if (cumPr<u01) iRtn=i; else break;
      cumPr=cumPr + pr[i];
   }
   return iRtn;
}

// random integer from 0,1,...,n-1 for any prob[] sum
int RdIntPrAny(double *pr, int n)
{
   int i, iRtn; double sumPr,cumPr, u01;

   sumPr=Sum(pr, n); iRtn=0; cumPr=pr[0]; u01=Unfm()*sumPr;
   for (i=1; i<n; i++) 
   {
      if (cumPr < u01) iRtn=i; else break;
      cumPr = cumPr + pr[i];
   }
   return iRtn;
}

// random integer from 0,1,...,n-1 c.o pr[]+sumPr (50% faster than xxxAny(...))
int RdIntPrSum(double *pr, double sumPr, int n)
{
   int i, iRtn; double cumPr, u01;

   iRtn=0; cumPr=pr[0]; u01=Unfm()*sumPr;
   for (i=1; i<n; i++) 
   {
      if (cumPr < u01) iRtn=i; else break;
      cumPr = cumPr + pr[i];
   }
   return iRtn;
}

// sample m random integers from 0:n-1 w/ prob[] & replacement
void RdIntPrWR(double *prob, int n, int *iRtn, int m)
{
   int i, j;
   double u01; 
   double *cumProb; Array1(cumProb, n); // cum prob density
   
   // check the probability density values
   for (i=0; i<n; i++) if (prob[i] < -.0 || prob[i] > 1.0)
   { cerr<<"\nERROR: Double value out of range: RdIntPrWR().\n\n"; exit(1); }

   // get the cumulative probabilities at 0,1,...n-1
   cumProb[0] = prob[0];
   for (i=1; i<n; i++) cumProb[i] = cumProb[i-1] + prob[i];
   if (cumProb[n-1] < 0.99999 || cumProb[n-1] > 1.00001)
      { cerr<<"\nERROR: Double value out of range: RdIntPrWR().\n\n"; exit(1); }

   // return m integers
   for (j=0; j<m; j++)
   {
      u01 = Unfm();
      if (u01<cumProb[0]) iRtn[j] = 0;
      else
      {
         for (i=1; i<n; i++) 
         {
            if (u01>=cumProb[i-1] && u01<cumProb[i]) { iRtn[j]=i; break; }
         }
      }
   }

   Drray1(cumProb, n);
}

// get 2 random integers w/ replacement
void get2RdIntWR(int iMin, int iMax, int &i1, int &i2)
{
   i1=RdInt(iMin, iMax);
   i2=RdInt(iMin, iMax);
}

// get 2 random integers w/o replacement
void get2RdIntNR(int iMin, int iMax, int &i1, int &i2)
{
   i1=RdInt(iMin, iMax);
   i2=RdInt(iMin, iMax-1);
   if (i2>=i1) i2++;
}

// get 4 random integers w/o replacement
void get4RdIntNR(int iMin, int iMax, int &i1, int &i2, int &i3, int &i4)
{
   int *iRd; Array1(iRd,4);
   getMltRdIntNR(iMin, iMax, iRd, 4);
   i1=iRd[0]; i2=iRd[1]; i3=iRd[2]; i4=iRd[3];
   Drray1(iRd,4);
}

//get sz ramdom integers w/ replacement
void getMltRdIntWR(int iMin, int iMax, int *iRtn, int sz)
{
   for (int i=0; i<sz; i++) iRtn[i]=RdInt(iMin, iMax);
}

// get bootstrapped sample vb[s] from v[s]
void bootstrap(double *v,int s, double *vb)
{
   int i,iRd;
   for (i=0;i<s;i++) { iRd=RdInt(0,s-1); vb[i]=v[iRd]; }
}


//get sz ramdom integers within (iMin, iMax) w/o replacement
void getMltRdIntNR(int iMin, int iMax, int *iRtn, int sz)
{
   int i, iRg=iMax-iMin+1;
   if (iMin>iMax) { cerr<<"\nERROR: Integer value out of range.\n\n"; exit(1); }
   if (sz>iRg) { cerr<<"\nERROR: Integer value out of range.\n\n"; exit(1); }
   int *perm; Array1(perm, iRg);

   asgnArray012Perm(perm, iRg);
   for (i=0; i<sz; i++) iRtn[i]=perm[i]+iMin;
   Drray1(perm, iRg);
}

// assign integer array with 0, 1, ...,s-1
void asgnArray012(int *a, int s)
{
   for (int i=0;i<s;i++) a[i]=i;
}

// assign integer array w/ permutation of 0, 1, ...,s-1
void asgnArray012Perm(int *perm, int s)
{
   int i, iRd;

   for (i=0; i<s; i++) perm[i]=i;
   for (i=0; i<s; i++)
   {
      iRd=RdInt(0, s-1-i);
      if (iRd!=0) Swap(&perm[i], &perm[i+iRd]);
   }
}

// poisson random variable with mean lambda
int Pois(double lambda) 
{
   int iRtn = 0;
   double sumNLU = 0; // sum of negative log(U(0,1))
   while (sumNLU < lambda) 
   {
      sumNLU = sumNLU-log(Unfm());
      if (sumNLU < lambda) iRtn++;
   }
   return iRtn;
}

// truncated poisson random function
int PoisT(double lambda, int iMin, int iMax)
{
   int iRtn=Pois(lambda);
   if (iRtn<iMin) iRtn=iMin;
   else if (iRtn>iMax) iRtn=iMax;
   return iRtn;
}

// restructed poisson random function
int PoisR(double lambda, int iMin, int iMax)
{
   int iRtn;
   do
   {
      iRtn=Pois(lambda);
   }
   while (iRtn>iMax || iRtn<iMin);
   return iRtn;
}

// binomial generator w/ probability p1
int Binm(double p1, int i1, int i2)
{
   if (p1 < -.00001 || p1 > 1.00001) {
      cerr<<"\nERROR: Double value out of range: Binm().\n\n"; exit(1); }
   if (Unfm() < p1) return i1; else return i2;
}

// binomial generator of i1 and i2 w/ any probability desity p1 and p2
int BinmAnyPr(double p1, double p2, int i1, int i2)
{
   if (p1 < 0 || p2 < 0) {
      cerr<<"\nERROR: Double value out of range: BinmAnyPr().\n\n"; exit(1); }
   if (Unfm() < p1/(p1+p2)) return i1; else return i2;
}

// geometric random number
int Geom(double p)
{
   int iRtn = 0;
   double sumPdf = p; // sum of geometric pdf(p)
   double u01 = Unfm();

   while (sumPdf < u01) { iRtn++; sumPdf += p * pow(1-p, iRtn); }
   return iRtn;
}

// restricted geometric random number
int GeomR(double p, int iMin, int iMax)
{
   int iRtn;

   do { iRtn=Geom(p); }
   while (iRtn>iMax || iRtn<iMin);
   return iRtn;
}

// poisson density function
double PoisDF(int x, double lambda)
{
   int i; double p_jk; double p;
   p_jk=exp(-lambda); for (i=1; i<=x; i++) p_jk *= (lambda/i); p=p_jk;
   return p;
}

// poisson cumulative distribution function
double PoisCDF(int x, double lambda)
{
   int i; double p;
   p=0; for (i=0; i<=x; i++) p += PoisDF(i, lambda);
   return p;
}

// poisson CDF trucated at pmax
double PoisCDFT(int x, double lambda, double pmax)
{
   int i; double p;
   p=0; for (i=0; i<=x; i++) { p += PoisDF(i, lambda); if (p>=pmax) break; }
   return p;
}

// multinomial random density function (i.e., probability mass function)
void MltnmRdDF(int N, double *df)
{
   int i; double *cdf; Array1(cdf,N);
   for (i=0;i<N-1;i++) cdf[i]=Unfm(); sortArrayAsc(cdf,N-1);
   for (i=0;i<N;i++)
   {
      if      (i==0)           df[i]=cdf[i];
      else if (i>=1 && i<=N-2) df[i]=cdf[i]-cdf[i-1];
      else if (i==N-1)         df[i]=1.-cdf[i-1];
   }
   Drray1(cdf,N);
}



#endif

// Line(date): 506(9/10/2009 for Windows & Linux),504(12/28)
// 535(4/8/2011),550(6/30)
