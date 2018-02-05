/*
PROGRAM: testpdad.h, DOB=1/7/2013                                    LAST=1/20/2014
PURPOSE: Get p-value accuracy by method, and run other analyses
*/

// UTILITY FUNCTIONS
void getxta(double mp,double vp,double *pta,int nDI,int nT, double **xta);       // f0a
void getxt2a_w2d(double mp,double vp,double *pta,int nDI,int nT,int nLP,         // f0a_1
   double w12,double r12,int iPST, double ***xt2a);
void getab_mv(int di,double m,double v, double &a,double &b);                    // f0b
void getwb_mv(double m,double v, double &a,double &b);                           // f0b_1
double dmv_wb(double *par,int n1,double *sta,int n2);                            // f0b_1a
void simxs_pd(int di,double a,double b,int n, double *xs);                       // f0c
void simxs_ws2d(int di1,double a1,double b1,int di2,double a2,double b2,int n,   // f0c_1
   double w12,double r12, double *xs);
void simxs_bmd(int di1,double a1,double b1,int di2,double a2,double b2,int n,    // f0c_2
   double pfstd,double dbm, double *xs);
double bisect_xt_minb(double *xs,int n,double pt, double ef,double ep);          // f0d
double bias_pfpt(double pf,double pt);                                           // f0e

// RESULTS FUNCTIONS
void pw_try(int pd,double mp,double vp,int n,int nR,double xmi,double xma,       // try
   int iPST,char *fn);
void pbbm(double mp,double vp,char *pdg,char *ng,int nT,int nR,int nLP,          // f1
    int iPST,char *fn);

void tsbm(double mp,double vp,char *pdg,char *ng,int nR, int iPST,char *fn);     // f2

void nbm(int mtd,int pd,int ana,int nT,double mp,double vp,int nR,               // f3
   int iPST,char *fn,int o1a2,int aTot,int &aDn);
void bisect_nbm(int mtd,double *xs,int nSM,int nST,int nR,double xt,double qt,   // f3a
   double bpctg,double bu,double ef, int &nMi,int &nMa,int &pn);
double getblmt(int mtd,int n,double *xs,int nSM,int nR,double xt,double qt,      // f3a_1
   double bpctg, double *xw,double *bw);
void nbm_mrun(double mp,double vp,int nR,int nA, int iPST,char *fn);             // f3b


// f0a: get xta[nDI][nT]=trueStat co mp+vp+pta[nT]
void getxta(double mp,double vp,double *pta,int nDI,int nT, double **xta)
{
   int j,t; double a,b;
   for (j=0;j<nDI;j++)
   {
      getab_mv(j+1,mp,vp, a,b); // f0b
      if      (j==0) for (t=0;t<nT;t++) xta[j][t]=qchisq(1.-pta[t],a,b);
      else if (j==1) for (t=0;t<nT;t++) xta[j][t]=qgamma(1.-pta[t],a,b);
      else if (j==2) for (t=0;t<nT;t++) xta[j][t]=qlnorm(1.-pta[t],a,b);
      else if (j==3) for (t=0;t<nT;t++) xta[j][t]=qweibull(1.-pta[t],a,b);
      else if (j==4) for (t=0;t<nT;t++) xta[j][t]=qgumbel(1.-pta[t],a,b);
      else { cerr<<"\nERROR: getxta() failed for unknown pd="<<j+1<<".\n\n"; exit(1); }
   }
}

// f0a_1: getTrueStats xt2a[nDI,nDI,nT] co mp+vp+nDI+pta[nT] for wtSum2Dist_nLP+w12+r12
void getxt2a_w2d(double mp,double vp,double *pta,int nDI,int nT,int nLP,
   double w12,double r12,int iPST, double ***xt2a)
{
   int k,l,t; double *aa,*ba,*xs;
   Array1(aa,nDI); Array1(ba,nDI); Array1(xs,nLP);

   for (k=0;k<nDI;k++) getab_mv(k+1,mp,vp, aa[k],ba[k]); // f0b
   for (k=0;k<nDI;k++) for (l=0;l<nDI;l++)
   {
      cout<<"TIME: getxt2a_w2d() for pd1="<<k+1<<", pd2="<<l+1<<",   time="<<flush;
      simxs_ws2d(k+1,aa[k],ba[k],l+1,aa[l],ba[l],nLP,w12,r12, xs); // f0c_1
      for (t=0;t<nT;t++) xt2a[k][l][t]=bisect_xt_minb(xs,nLP,pta[t],1e-4,1e-4); // f0d
      prtHmsDttm(cout,iPST); // cout<<setw(11)<<""; prtMdyHms(cout);
   }
   Drray1(aa,nDI); Drray1(ba,nDI); Drray1(xs,nLP);
}


// f0b: get dist par a/b for a specified di=distId
void getab_mv(int di,double m,double v, double &a,double &b)
{
   if (m>0. && v>0.)
   {
      if      (di==1) { if (m>v/4 && m<=v/2) { a=2*m-v/2; b=v/2-m; } else { a=m; b=0.; } }
      else if (di==2) { a=m*m/v; b=v/m; }
      else if (di==3) { a=log(m/sqrt((v+m*m)/(m*m))); b=sqrt(log((v+m*m)/(m*m))); }
      else if (di==4) { getwb_mv(m,v, a,b); } // f0b_1
      else if (di==5) { b=sqrt(6*v)/PI; a=m-0.5772156649*b; }
      else { cerr<<"\nERROR: getab_mv() failed for unknown di="<<di<<".\n\n"; exit(1); }
   }
   else { cerr<<"\nERROR: getab_mv() failed for invalid m="<<m<<", v="<<v<<".\n\n"; exit(1); }
}

// f0b_1: get 2 para a+b of weibull dist co me+var. note: hookeOftenReq>=1000Iter
void getwb_mv(double m,double v, double &a,double &b)
{
   double *par,*sta,*dlt,*pmi,*pma;
   Array1(par,2); Array1(sta,2); Array1(dlt,2); Array1(pmi,2); Array1(pma,2);

   par[0]=-log(2.)/log(1.-sqrt(v/3)/m); par[1]=m/exp(lgamma(1.+1./par[0]));
   sta[0]=m; sta[1]=v; dlt[0]=1; dlt[1]=1; pmi[0]=0; pmi[1]=0;
   pma[0]=Max(2*par[0],100.); pma[1]=Max(2*par[1],100.);
   hooke_minf(par,2,sta,2,dlt,pmi,pma,0.5,1e-6,1000,dmv_wb); // pdad_f3a
   a=par[0]; b=par[1];

   Drray1(par,2); Drray1(sta,2); Drray1(dlt,2); Drray1(pmi,2); Drray1(pma,2);
}

// f0b_1a: get fo=difOfMeanAndVar btw previous sta[] and m+v from weibull(par[])
double dmv_wb(double *par,int n1,double *sta,int n2)
{
   double a,b,m,v;
   a=par[0]; b=par[1]; m=b*exp(lgamma(1.+1./a)); v=b*b*exp(lgamma(1.+2./a))-m*m;
   return pow(m-sta[0],2.)+pow(v-sta[1],2.);
}


// f0c: simu xs[]=statData for di=1~5=chisq/gamma/lnorm/weibull/gumbel
void simxs_pd(int di,double a,double b,int n, double *xs)
{
   int i;
   if      (di==1) for (i=0;i<n;i++) xs[i]=rchisq(a,b);
   else if (di==2) for (i=0;i<n;i++) xs[i]=rgamma(a,b);
   else if (di==3) for (i=0;i<n;i++) xs[i]=rlnorm(a,b);
   else if (di==4) for (i=0;i<n;i++) xs[i]=rweibull(a,b);
   else if (di==5) for (i=0;i<n;i++) xs[i]=rgumbel(a,b);
   else { cerr<<"\nERROR: simxs_pd() failed for unknown di="<<di<<".\n\n"; exit(1); }
}

// f0c_1: simu xs[]=statData for wtSum2Dist
void simxs_ws2d(int di1,double a1,double b1,int di2,double a2,double b2,int n,
   double w12,double r12, double *xs)
{
   int i; double w1,w2,x1,x2;

   w1=w12+(1-w12)*r12; w2=(1-w12)*sqrt(1-r12*r12);
   for (i=0;i<n;i++)
   {
      x1=( di1==1?rchisq(a1,b1):( di1==2?rgamma(a1,b1):( di1==3?rlnorm(a1,b1):
         ( di1==4?rweibull(a1,b1):rgumbel(a1,b1) ) ) ) );
      x2=( di2==1?rchisq(a2,b2):( di2==2?rgamma(a2,b2):( di2==3?rlnorm(a2,b2):
         ( di2==4?rweibull(a2,b2):rgumbel(a2,b2) ) ) ) );
      xs[i]=w1*x1+w2*x2;
   }
}

// f0c_2: simu xs[]=statData for bimodalDist (p1st=pFstDist, dbm=difBtwMeanOf2Dist)
void simxs_bmd(int di1,double a1,double b1,int di2,double a2,double b2,int n,
   double pfstd,double dbm, double *xs)
{
   int i;

   for (i=0;i<n;i++)
   {
      if (Unfm()<pfstd)
      {
         xs[i]=( di1==1?rchisq(a1,b1):( di1==2?rgamma(a1,b1):( di1==3?rlnorm(a1,b1):
            ( di1==4?rweibull(a1,b1):rgumbel(a1,b1) ) ) ) );
      }
      else
      {
         xs[i]=dbm+( di2==1?rchisq(a2,b2):( di2==2?rgamma(a2,b2):( di2==3?rlnorm(a2,b2):
            ( di2==4?rweibull(a2,b2):rgumbel(a2,b2) ) ) ) );
      }
   }
}


// f0d: use bisection to get pn=(root xt) co xs[]+pt+minBias
double bisect_xt_minb(double *xs,int n,double pt, double ef,double ep)
{
   int it; double pa,pb,pn,fa,fb,fn,d1; // pa/pb/pn=parLt/Rt/New, ep/ef=epsilonPar/Fun

   pa=Min(xs,n); d1=1.-0.5/n; fa=sign(-log10(d1)+log10(pt))*bias_pfpt(d1,pt); // f0e
   pb=Max(xs,n); d1=0.5/n;    fb=sign(-log10(d1)+log10(pt))*bias_pfpt(d1,pt); // f0e

   if (fa>0 || fb<0) { cerr<<"\nERROR: bisect_xt_minb() failed for invalid increasing"
      <<" function: f("<<pa<<")="<<fa<<", f("<<pb<<")="<<fb<<".\n\n"; exit(1); }
   else if (fa>=-ef || fb<=ef) { if (-fa<=fb) { pn=pa; fn=fa; } else { pn=pb; fn=fb; } }
   else
   {
      it=0; while (fabs(pb-pa)>ep && it<100)
      {
         pn=(pa+pb)/2; d1=1.-Percentage(xs,n,pn);
         fn=sign(-log10(d1)+log10(pt))*bias_pfpt(d1,pt); it++;
         if      (fn<-ef)  { pa=pn; fa=fn; }
         else if (fn> ef)  { pb=pn; fb=fn; }
         else if (fn>0)    { if ( fn>-fa) { pn=pa; fn=fa; } break; }
         else              { if (-fn> fb) { pn=pb; fn=fb; } break; }
      }
   }
   return pn;
}


// f0e: get scaled log10 bias co pf+pt
double bias_pfpt(double pf,double pt)
{
   return fabs(log10(pf/pt))/(-log10(Min(pt,0.9)));
}


// try: try print of fitted p-values and weight details
void pw_try(int pd,double mp,double vp,int n,int nR,double xmi,double xma, int iPST,char *fn)
{
   int nDI,r,k,l,nd,i, *di; double xt,pt,a1,b1,a2,b2,pf,fbst, **dp2,*xs,*w; nDI=5;
   ofstream out,ou2; char fn2[100]=""; strcpy_newE(fn2,fn,"dat");
   openOF(out,fn); out<<setiosflags(ios::fixed|ios::showpoint);
   openOF(ou2,fn2); ou2<<setiosflags(ios::fixed|ios::showpoint);
   Array2(dp2,nDI,4); Array1(xs,n); Array1(di,nDI); Array1(w,nDI);

   for (k=0;k<nDI;k++) getab_mv(k+1,mp,vp, dp2[k][0],dp2[k][1]); // f0b
   out<<"(1) Input parameters: pd="<<pd<<", mp="<<setprecision(3)<<mp<<", vp="<<vp
      <<", xmi="<<xmi<<", xma="<<xma<<"\n\n"
      <<"(2) Estimated parameters:\n"<<setprecision(4);
   for (k=0;k<nDI;k++)
   {
      out<<"   "<<(k==0?"chisq":(k==1?"gamma":(k==2?"lnorm":(k==3?"weibull":"gumbel"))))
         <<"("<<dp2[k][0]<<", "<<dp2[k][1]<<")\n";
   }

   for (r=0;r<nR;r++)
   {
      xt=xmi+((double)r/(nR-1))*(xma-xmi);
      pt=1.-getp_dp(xt,dp2[pd-1],pd); // pdad_f1e
      if (pd>=1 && pd<=nDI) simxs_pd(pd,dp2[pd-1][0],dp2[pd-1][1],n, xs); // f0c
      else
      {
         k=RdInt(0,nDI-1); l=RdInt(0,nDI-1);
         a1=dp2[k][0]; b1=dp2[k][1]; a2=dp2[l][0]; b2=dp2[l][1];
         simxs_ws2d(k+1,a1,b1,l+1,a2,b2,n,0.5,0.5, xs); // f0c_1
      }
      if (r==0)
      {
         out<<"\n(3) For replicate "<<r+1<<", xs["<<n<<"]=\n"; prtArray1BR(out,xs,n,10);
      }

      pf=p_dad_detail(xs,n,xt,"1,2,3,4,5", nd,di,w,fbst); // pdad_f1a
      if (r==0)
      {
         out<<"\n"<<setw(4)<<"rep"<<setw(11)<<"xt"<<setw(12)<<"pt";
         for (k=0;k<nDI;k++) out<<setw(3)<<"di"<<k+1;
         for (k=0;k<nDI;k++) out<<setw(6)<<"w"<<k+1;
         out<<setw(10)<<"fbst"<<setw(12)<<"p_dad"<<setw(8)<<"bias"<<endl;
      }
      out<<setw(4)<<r+1<<setw(11)<<setprecision(6)<<xt<<setw(12)<<setprecision(8)<<pt
         <<setprecision(3);
      for (k=0;k<nDI;k++) { if (k<nd) out<<setw(4)<<di[k]; else out<<setw(4)<<"na"; }
      for (k=0;k<nDI;k++) { if (k<nd) out<<setw(7)<<w[k]; else out<<setw(7)<<"na"; }
      out<<setw(10)<<fbst<<setw(12)<<setprecision(8)<<pf
         <<setw(8)<<setprecision(3)<<bias_pfpt(pf,pt)<<endl;

      ou2<<setprecision(6)<<xt; for (i=0;i<n;i++) ou2<<","<<xs[i]; ou2<<endl;
   }

   Drray2(dp2,nDI,4); Drray1(xs,n); Drray1(di,nDI); Drray1(w,nDI); out.close(); ou2.close();
   cout<<"NOTE: pw_try() done, results saved to "<<fn<<" and "<<fn2<<", time used ";
   prtHms(cout,iPST); cout<<endl;
}


// f1: compute P-value & bias By 4 Methods and n+pd+stat+rep as tabulated results
void pbbm(double mp,double vp,char *pdg,char *ng,int nT,int nR,int nLP, int iPST,char *fn)
{
   cout<<"TIME: pbbm() started:"<<setw(19)<<""<<"time="; prtHmsDttm(cout,iPST);
   int nPD,nN,nDI,nM,nRT,nMTD, iN,j,k,l,t,r,n,iRT,m,i1, *pda,*na;
   double xt,p1,d1, *aa,*ba,*pta,**xta,***xt2a,*da, *xs,****ps,****bs;
   int nd,*di; double fbst,*w;
   ofstream out,ou2; char fn2[100]=""; strcpy_newE(fn2,fn,"dat");
   openOF(out,fn); out<<setiosflags(ios::fixed|ios::showpoint);
   openOF(ou2,fn2); ou2<<setiosflags(ios::fixed|ios::showpoint);

   Array1(pda,10); Array1(na,10); ipt_ia(pdg,pda,nPD); ipt_ia(ng,na,nN);
   nDI=5; nM=na[nN-1]; nRT=nR*nT; nMTD=4;
   Array1(aa,nDI); Array1(ba,nDI); Array1(pta,nT);
   Array2(xta,nDI,nT); Array3(xt2a,nDI,nDI,nT); Array1(da,nDI*nDI*nT);
   Array1(xs,nM); Array4(ps,nMTD,nN,nPD,nRT); Array4(bs,nMTD,nN,nPD,nRT);
   Array1(di,nDI); Array1(w,nDI);

   // 1: get na[nN]+pda[nPD]+aa/ba[nDI]+pta[nT]+xta[nDI,nT]+xt2a[nDI,nDI,nT] asDescStats
   for (j=0;j<nPD;j++) pda[j]=(j==nPD-1?11:j+1);
   for (k=0;k<nDI;k++) getab_mv(k+1,mp,vp, aa[k],ba[k]); // f0b
   for (t=0;t<nT;t++) pta[t]=(t==0?0.5:(t==1?0.05:pow(10.,-t)));
   getxta(mp,vp,pta,nDI,nT, xta); // f0a
   getxt2a_w2d(mp,vp,pta,nDI,nT,nLP,0.5,0.5,iPST, xt2a); // f0a_1

   // 2: get and print p-values and biases by n+pd+xt
   // print variable names for out+ou2:
   out<<"(1) P-values & biases for "<<nN<<" sample size"<<(nN>=2?"s":"")<<", "<<nDI+1
      <<" population distributions and "<<nT<<" true p-values\n"
      <<"(each row contains estimated p-values and biases for a single analysis):\n";
   out<<setw(4)<<"n"<<setw(4)<<"pd"<<setw(8)<<"p_true"<<setw(8)<<"stat"
      <<setw(9)<<"p_pm"<<setw(9)<<"p_gum"<<setw(9)<<"p_gpd"<<setw(9)<<"p_dad"
      <<setw(7)<<"b_pm"<<setw(7)<<"b_gum"<<setw(7)<<"b_gpd"<<setw(7)<<"b_dad"<<endl;

   ou2<<setw(4)<<"n"<<setw(4)<<"pd"<<setw(8)<<"stat"<<setw(5)<<"rep"<<setw(7)<<"q_true"
      <<setw(7)<<"q_pm"<<setw(7)<<"q_gum"<<setw(7)<<"q_gpd"<<setw(7)<<"q_dad"
      <<setw(7)<<"b_pm"<<setw(7)<<"b_gum"<<setw(7)<<"b_gpd"<<setw(7)<<"b_dad";
   for (k=0;k<nDI;k++) ou2<<setw(6)<<"w"<<k+1; ou2<<setw(7)<<"w_sum"<<endl;

   // compute p-values and biases
   for (j=0;j<nPD;j++) for (r=0;r<nR;r++)
   {
      if (j<=nPD-2) simxs_pd(pda[j],aa[j],ba[j],nM, xs); // f0c
      else
      {
         k=RdInt(0,nDI-1); l=RdInt(0,nDI-1);
         simxs_ws2d(k+1,aa[k],ba[k],l+1,aa[l],ba[l],nM,0.5,0.5, xs); // f0c_1
      }
      for (iN=0;iN<nN;iN++) for (t=0;t<nT;t++)
      {
         n=na[iN]; iRT=r*nT+t; xt=(j==nPD-1?xt2a[k][l][t]:xta[j][t]);
         for (m=0;m<nMTD;m++)
         {
            if      (m==0) p1=1.-Percentage(xs,n,xt);
            else if (m==1) p1=p_gum(xs,n,xt); // pdad_f2a
            else if (m==2) p1=p_gpd(xs,n,xt); // pdad_f2b
            else if (m==3) p1=p_dad_detail(xs,n,xt,"1,2,3,4,5", nd,di,w,fbst); // pdad_f1a
            ps[m][iN][j][iRT]=p1;
            if (p1>0. && p1<1.) bs[m][iN][j][iRT]=bias_pfpt(p1,pta[t]); // f0e
            else                bs[m][iN][j][iRT]=(p1==0?-1:(p1==1?-2:-9));
         }

         if (r==0)
         {
            out<<setw(4)<<n<<setw(4)<<pda[j]<<setw(8)<<dToS_sci(pta[t],0)
               <<setw(8)<<setprecision(3)<<xt;
            for (m=0;m<nMTD;m++)
            {
               d1=ps[m][iN][j][iRT]; out<<setw(9)<<((d1>=0 && d1<=1.)?dToS_sci(d1,2):"IC");
            }
            for (m=0;m<nMTD;m++)
            {
               d1=bs[m][iN][j][iRT];
               if (d1>0) out<<setw(7)<<setprecision(3)<<d1;
               else out<<setw(7)<<(d1==-1?"0":(d1==-2?"1":"IC"));
            }
            out<<endl;
         }

         ou2<<setw(4)<<n<<setw(4)<<pda[j]<<setw(8)<<setprecision(3)<<xt<<setw(5)<<r+1
            <<setw(7)<<-log10(pta[t]);
         for (m=0;m<nMTD;m++) ou2<<setw(7)<<-log10(ps[m][iN][j][iRT]);
         for (m=0;m<nMTD;m++) ou2<<setw(7)<<bs[m][iN][j][iRT];
         for (i1=0;i1<nDI;i1++) ou2<<setw(7)<<w[i1]; ou2<<setw(7)<<Sum(w,nDI)<<endl;
      } // t
      if (r==nR-1) { i1=digit(Max(na,nN)); cout<<"TIME: pbbm() for n="<<setw(i1)<<n
         <<", pd="<<setw(2)<<pda[j]<<","<<setw(13-i1)<<""<<"time="; prtHms(cout,iPST); }
   } // j+r

   // print bias by n+pd across all stats xt
   out<<"\n"<<"(2) P-value bias for "<<nN<<" sample sizes and "<<nDI+1
      <<" population distributions\n"<<"(each row contains bias statistics for "<<nT
      <<" distinct p-values in "<<nR<<" replicated analyses):\n";
   out<<setw(4)<<"n"<<setw(4)<<"pd"<<setw(7)<<"n_ana"<<setw(13)<<"bias_pm"
      <<setw(13)<<"bias_gum"<<setw(13)<<"bias_gpd"<<setw(13)<<"bias_dad"<<endl;
   for (iN=0;iN<nN;iN++) for (j=0;j<nPD;j++)
   {
      out<<setw(4)<<na[iN]<<setw(4)<<pda[j]<<setw(7)<<nRT<<setprecision(3);
      for (m=0;m<nMTD;m++)
      {
         if (Min(bs[m][iN][j],nRT)<0)
         {
            d1=(double)CountLe(bs[m][iN][j],nRT,0.)*100/nRT;
            out<<setw(7)<<"IC"<<setw(5)<<setprecision(1)<<d1<<"%";
         }
         else out<<setw(7)<<setprecision(3)<<Mean(bs[m][iN][j],nRT)<<"+"
            <<setw(5)<<SD(bs[m][iN][j],nRT);
      }
      out<<endl;
   }
   for (iN=0;iN<nN;iN++)
   {
      out<<(iN==0?"\n":"")<<setw(4)<<na[iN]<<setw(11)<<nPD*nRT<<setprecision(3);
      for (m=0;m<nMTD;m++)
      {
         if (Min(bs[m][iN],nPD,nRT)<0)
         {
            d1=(double)CountLe(bs[m][iN],nPD,nRT,0.)*100/(nPD*nRT);
            out<<setw(7)<<"IC"<<setw(5)<<setprecision(1)<<d1<<"%";
         }
         else out<<setw(7)<<setprecision(3)<<Mean(bs[m][iN],nPD,nRT)<<"+"
            <<setw(5)<<SD(bs[m][iN],nPD,nRT);
      }
      out<<endl;
   }

   Drray1(na,10); Drray1(pda,10); Drray1(aa,nDI); Drray1(ba,nDI); Drray1(pta,nT);
   Drray2(xta,nDI,nT); Drray3(xt2a,nDI,nDI,nT); Drray1(da,nDI*nDI*nT);
   Drray1(xs,nM); Drray4(ps,nMTD,nN,nPD,nRT); Drray4(bs,nMTD,nN,nPD,nRT);
   Drray1(di,nDI); Drray1(w,nDI); out.close(); ou2.close();
   cout<<"NOTE: pbbm() done, results saved to "<<fn<<" and "<<fn2<<", time used ";
   prtHms(cout,iPST); cout<<endl;
}


// f2: get fn=(test size by mtd+pd+n+sl), fn2=(exp p and est p by mtd+pd)
void tsbm(double mp,double vp,char *pdg,char *ng,int nR, int iPST,char *fn)
{
   cout<<"TIME: tsbm() started:"<<setw(19)<<""<<"time="; prtHmsDttm(cout,iPST);
   int nPD,nN,nDI,nM,nMTD, j,k,l,pd,r,iN,n,m,t, *pda,*na; double x,p, *aa,*ba,*xs;
   int nSL,***nsg; double *sla,****t1e; double ***ps,*xo,*pexp;
   ofstream out,ou2; char fn2[100]=""; strcpy_newE(fn2,fn,"dat");
   openOF(out,fn); out<<setiosflags(ios::fixed|ios::showpoint);
   openOF(ou2,fn2); ou2<<setiosflags(ios::fixed|ios::showpoint);
   cout<<setiosflags(ios::fixed|ios::showpoint);

   Array1(pda,10); Array1(na,10); ipt_ia(pdg,pda,nPD); ipt_ia(ng,na,nN);
   nDI=5; nM=na[nN-1]+1; nMTD=2;
   Array1(aa,nDI); Array1(ba,nDI); Array1(xs,nM);

   // for tsbm=(tabulated results):
   nSL=(int)(log10((double)nR))-2;
   Array1(sla,nSL); Array3(nsg,nMTD,nN,nSL); Array4(t1e,nMTD,nPD,nN,nSL);
   for (t=0;t<nSL;t++) sla[t]=(t==0?0.05:pow(10.,-(t+1)));

   // for tsqq=(qqplot results):
   Array3(ps,nMTD,nN,nR); Array1(xo,nR); Array1(pexp,nR);
   ou2<<setw(2)<<"pd"<<setw(8)<<"rep"<<setw(10)<<"stat"<<setw(10)<<"pexp";
   for (m=0;m<nMTD;m++) for (iN=0;iN<nN;iN++) ou2<<setw(9)<<(m==0?"pgum_n":"pdad_n")<<iN+1;
   ou2<<endl;

   for (k=0;k<nDI;k++) getab_mv(k+1,mp,vp, aa[k],ba[k]); // f0b
   for (j=0;j<nPD;j++)
   {
      pd=pda[j];
      asgnArrayConst(nsg,nMTD,nN,nSL,0); // for tsbm
      for (r=0;r<nR;r++)
      {
         if (pd>=1 && pd<=5) simxs_pd(pd,aa[j],ba[j],nM, xs); // f0c
         else if (pd==11)
         {
            k=RdInt(0,nDI-1); l=RdInt(0,nDI-1);
            simxs_ws2d(k+1,aa[k],ba[k],l+1,aa[l],ba[l],nM,0.5,0.5, xs); // f0c_1
         }
         x=xs[nM-1];
         xo[r]=x; // for tsqq
         for (m=0;m<nMTD;m++) for (iN=0;iN<nN;iN++)
         {
            n=na[iN];
            if      (m==0) p=p_gum(xs,n,x); // pdad_f2a
            else if (m==1) p=p_dad(xs,n,x); // pdad_f1
            for (t=0;t<nSL;t++) nsg[m][iN][t]=nsg[m][iN][t]+(p<=sla[t]); // for tsbm
            ps[m][iN][r]=p; // for tsqq
         } // m,iN
         if (r==0 || r==nR-1 || (r+1)%(Max(nR,1000)/1000)==0)
         {
            cout<<"\rTIME: tsbm() for pd="<<setw(2)<<pd<<","<<setw(17)<<""<<"time=";
            prtHmsDttmNR(cout,iPST);
            cout<<setw(6)<<setprecision(1)<<((double)(j*nR+r+1)*100/(nR*nPD))<<"%"<<flush;
         }
      } // r

      for (m=0;m<nMTD;m++) for (iN=0;iN<nN;iN++) for (t=0;t<nSL;t++)
         t1e[m][j][iN][t]=(double)nsg[m][iN][t]/nR; // for tsbm

      // for tsqq: get+prt pExp+pEstByMtd
      pctgVec(xo,nR,pexp); for (r=0;r<nR;r++) pexp[r]=1.-pexp[r];
      for (r=0;r<nR;r++)
      {
         ou2<<setw(2)<<pda[j]<<setw(8)<<r+1<<setw(10)<<setprecision(4)<<xo[r]
            <<setw(10)<<setprecision(6)<<pexp[r];
         for (m=0;m<nMTD;m++) for (iN=0;iN<nN;iN++) ou2<<setw(10)<<ps[m][iN][r]; ou2<<endl;
      }
      cout<<endl;
   } // j

   // for tsbm: print T1E rate and ratio of T1E rate
   out<<setw(4)<<"n"<<setw(4)<<"pd";
   for (m=0;m<nMTD;m++) for (t=0;t<nSL;t++) out<<" t1e"
      <<(m==0?"gum_s":(m==1?"dad_s":"xxx_s"))<<t+1;
   for (m=0;m<nMTD;m++) for (t=0;t<nSL;t++) out<<" rat"
      <<(m==0?"gum_s":(m==1?"dad_s":"xxx_s"))<<t+1; out<<endl;
   for (iN=0;iN<nN;iN++) for (j=0;j<nPD;j++)
   {
      out<<setw(4)<<na[iN]<<setw(4)<<pda[j];
      out<<setprecision(6);
      for (m=0;m<nMTD;m++) for (t=0;t<nSL;t++) out<<setw(10)<<t1e[m][j][iN][t];
      out<<setprecision(3);
      for (m=0;m<nMTD;m++) for (t=0;t<nSL;t++) out<<setw(10)
         <<(t1e[m][j][iN][t]/sla[t]); out<<endl;
   }

   Drray1(pda,10); Drray1(na,10); Drray1(aa,nDI); Drray1(ba,nDI); Drray1(xs,nM);
   Drray1(sla,nSL); Drray3(nsg,nMTD,nN,nSL); Drray4(t1e,nMTD,nPD,nN,nSL); // for tsbm
   Drray3(ps,nMTD,nN,nR); Drray1(xo,nR); Drray1(pexp,nR); // for tsqq
   cout<<"NOTE: tsbm() done, results saved to "<<fn<<" and "<<fn2<<", time used ";
   prtHms(cout,iPST); cout<<endl;
}


// f3: Get sample size Number for acceptable bias limit By Method (nbm)
void nbm(int mtd,int pd,int ana,int nT,double mp,double vp,int nR,
   int iPST,char *fn,int o1a2,int aTot,int &aDn)
{
   int nLP,nSM,nST,nBL, t,k, *nmi,*nma,*ss;
   double a,b,qt,pt,xt,pDn, *bu,*xs; ofstream out;

   // should: a) nLP>=nSM, b) since nma_maxPmBu0.1~=0.8e6, so nST=1.6e6>2*nma, nSM=5e6>5*nma
   nT=6; nLP=10000000; nST=(mtd==1?1600000:30000); nSM=Min(nST*5,5000000); nBL=2;
   Array1(bu,nBL); Array1(nmi,nBL); Array1(nma,nBL); Array1(ss,nBL); Array1(xs,nLP);

   // 1: get sample size for bx=abs[(log10(px)-log10(pt))/log10(pt)] below a bias limit
   openApndOF(out,fn,o1a2); out<<setiosflags(ios::fixed|ios::showpoint);
   cout<<setiosflags(ios::fixed|ios::showpoint);
   if (nBL==2) { bu[0]=0.3; bu[1]=0.1; }
   if (mtd<1 || mtd>3) { cerr<<"\nERROR: nbm() w/ invalid mtd="<<mtd<<".\n\n"; exit(1); }
   if (o1a2==1)
   {
      out<<"mtd  pd  ana"<<setw(9)<<"stat"<<setw(8)<<"p_true";
      for (k=0;k<nBL;k++) out<<setw(8)<<"ss"<<k+1;
      out<<"  "; for (k=0;k<nBL;k++) out<<setw(5)<<"bu"<<k+1;
      for (k=0;k<nBL;k++) out<<setw(5)<<"nmi"<<k+1<<setw(8)<<"nma"<<k+1;
      out<<setw(11)<<"time_used"<<endl;

      cout<<"mtd pd ana"<<setw(7)<<"qt";
      for (k=0;k<nBL;k++) cout<<setw(8)<<"ss"<<k+1;
      for (k=0;k<nBL;k++) cout<<setw(5)<<"nmi"<<k+1<<setw(8)<<"nma"<<k+1;
      cout<<setw(11)<<"time_used"<<setw(8)<<"pct_dn"<<endl;
   }
   if      (pd>=1 && pd<=5)
   {
      getab_mv(pd,mp,vp, a,b); // f0b
      simxs_pd(pd,a,b,nSM, xs); // f0c
   }
   else if (pd==11)
   {
      getab_mv(1,mp,vp, a,b); // f0b
      simxs_ws2d(1,a,b,1,a,b,nLP,0.5,0.5, xs); // f0c_1
   }
   for (t=0;t<nT;t++)
   {
      qt=(t==0?-log10(0.5):(t==1?-log10(0.05):(double)t)); pt=pow(10.,-qt);
      if      (pd==1) xt=qchisq(1.-pt,a,b);
      else if (pd==2) xt=qgamma(1.-pt,a,b);
      else if (pd==3) xt=qlnorm(1.-pt,a,b);
      else if (pd==4) xt=qweibull(1.-pt,a,b);
      else if (pd==5) xt=qgumbel(1.-pt,a,b);
      else if (pd==11) xt=bisect_xt_minb(xs,nLP,pow(10.,-qt),1e-4,1e-4); // f0d
      for (k=0;k<nBL;k++)
      {
         nmi[k]=10;
         if      (mtd==1) nma[k]=Min(Max(5*(int)(pow(10.,t)),1000), nST);
         else if (mtd==2) nma[k]=Min((t<=1?500*(t+1):2500*(int)(pow(2.,t-2))), nST);
         bisect_nbm(mtd,xs,nSM,nST,nR,xt,qt,0.95,bu[k],1e-3, nmi[k],nma[k],ss[k]); // f3a
      }

      out<<setw(3)<<mtd<<setw(4)<<pd<<setw(5)<<ana<<setw(9)<<setprecision(3)<<xt;
      if (t<=1) out<<setw(8)<<setprecision(2)<<pt; else out<<setw(8)<<dToS_sci(pt,0);
      for (k=0;k<nBL;k++) out<<setw(9)<<ss[k];
      out<<"  "<<setprecision(2); for (k=0;k<nBL;k++) out<<setw(6)<<bu[k];
      for (k=0;k<nBL;k++) out<<setw(6)<<nmi[k]<<setw(9)<<nma[k];
      out<<"   "; prtHms(out,iPST);

      cout<<setw(3)<<mtd<<setw(3)<<pd<<setw(4)<<ana<<setw(7)<<setprecision(3)<<qt;
      for (k=0;k<nBL;k++) cout<<setw(9)<<ss[k];
      for (k=0;k<nBL;k++) cout<<setw(6)<<nmi[k]<<setw(9)<<nma[k];
      cout<<"   "; prtHmsNR(cout,iPST);
      aDn+=(int)pow(3.,t); pDn=(double)aDn*100./aTot;
      cout<<setw(7)<<setprecision(1)<<pDn<<"%"<<endl;
   } // t

   Drray1(bu,nBL); Drray1(nmi,nBL); Drray1(nma,nBL); Drray1(ss,nBL); Drray1(xs,nSM);
   out.close(); cout<<"NOTE: nbm() done, results "<<(o1a2==1?"saved":"appended")
      <<" to "<<fn<<".\n"<<endl;
}


// f3a: bisection mtd for sample size co f(xs[nSM],xt,qt)=bpctl-bu=0
void bisect_nbm(int mtd,double *xs,int nSM,int nST,int nR,double xt,double qt,
   double bpctg,double bu,double ef, int &nMi,int &nMa,int &pn)
{
   int it,pa,pb,i1; double fa,fb,fn,wa,wb, *xw,*bw; Array1(xw,nSM); Array1(bw,nR);

   pa=nMi; fa=getblmt(mtd,pa,xs,nSM,nR,xt,qt,bpctg, xw,bw)-bu; // f3a_1
   pb=nMa; fb=getblmt(mtd,pb,xs,nSM,nR,xt,qt,bpctg, xw,bw)-bu; // f3a_1
   while (fa<0 && nMi>2)
   {
      nMi=Max(nMi/2,5); pa=nMi; fa=getblmt(mtd,pa,xs,nSM,nR,xt,qt,bpctg, xw,bw)-bu;
      if (nMi==5 && fa<0) { pn=-5; goto end_bisect_nbm; }
   }
   if (fb>0)
   {
      if (pb>=nST) { pn=-pb-1; goto end_bisect_nbm; }
      else
      {
         while (fb>0)
         {
            nMa=Min(2*nMa,nST); pb=nMa; fb=getblmt(mtd,pb,xs,nSM,nR,xt,qt,bpctg, xw,bw)-bu;
            if (pb==nST && fb>0) { pn=-pb-1; goto end_bisect_nbm; }
         }
      }
   }

   if (fa<-ef || fb>ef) { cerr<<"\nERROR: bisect_nbm() has no root solution.\n\n"; exit(1); }
   if (fa<=ef || fb>=-ef) { if (fa<=-fb) { pn=pa; fn=fa; } else { pn=pb; fn=fb; } }
   it=0; pn=0;
   while (fa>0 && fb<0 && pb-pa>=2 && it<100)
   {
      pn=Round((pa+pb)/2); if (pb-pa>20) { i1=Min((pb-pa)/10,5); pn=pn+RdInt(-i1,i1); }
      fn=getblmt(mtd,pn,xs,nSM,nR,xt,qt,bpctg, xw,bw)-bu;
      if      (fn<-ef && fn>fb) { pb=pn; fb=fn; }
      else if (fn> ef && fn<fa) { pa=pn; fa=fn; }
      else
      {
         wa=fabs(fb)/(fabs(fa)+fabs(fb)); wb=fabs(fa)/(fabs(fa)+fabs(fb));
         pn=Round(wa*pa+wb*pb); fn=wa*fa+wb*fb;
         break;
      }
      it++;
   }

end_bisect_nbm: Drray1(xw,nSM); Drray1(bw,nR);
}


// f3a_1: get bpctl of bias boundary co n+xs[nSM]+nR+xt+qt+bpctg
double getblmt(int mtd,int n,double *xs,int nSM,int nR,double xt,double qt,
   double bpctg, double *xw,double *bw)
{
   int r; double d1,bb;

   for (r=0;r<nR;r++)
   {
      passArray_rdSubseg(xs,nSM,xw,n);
      if      (mtd==1) { d1=1.-Percentage(xw,n,xt); bw[r]=bias_pfpt(d1,pow(10.,-qt)); }
      else if (mtd==2) { d1=p_dad(xw,n,xt); bw[r]=bias_pfpt(d1,pow(10.,-qt)); }
   }
   bb=Percentile(bw,nR,bpctg);

   return bb;
}


// f3b: run multiple nbm() by mtd+pd+ana
void nbm_mrun(double mp,double vp,int nR,int nA, int iPST,char *fn)
{
   int mtd,j,pd,iA,t,mMi,mMa,nDI,nT, o1a2,aTot,aDn;

   mMi=1; mMa=2; nDI=5; nT=6; aDn=0;
   aTot=0; for (mtd=mMi;mtd<=mMa;mtd++) for (j=0;j<=nDI;j++) for (iA=0;iA<nA;iA++)
      for (t=0;t<nT;t++) aTot+=(int)pow(3.,t);
   for (mtd=mMi;mtd<=mMa;mtd++) for (j=0;j<=nDI;j++) for (iA=0;iA<nA;iA++)
   {
      pd=(j<nDI?j+1:11); o1a2=((mtd==mMi && pd==1 && iA==0)?1:2);
      nbm(mtd,pd,iA+1,nT,mp,vp,nR, iPST,fn,o1a2,aTot,aDn); // f3
   }
}
