/*
PROGRAM: f4_util.h
PURPOSE: All UTILITY AND CALCULATION FUNCTIONS
*/


#ifndef F4_UTIL_H
#define F4_UTIL_H


////////////////////////////////////////////////////////////
// [4] UTILITY FUNCTIONS
////////////////////////////////////////////////////////////

// [4.1] PRINT MESSAGE, STRING/MATRIX, TIME, SOFTWARE HEADER, ETC
// print error message and exit program
void prtMsg_exit(const char *msg)
{
   cerr<<"\nERROR: "<<msg<<".\n\n"; exit(1);
}
void prtMsg1i_exit(const char *msg,int i1)
{
   cerr<<"\nERROR: "<<msg<<" "<<i1<<".\n\n"; exit(1);
}

void errata1(int errCode)
{
   cerr<<"\nERROR";
   if      (errCode==1)  {cerr<<"(001): Integer value out of range.";}
   else if (errCode==2)  {cerr<<"(002): Double value out of range.";}
   else if (errCode==3)  {cerr<<"(003): Unable to allocate memory.";}
   else if (errCode==4)  {cerr<<"(004): Unable to deallocate memory.";}
   else if (errCode==5)  {cerr<<"(005): Unable to open input file.";}
   else if (errCode==6)  {cerr<<"(006): Unable to open output file."; }
   else if (errCode==7)  {cerr<<"(007): Variable is busy, cannot be modified.";}
   else if (errCode==8)  {cerr<<"(008): Error(s) found in inputted command file."; }
   else if (errCode==11) {cerr<<"(011): Mendelian inconsistancy.";}
   else if (errCode==999){cerr<<"(999): Program terminated due to unknown error(s).";}
   cerr<<"\n\n"; exit(1);
}
void errata2(int errCode, char *vNm)
{
   cerr<<"\nERROR";
   if      (errCode==1)  {cerr<<"(001): Integer value is out of range ("<<vNm<<").";}
   else if (errCode==2)  {cerr<<"(002): Double value is out of range ("<<vNm<<").";}
   else if (errCode==7)  {cerr<<"(007): Variable is busy, cannot be modified ("<<vNm<<").";}
   cerr<<"\n\n"; exit(1);
}
void prt_err_msg(char *c)
{
   cerr<<"\nERROR: "<<c<<"\n\n"; exit(1);
}
void errata4(char *c1, int i1, char *c2)
{
   cerr<<"\nERROR: "<<c1<<i1<<c2<<"\n\n"; exit(1);
}
void errata5(char *c1, int i1, char *c2, int i2, char *c3)
{
   cerr<<"\nERROR: "<<c1<<i1<<c2<<i2<<c3<<"\n\n"; exit(1);
}
void erratai(char *vNm, int v)
{
   cerr<<"\nERROR: Integer value is out of range ("<<vNm<<" = "<<v<<")."<<"\n\n"; exit(1);
}
void erratad(char *vNm, double v)
{
   cerr<<"\nERROR: Double value is out of range ("<<vNm<<" = "<<v<<").\n\n"; exit(1);
}

// print error message "cannot open input file" & exit program
void errOpenIF_exit(char *fn)
{
   if (strcmp(fn,"")==0) { cerr<<"\nERROR: Cannot find valid input file.\n\n"; exit(1); }
   else                  { cerr<<"\nERROR: Cannot open input file "<<fn<<".\n\n"; exit(1); }
}
void errOpenIF_exit(const char *fn)
{
   if (strcmp(fn,"")==0) { cerr<<"\nERROR: Cannot find valid input file.\n\n"; exit(1); }
   else                  { cerr<<"\nERROR: Cannot open input file "<<fn<<".\n\n"; exit(1); }
}

// print error message "cannot open output file" & exit program
void errOpenOF_exit(char *fn)
{
   if (strcmp(fn,"")==0) { cerr<<"\nERROR: Cannot find valid output file.\n\n"; exit(1); }
   else                  { cerr<<"\nERROR: Cannot open output file "<<fn<<".\n\n"; exit(1); }
}
void errOpenOF_exit(const char *fn)
{
   if (strcmp(fn,"")==0) { cerr<<"\nERROR: Cannot find valid output file.\n\n"; exit(1); }
   else                  { cerr<<"\nERROR: Cannot open output file "<<fn<<".\n\n"; exit(1); }
}


// get cSz=(longest line length in all lines ofIptFile) & ncDlm=nColDlm
void getCsz(char *fn, int &cSz,int &ncMi,int &ncMa)
{
   int i1; string sLn;
   ifstream din(fn,ios::in); if (!din) errOpenIF_exit(fn);

   ncMi=-1; ncMa=-1;
   while (getline(din,sLn))
   {
      cSz=Max(cSz,(int)sLn.length());
      i1=wordCount(sLn);
      if (ncMi==-1) ncMi=i1; else ncMi=Min(ncMi,i1);
      if (ncMa==-1) ncMa=i1; else ncMa=Max(ncMa,i1);
   }
}

// get rSz=(num of lines) of input file fn
void getRsz(char *fn, int &rSz,int &rEmpty)
{
   string sLn;
   ifstream din(fn,ios::in); if (!din) errOpenIF_exit(fn);

   rSz=0; rEmpty=0;
   while (getline(din,sLn)) { rSz++; if ((int)sLn.length()==0) rEmpty++; }
}

//  get dlm=(' ', ',', etc) from top 5 row ofIptFile
char getDlm(char *fn)
{
   char dlm; int nRow,i1,r1,r2,r,i2,j1,j2; string sLn; // nrDlm=nRowForDlm
   ifstream din(fn,ios::in); if (!din) errOpenIF_exit(fn);

   dlm='.'; getRsz(fn, nRow,i1); r1=0; r2=0;
   for (r=1;r<=nRow;r++)
   {
      getline(din,sLn);  i1=index(sLn,"  "); i2=index(sLn," "); j1=index(sLn,","); j2=index(sLn,"	");
      if      (i1>=1 && j1+j2==0)         { if (dlm=='.') { dlm='c'; r1=r; } else if (dlm!='c') r2=r; }
      else if (i2>=2 && j1+j2==0)         { if (dlm=='.') { dlm=' '; r1=r; } else if (dlm!=' ') r2=r; }
      else if (j1>=2 && (i2==0 || j1<i2)) { if (dlm=='.') { dlm=','; r1=r; } else if (dlm!=',') r2=r; }
      else if (j2>=2 && (i2==0 || j2<i2)) { if (dlm=='.') { dlm='	'; r1=r; } else if (dlm!='	') r2=r; }
   }
   if (dlm=='.' || r2>=1)
   {
      cerr<<"\nERROR: "<<fn<<": ";
      if      (dlm='.' && r2==0) cerr<<"Cannot determine file format from ";
      else if (r2>=1)            cerr<<"File format is different between rows "<<r1<<" and "<<r2;
      cerr<<".\n\n"; exit(1);
   }

   return dlm;
}

// prt file info of rowSz+colSz+fmt ofIptFile
void prtFileInfo(ostream &out,char *fn)
{
   int cSz,ncMi,ncMa,rSz,rEmpty; char dlm;

   getCsz(fn, cSz,ncMi,ncMa);
   getRsz(fn, rSz,rEmpty);
   dlm=getDlm(fn);

   out<<"***File information of "<<fn<<"***\n";
   out<<"   No. of lines               = "<<rSz<<"\n"
      <<"   Maximum line size          = "<<cSz<<"\n";
   out<<"   No. of delimited columns   = "<<ncMi; if (ncMa>ncMi) out<<" to "<<ncMa; out<<"\n";
   if (ncMi==ncMa && ncMi>=2)
   {
      out<<"   File format                = ";
      if      (dlm=='c') out<<"column-delimited values (cdv)";
      else if (dlm==' ') out<<"space-separated values (ssv)";
      else if (dlm==',') out<<"comma-separated values (csv)";
      else if (dlm=='	') out<<"tab-separated values (tsv)";
      else               out<<"unknown";
      out<<"\n";
   }

   if (rEmpty>=1) cout<<"WARNING: "<<rEmpty<<" line"<<(rEmpty>=2?"s are":" is")<<" empty.\n";
   out<<endl;
}


// set ios fixed show point & precision
void setFixedShowPoint(ostream &out,int prec)
{
   out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(prec);
}

// print double value in non-fixed format
void prtDbl(ostream &out,double d)
{
   int w1,dec;
   w1=digit(d); if (d<0) w1++;
   dec=decnnz(d,3);
   if (dec>=1) w1=w1+dec+1;
   if (dec>=1) out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(dec);
   if (dec==0) out<<(int)d; else if (dec>=1) out<<setw(w1)<<d;
}
void prtDblFW(ostream &out,double d,int w)
{
   int dec;

   dec=decnnz(d,3);
   if (dec>=1) out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(dec);
   if      (dec==0) out<<setw(w)<<(int)d;
   else if (dec>=1) out<<setw(w)<<d;
}
void prtDblFWLt(ostream &out,double d,int w)
{
   int dec,w1;

   dec=decnnz(d,3); w1=digit(d)+(dec>=1?dec+1:0);
   if (dec>=1) out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(dec);
   if      (dec==0) out<<(int)d<<setw(w-w1)<<"";
   else if (dec>=1) out<<d<<setw(w-w1)<<"";
}

// print characters in specified length
void prtCh(ostream &out, char c, int length)
{
   for (int i=0; i<length; i++) out<<c; out<<endl;
}

// print characters in specified length (w/o return)
void prtChNR(ostream &out, char c, int length)
{
   for (int i=0; i<length; i++) out<<c;
}


// print time for int sec in hh:mm:ss format
void prtHmsSecNR(ostream &out,int sec)
{
   int h,m,s; h = sec/3600; m = (sec/60)%60; s = sec%60;
   out<<(h<=9?"0":"")<<h<<":"<<(m<=9?"0":"")<<m<<":"<<(s<=9?"0":"")<<s<<flush;
}
void prtHmsSec(ostream &out,int sec)
{
   prtHmsSecNR(out,sec); out<<endl;
}

// get unit of millisecond in C++ co computer system (i.e., 1=Windows, 1000=Linux)
int getMilliSecUnitCoSystem()
{
   int i,j,jk,msUnit; clock_t msF,msL;

   msUnit=1000;
   for (j=0;j<10;j++)
   {
      msF=clock(); for (i=0;i<1000000;i++) jk=rand(); msL=clock();
      if (msF%1000!=0 || msL%1000!=0) msUnit=1;
   }
   return msUnit;
}

// print time from milliseconds msec in hh:mm:ss.sss format
void prtHmmMsecNR(ostream &out,int msec)
{
   int sec, res;
   sec=msec/1000; res=msec%1000;
   prtHmsSecNR(out,sec); out<<"."<<(res<=9?"00":(res<=99?"0":""))<<res;
}
void prtHmmMsec(ostream &out,int msec)
{
   prtHmmMsecNR(out,msec); out<<endl;
}

// print hh:mm:ss since int iBeg w/o & w/ return
void prtHmsNR(ostream &out,int iBeg)
{
   time_t tCur; time(&tCur);
   prtHmsSecNR(out, (int)tCur-iBeg);
}
void prtHms(ostream &out,int iBeg)
{
   prtHmsNR(out,iBeg); out<<endl;
}


// print date+time in standard 19-letter-format mm/dd/yyyy hh:mm:ss
void prtMdyHmsNR(ostream &out,time_t tIn)
{
   char c[20];
   strftime(c, 20, "%m/%d/%Y %H:%M:%S", localtime(&tIn)); 
   out<<c<<flush;
}
void prtMdyHms(ostream &out,time_t tIn)
{
   prtMdyHmsNR(out,tIn); out<<endl;
}
void prtMdyHmsNR(ostream &out) // current date+time wo rtn
{
   time_t tCur; time(&tCur); prtMdyHmsNR(out,tCur);
}
void prtMdyHms(ostream &out) // current date+time w/ rtn
{
   prtMdyHmsNR(out); out<<endl;
}


// print hh:mm:ss since iBeg and current date+time woRtn & wRtn
void prtHmsDttmNR(ostream &out,int iBeg)
{
   time_t tCur; time(&tCur);
   prtHmsNR(out,iBeg); out<<"   "; prtMdyHmsNR(out,tCur);
}
void prtHmsDttm(ostream &out,int iBeg)
{
   prtHmsDttmNR(out,iBeg); out<<endl;
}


// [4.2] SEARCH, COUNT, COPY/MODIFY/COMPARE & CONVERSION OF CHARACTER ARRAYS

// search string s1 for first occurrence of string s2
int index(string s1,const char *s2)
{
   int i, j, sz1, sz2, dif, col;
   sz1=(int)s1.length(); sz2=(int)strlen(s2); col=0;
   if (sz2<=sz1)
   {
      for (i=0; i<sz1-sz2+1; i++)
      {
         dif=0; for (j=0; j<sz2; j++) if (s1[i+j]!=s2[j]) dif=1;
         if (dif==0) { col=i+1; break; }
      }
   }
   return col;
}
int index(char *s1,char *s2)
{
   int i, j, sz1, sz2, dif, col;
   sz1=(int)strlen(s1); sz2=(int)strlen(s2); col=0;
   if (sz2<=sz1)
   {
      for (i=0; i<sz1-sz2+1; i++)
      {
         dif=0; for (j=0; j<sz2; j++) if (s1[i+j]!=s2[j]) dif=1;
         if (dif==0) { col=i+1; break; }
      }
   }
   return col;
}
int index(const char *s1,char *s2)
{
   char c[512]=""; strcpy(c,s1); return index(c,s2);
}
int index(char *s1,const char *s2)
{
   char c[512]=""; strcpy(c,s2); return index(s1,c);
}
int index(const char *s1,const char *s2)
{
   char c1[512]="",c2[512]=""; strcpy(c1,s1); strcpy(c2,s2); return index(c1,c2);
}

// search string s1 for nth occurrence of string s2
int index(char *s1,const char *s2,int n)
{
   int i, j, sz1, sz2, dif, col, nOcc;
   sz1=(int)strlen(s1); sz2=(int)strlen(s2); col=0; nOcc=0;
   if (sz2<=sz1)
   {
      for (i=0; i<sz1-sz2+1; i++)
      {
         dif=0; for (j=0; j<sz2; j++) if (s1[i+j]!=s2[j]) dif=1;
         if (dif==0) nOcc++;
         if (nOcc==n) { col=i+1; break; }
      }
   }
   return col;
}

// search upcase & compressed string s1 for the 1st occurrence of string s2
int indexUC(char *s1, char *s2)
{
   char s1UC[512]=""; strcpy_upcp(s1UC,s1);
   char s2UC[512]=""; strcpy_upcp(s2UC,s2);
   return index(s1UC, s2UC);
}

// search string s1 as a series of words for first occurrence of the whole string s2
int indexWd(char *s1,const char *s2)
{
   int col,nwd1,i,iwd1,iwd,dlmPre,sz1,j; char wd[20];

   col=0; nwd1=wordCount(s1);
   for (i=0;i<nwd1;i++)
   {
      ipt_nth_word(s1,i+1,wd);
      if (strcmp(wd,s2)==0)
      {
         iwd1=i+1;
         sz1=(int)strlen(s1); iwd=0; dlmPre=1;
         for (j=0;j<sz1;j++)
         {
            if      (dlmPre==1 && (s1[j]!=' ' && s1[j]!=',' && s1[j]!='	')) { iwd++; dlmPre=0; }
            else if (dlmPre==0 && (s1[j]==' ' || s1[j]==',' || s1[j]=='	')) dlmPre=1;
            if (iwd==iwd1) { col=j+1; break; }
         }
      }
   }
   return col;
}

// search string s1 for first occurrence of any word in string s2
int indexAnyWdS2(char *s1,const char *s2)
{
   int nwd,col,i,i1; char wd[20];

   nwd=wordCount(s2); col=0;
   for (i=0;i<nwd;i++)
   {
      ipt_nth_word(s2,i+1,wd);
      i1=index(s1,wd);
      if (col==0 && i1>=1) col=i1;
   }
   return col;
}

// search string s1 for first occurrence of any character in string s2
int indexc(char *s1,char *s2)
{
   int i, j, i1=(int)strlen(s1), i2=(int)strlen(s2), col=0;
   for (i=0;i<i1;i++) for (j=0;j<i2;j++)
   {
      if (s1[i]==s2[j]) {col=i+1; goto end_i;}
   }
end_i:;
   return col;
}
int indexc(char *s1,const char *s2)
{
   char c[512]=""; strcpy(c,s2); return index(s1,c);
}

// search string s1 for nth occurrence of any character in string s2
int indexc(char *s1,const char *s2, int n)
{
   int i,j,i1,i2,col,nOcc;
   i1=(int)strlen(s1); i2=(int)strlen(s2); col=0; nOcc=0;
   for (i=0;i<i1;i++)
   {
      for (j=0;j<i2;j++) { if (s1[i]==s2[j]) { nOcc++; break; } }
      if (nOcc==n) { col=i+1; break; }
   }
   return col;
}

// search string s1 for 1st occurrence of integer i1
int indexci(char *s1,int i1)
{
   int col; char c1[10]="";
   iToC(i1,c1); col=index(s1,c1);
   return col;
}

// search integer i1 for first occurrence of integer i2
int indexi(int i1,int i2)
{
   int col; char c1[10]="",c2[10]="";
   iToC(i1,c1); iToC(i2,c2); col=index(c1,c2);
   return col;
}

// count the number of occurrences of string s2 in string s1
int strCount(char *s1,const char *s2)
{
   int i,j,i1,i2,iJk, dif, nOcc;
   i1=(int)strlen(s1); i2=(int)strlen(s2); nOcc=0;
   if (i2<=i1)
   {
      for (i=0;i<i1-i2+1;i++)
      {
         dif=0; for (j=0;j<i2;j++) if (s1[i+j]!=s2[j]) dif=1;
         if (dif==0)
         {
            nOcc++; i=i+(i2-1); iJk=i;
            for (j=iJk;j<i1;j++) { if (s1[i+1]==' ') i++; else break; }
         }
      }
   }
   return nOcc;
}

// count the number of occurrences of any char of string s2 in string s1
int strCountAny(char *s1, char *s2)
{
   int i,j,i1,i2,nOcc,dblBlk;
   char c[3]="";
   i1=(int)strlen(s1); i2=(int)strlen(s2); nOcc=0;
   for (i=0;i<i1;i++)
   {
      dblBlk=0;
      if (i>=1)
      {
         strncpy(c,s1+i-1,2); c[2]='\0';
         if (strcmp(c,"  ")==0) dblBlk=1;
      }
      if (dblBlk==0)
      {
         for (j=0;j<i2;j++) { if (s1[i]==s2[j]) { nOcc++; break; } }
      }
   }
   return nOcc;
}

// count string s1 for any non-adjacent delimiters in string s2
int dlmCount(char *s1,const char *s2)
{
   int n,i1,i2,n1,n2,fd,dlmPre;

   n1=(int)strlen(s1); n2=(int)strlen(s2);
   n=0; dlmPre=0;
   for (i1=0;i1<n1;i1++)
   {
      fd=0; for (i2=0;i2<n2;i2++) { if (s1[i1]==s2[i2]) { fd=1; break; } }
      if      (i1==0 && fd==1) { n++; dlmPre=1; }
      else if (i1==0 && fd==0) { dlmPre=0; }
      else if (i1>=1 && fd==1 && dlmPre==0) { n++; dlmPre=1; }
      else if (i1>=1 && fd==0 && dlmPre==1) { dlmPre=0; }
   }
   return n;
}

// count number of words (i.e., consecutive non-space chars) in string sIn
int wordCount(string sIn)
{
   int i,sSz,nWd,dlmPre;

   sSz=(int)sIn.length(); nWd=0; dlmPre=1;
   for (i=0;i<sSz;i++)
   {
      if      (dlmPre==1 && (sIn[i]!=' ' && sIn[i]!=',' && sIn[i]!='	')) { nWd++; dlmPre=0; }
      else if (dlmPre==0 && (sIn[i]==' ' || sIn[i]==',' || sIn[i]=='	')) dlmPre=1;
   }
   return nWd;
}
int wordCount(char *sIn)
{
   int i,sSz,nWd,dlmPre;

   sSz=(int)strlen(sIn); nWd=0; dlmPre=1;
   for (i=0;i<sSz;i++)
   {
      if      (dlmPre==1 && (sIn[i]!=' ' && sIn[i]!=',' && sIn[i]!='	')) { nWd++; dlmPre=0; }
      else if (dlmPre==0 && (sIn[i]==' ' || sIn[i]==',' || sIn[i]=='	')) dlmPre=1;
   }
   return nWd;
}
int wordCount(const char *sIn)
{
   char c[512]=""; strcpy(c,sIn); return wordCount(c);
}


// string copy by removing any characters in string s1
void compressc(char *sOut,char *sIn,const char *s1)
{
   int i, j, i1=(int)strlen(sIn), i2=(int)strlen(s1), col=0, remove;
   for (i=0;i<i1;i++)
   {
      remove=0;
      for (j=0;j<i2;j++) { if (sIn[i]==s1[j]) { remove=1; break; } }
      if (remove==0) sOut[col++]=sIn[i];
   }
   sOut[col]='\0';
}

// string copy by removing begin & end blanks, and reducing multiple blanks to single ones
void strcpy_rmbk(char *sOut,char *sIn)
{
   int i,j,jFst,nPre,n,sSz;
   char s1[LS]="",s2[LS]="", cPre[2]="",c[2]="";

   // remove beginning, ending & multiple blanks
   sSz=(int)strlen(sIn); j=0; jFst=1;
   for (i=0;i<sSz;i++)
   {
      if      (jFst==1 && sIn[i]!=' ') { s1[j]=sIn[i]; jFst=0; j++; }
      else if (jFst==0 && i>=1)
      {
         cPre[0]=sIn[i-1]; cPre[1]='\0'; c[0]=sIn[i]; c[1]='\0';
         nPre=index(
"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ~!@#$%^&*()_+-={}[]|:;\"'<>,.?\\/",cPre);
         n=index(
"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ~!@#$%^&*()_+-={}[]|:;\"'<>,.?\\/",c);
         if      (nPre>=1 && n>=1) { s1[j]=sIn[i]; j++; }
         else if (nPre==0 && n>=1) { s1[j]=' '; j++; s1[j]=sIn[i]; j++; }
      }
   }
   s1[j]='\0';

   // remove space before "="
   sSz=(int)strlen(s1); j=0;
   for (i=0;i<sSz;i++)
   {
      if      (i<=sSz-2 && s1[i]==' ' && s1[i+1]=='=') { s2[j]=s1[i+1]; i++; j++; }
      else if (i<=sSz-2 && s1[i]==' ' && s1[i+1]==';') { s2[j]=s1[i+1]; i++; j++; }
      else                                             { s2[j]=s1[i]; j++; }
   }
   s2[j]='\0';
   strcpy(s1,s2);

   // remove space after "="
   sSz=(int)strlen(s1); j=0;
   for (i=0;i<sSz;i++)
   {
      if      (i<=sSz-2 && s1[i]=='=' && s1[i+1]==' ') { s2[j]=s1[i]; i++; j++; }
      else if (i<=sSz-2 && s1[i]==';' && s1[i+1]==' ') { s2[j]=s1[i]; i++; j++; }
      else                                             { s2[j]=s1[i]; j++; }
   }
   s2[j]='\0';
   strcpy(sOut,s2);
}

// string copy by removing begin+end+multiple blanks, and removing blanks franking operators "*|=;,"
void strcpy_rmbk_op(char *sOut,char *sIn)
{
   int i,j,sSz,opPre,opNxt; char s1[LS]="";

   strcpy_rmbk(s1,sIn);
   sSz=(int)strlen(s1); j=0;
   for (i=0;i<sSz;i++)
   {
      if      (i==0) { sOut[j]=s1[i]; j++; }
      else if (i>=1)
      {
         if (s1[i-1]=='*'||s1[i-1]=='|'||s1[i-1]=='='||s1[i-1]==';'||s1[i-1]==',') opPre=1;
         else opPre=0;
         if (i<=sSz-2 && s1[i+1]=='*'||s1[i+1]=='|'||s1[i+1]=='='||s1[i+1]==';'||s1[i+1]==',') opNxt=1;
         else opNxt=0;
         if (s1[i]!=' '||(opNxt==0 && opPre==0)) { sOut[j]=s1[i]; j++; }
      }
   }
   sOut[j]='\0';
}

// string copy with letters replaced by upper ones
void strcpy_up(char *sOut, char *sIn)
{
   int i,k=(int)strlen(sIn);
   for (i=0;i<k;i++)
   {
      if (isalpha(sIn[i])) sOut[i]=(char)toupper(sIn[i]); else sOut[i]=sIn[i];
   }
   sOut[k]='\0';
}
void strcpy_up(char *sOut,const char *sIn)
{
   char c[512]=""; strcpy(c,sIn); strcpy_up(sOut,c);
}

// string copy with letters replaced by upper & compressed ones
void strcpy_upcp(char *sOut, char *sIn)
{
   char sUp[512]="";
   strcpy_up(sUp,sIn);
   compressc(sOut,sUp," ");
}

// string copy with letters replaced by lower cases
void strcpy_lw(char *sOut, char *sIn)
{
   int i,k=(int)strlen(sIn);
   for (i=0; i<k; i++)
   {
      if (isalpha(sIn[i])) sOut[i]=(char)tolower(sIn[i]); else sOut[i]=sIn[i];
   }
   sOut[k]='\0';
}

// string copy after column n (n=0 for strcpy)
void strcpy_an(char *sOut, char *sIn, int n)
{
   int k=(int)strlen(sIn); strncpy(sOut,sIn+n,k-n); sOut[k-n]='\0';
}
void strcpy_an(char *sOut, const char *sIn, int n)
{
   int k=(int)strlen(sIn); strncpy(sOut,sIn+n,k-n); sOut[k-n]='\0';
}

// string copy before column n (n>=strlen+1 for strcpy)
void strcpy_bn(char *sOut, char *sIn, int n)
{
   int k1=(int)strlen(sIn), k2=Min(n-1, k1);
   if (k2>=1) { strncpy(sOut,sIn,k2); sOut[k2]='\0'; } else sOut[0]='\0';
}
void strcpy_bn(char *sOut, const char *sIn, int n)
{
   int k1=(int)strlen(sIn), k2=Min(n-1, k1);
   if (k2>=1) { strncpy(sOut,sIn,k2); sOut[k2]='\0'; } else sOut[0]='\0';
}

// string copy before column m and after column n (m <= n)
void strcpy_bman(char *sOut,char *sIn,int m,int n)
{
   int k1=(int)strlen(sIn), k2=Min(m-1, k1); char c[512]="";
   if (k2>=0 && m<=n)
   {
      strncpy(c,sIn,k2); strncpy(c+k2,sIn+n,k1-n); c[k1-n+k2]='\0'; strcpy(sOut,c);
   }
   else sOut[0]='\0';
}

// string copy after column m and before column n (m+2 <= n)
void strcpy_ambn(char *sOut,char *sIn,int m,int n)
{
   char *s1; int sSz; sSz=(int)strlen(sIn); Array1(s1,sSz+1);
   if (m<n-1) { strcpy_bn(s1,sIn,n); strcpy_an(sOut,s1,m); } else sOut[0]='\0';
   Drray1(s1,sSz+1);
}

// string copy from column m to column n (m <= n)
void strcpy_fmtn(char *sOut,char *sIn,int m,int n)
{
   strcpy_ambn(sOut,sIn,m-1,n+1);
}

// string copy after any delimiter in string s1
void strcpy_adlm(char *sOut,char *sIn,const char *s1)
{
   int i,j,k0,k1,b,e;

   if (indexc(sIn,s1)==0) sOut[0]='\0';
   else
   {
      k0=(int)strlen(sIn); k1=(int)strlen(s1); b=0; e=k0-1;
      if (k0>=1 && k1>=1)
      {
         for (i=0;i<k0;i++) for (j=0;j<k1;j++) { if (sIn[i]==s1[j]) { b=i+1; goto Nxt1; } }
      }
Nxt1:
      if (e>=b) { strncpy(sOut,sIn+b,e-b+1); sOut[e-b+1]='\0'; } else sOut[0]='\0';
   }
}

// string copy equal & after any delimiter in string s1
void strcpy_eadlm(char *sOut, char *sIn, char *s1)
{
   int n=indexc(sIn,s1);
   if (n==0) sOut[0]='\0'; else strcpy_an(sOut, sIn,n-1);
}

// string copy before any delimiter in string s1
void strcpy_bdlm(char *sOut, char *sIn, const char *s1)
{
   int szOut,szIn,i,j,k1,b,e;

   szOut=(int)strlen(sOut); szIn=(int)strlen(sIn); k1=(int)strlen(s1); b=0; e=szIn-1;
   if (szIn>=1 && k1>=1)
   {
      for (i=0;i<szIn;i++) for (j=0;j<k1;j++)
      {
         if (sIn[i]==s1[j]) { e=i-1; goto Nxt1; }
      }
   }
Nxt1:
   if (e>=b)
   {
      strncpy(sOut,sIn+b,e-b+1);
      for (i=e-b+1;i<Max(szIn,szOut);i++) sOut[i]='\0';
   }
   else { for (i=0;i<szOut;i++) sOut[i]='\0'; };
}

// string copy within and before the first delimiter in string s1
void strcpy_wbfd(char *sOut,char *sIn,const char *s1)
{
   int i,j,k0,k1,b,e;

   k0=(int)strlen(sIn); k1=(int)strlen(s1); b=0; e=k0-1;
   if (k0>=1 && k1>=1)
   {
      for (i=0;i<k0;i++) for (j=0;j<k1;j++) { if (sIn[i]==s1[j]) { e=i; goto Nxt1; } }
   }
Nxt1:
   if (e>=b) { strncpy(sOut,sIn+b,e-b+1); sOut[e-b+1]='\0'; } else sOut[0]='\0';
}

// string copy upcase and within/before the first delimiter in string s1
void strcpy_wbup(char *sOut, char *sIn,const char *s1)
{
   strcpy_up(sOut,sIn); strcpy_wbfd(sOut,sOut,s1);
}

// string copy after and before any delimiters in strings s1 and s2 respectively
void strcpy_adbd(char *sOut,char *sIn,const char *s1,const char *s2)
{
   strcpy_adlm(sOut,sIn,s1);
   strcpy_bdlm(sOut,sOut,s2);
}

// string copy of sIn after first occurrence of string s1
void strcpy_astr(char *sOut,char *sIn,char *s1)
{
   int k1, k2, k3;
   k1=(int)strlen(sIn); k2=(int)strlen(s1); k3=index(sIn,s1);
   if (k3>=1) { strncpy(sOut,sIn+(k2+k3-1),k1-k2-k3+1); sOut[k1-k2-k3+1]='\0'; }
   else sOut[0]='\0';
}
void strcpy_astr(char *sOut,char *sIn,const char *s1)
{
   char c[512]=""; strcpy(c,s1); strcpy_astr(sOut,sIn,c);
}

// string copy after upcase string s1
void strcpy_asup(char *sOut,char *sIn,const char *s1)
{
   int k1, k2, k3;
   char sInUp[512], s1Up[512];
   strcpy_up(sInUp,sIn);
   strcpy_up(s1Up,s1);
   k1=(int)strlen(sIn); k2=(int)strlen(s1); k3=index(sInUp,s1Up);
   if (k3>=1) { strncpy(sOut,sIn+(k2+k3-1),k1-k2-k3+1); sOut[k1-k2-k3+1]='\0'; }
   else sOut[0]='\0';
}

// string copy upcase & after upcase string s1
void strcpy_upas(char *sOut,char *sIn,const char *s1)
{
   char sInUp[512], s1Up[512]; strcpy_up(sInUp,sIn); strcpy_up(s1Up,s1);
   strcpy_astr(sOut, sInUp,s1Up);
}

// string copy before string s1
void strcpy_bstr(char *sOut,char *sIn,const char *s1)
{
   int k1; k1=index(sIn,s1);
   if (k1>=1)
   {
      strncpy(sOut,sIn,k1-1); sOut[k1-1]='\0';
   }
   else sOut[0]='\0';
}

// string copy before upcase string s1
void strcpy_bsup(char *sOut,char *sIn,const char *s1)
{
   int k1; char sInUp[512], s1Up[512];
   strcpy_up(sInUp,sIn); strcpy_up(s1Up,s1); k1=indexc(sInUp,s1Up);
   if (k1>=2) { strncpy(sOut,sIn,k1-1); sOut[k1-1]='\0'; }
   else sOut[0]='\0';
}

// string copy after string s1 and before string s2
void strcpy_a1b2(char *sOut,char *sIn,const char *s1,const char *s2)
{
   strcpy_astr(sOut,sIn,s1); strcpy_bstr(sOut,sOut,s2);
}

// string copy after upcase const string s1 and before upcase const string s2
void strcpy_abup(char *sOut,char *sIn,const char *s1,const char *s2)
{
   strcpy_asup(sOut,sIn,s1);
   strcpy_bsup(sOut,sOut,s2);
}

// string copy before upcase string s1 and after upcase string s2
void strcpy_baup(char *sOut,char *sIn,const char *s1,const char *s2)
{
   char sb1[512]="",sa2[512]=""; strcpy_bsup(sb1,sIn,s1); strcpy_asup(sa2,sIn,s2);
   strcpy_addC(sOut,sb1,sa2);
   strcpy_rmbk(sOut,sOut);
}

// string copy after string s1 and before any delimiter in string s2
void strcpy_asbd(char *sOut,char *sIn,const char *s1,const char *s2)
{
   int i0,i1,i2,k1,k2,n,jk,b,e;
   i0=(int)strlen(sIn); i1=(int)strlen(s1); i2=(int)strlen(s2); b=i0-1; e=0;
   k1=index(sIn,s1);
   if (k1>=1)
   {
      k2=0;
      for (n=1;n<4;n++) { jk=indexc(sIn,s2,n); if (jk>k1+i1) { k2=jk; break; } }
   }
   else k2=index(sIn,s2);
   if      (k1>=1 && k2>k1+i1) { b=k1+i1-1; e=k2-2; }
   else if (k1==0 && k2>=1) { b=0; e=k2-2; }
   else if (k1>=1 && k2==0) { b=k1+i1-1; e=i0-1; }
   if (e>=b) { strncpy(sOut,sIn+b,e-b+1); sOut[e-b+1]='\0'; } else sOut[0]='\0';
}



// string copy and add characters at the end (note: strcpy_addC(sIn,sIn,sAdd) is string_append)
void strcpy_addC(char *sOut,char *sIn,char *sExt)
{
   int i1=(int)strlen(sIn), i2=(int)strlen(sExt);
   strncpy(sOut, sIn, i1); sOut[i1]='\0';
   strncpy(sOut+i1, sExt, i2); sOut[i1+i2]='\0';
}
void strcpy_addC(char *sOut,char *sIn,const char *sExt)
{
   char c[512]=""; strcpy(c,sExt); strcpy_addC(sOut,sIn,c);
}
void strcpy_addC(char *sOut, const char *sIn, char *sExt)
{
   int i1=(int)strlen(sIn), i2=(int)strlen(sExt);
   strncpy(sOut, sIn, i1); sOut[i1]='\0';
   strncpy(sOut+i1, sExt, i2); sOut[i1+i2]='\0';
}

// string append w/ a const char & sLk in btw if sIO is non-empty
void strcpy_apd_sLk(string &sIO,string sAdd,const char *sLk)
{
   if ((int)sIO.length()>=1) sIO=sIO+sLk; sIO=sIO+sAdd;
}

// string copy and renew ext characters at the end
void strcpy_newE(char *sOut, char *sIn, const char *sExt)
{
   int i0=index(sIn,"."); char sPre[100]="";
   if      (i0>=1) { strcpy_bn(sPre,sIn,i0+1); strcpy_addC(sOut,sPre,sExt); }
   else if (i0==0) { strcpy_addC(sOut,sIn,"."); strcpy_addC(sOut,sOut,sExt); }
}

// string copy, delete last 4-char extension and add new characters
void strcpy_del4CAddC(char *sOut, char *sIn, char *sExt)
{
   int i1=(int)strlen(sIn), i2=(int)strlen(sExt);
   strncpy(sOut, sIn, i1-4);
   strncpy(sOut+i1-4, sExt, i2);
   sOut[i1-4+i2]='\0';
}

// string copy and renew 4-character extension from integer
void strcpy_newiExt(char *sOut, char *sIn, int iExt)
{
   if (iExt>999) errata1(1);
   char sExt[5]=".000";
   sExt[3]=(char)(48+iExt%10);
   sExt[2]=(char)(48+(iExt/10)%10);
   sExt[1]=(char)(48+(iExt/100)%10);
   strcpy_newE(sOut,sIn,sExt);
}

// string copy and add 4-character extension from interger
void strcpy_addiExt(char *sOut, char *sIn, int iExt)
{
   if (iExt>999) errata1(1);
   char sExt[5]=".000";
   sExt[3]=(char)(48+iExt%10);
   sExt[2]=(char)(48+(iExt/10)%10);
   sExt[1]=(char)(48+(iExt/100)%10);
   strcpy_addC(sOut,sIn,sExt);
}

// convert integer to charater array
void iToC(int iIn, char *sOut)
{
   int i, dgt=digit(iIn);
   for (i=0; i<dgt; i++) sOut[i]=(char)(48+(iIn/(int)pow(10.,dgt-1-i))%10);
   sOut[dgt]='\0';
}

// string copy and add integer at the end
void strcpy_addI(char *sOut, char *sIn, int i)
{
   char s[100]="";
   iToC(i, s); strcpy_addC(sOut,sIn,s);
}

// string copy and add integer with padded zeros
void strcpy_addI_pad0(char *sOut, char *sIn, int i,int i0Sz)
{
   char s0[100]="000000000000"; int iSz; iSz=digit(i);
   if (iSz<i0Sz)
   {
      s0[i0Sz-iSz]='\0'; strcpy_addC(sOut,sIn,s0); strcpy_addI(sOut,sOut,i);
   }
   else strcpy_addI(sOut,sIn,i);
}

// string copy and add integer & characters at the end
void strcpy_addIC(char *sOut, char *sIn, int i, char *s)
{
   char sJk[100]=""; strcpy_addI(sJk,sIn,i); strcpy_addC(sOut,sJk,s);
}

// string copy and add characters & integer at the end
void strcpy_addCI(char *sOut, char *sIn, char *s, int i)
{
   char sJk[100]=""; strcpy_addC(sJk,sIn,s); strcpy_addI(sOut,sJk,i);
}

// string copy, add characters & integer with padded zeros
void strcpy_addCI_pad0(char *sOut, char *sIn, char *s, int i, int i0Sz)
{
   char sJk[100]="";
   strcpy_addC(sJk,sIn,s); strcpy_addI_pad0(sOut,sJk,i,i0Sz);
}


// string copy, add integer before ".", and keep extension at the end
void strcpy_addI_keepE(char *sOut, char *sIn, int i)
{
   int i0=index(sIn,"."); char sPre[100]="", sExt[5]="";
   if (i0>=1)
   {
      strcpy_bn(sPre,sIn,i0); strcpy_an(sExt,sIn,i0-1); strcpy_addIC(sOut,sPre,i,sExt);
   }
   else if (i0==0) { strcpy_addI(sOut,sIn,i); }
}

// string copy, add chars before ".", and keep extension at the end
void strcpy_addC_keepE(char *sOut, char *sIn, char *sNew)
{
   int i0=index(sIn,"."); char sPre[100]="", sExt[5]="";
   if (i0>=1)
   {
      strcpy_bn(sPre,sIn,i0); strcpy_an(sExt,sIn,i0-1);
      strcpy_addC(sOut,sPre,sNew); strcpy_addC(sOut,sOut,sExt);
   }
   else if (i0==0)
   {
      strcpy_addC(sOut,sIn,sNew); strcpy_addC(sOut,sOut,"."); strcpy_addC(sOut,sOut,sExt);
   }
}
void strcpy_addC_keepE(char *sOut, const char *sIn, char *sNew)
{
   int i0=index(sIn,"."); char sPre[100]="", sExt[5]="";
   if (i0>=1)
   {
      strcpy_bn(sPre,sIn,i0); strcpy_an(sExt,sIn,i0-1);
      strcpy_addC(sOut,sPre,sNew); strcpy_addC(sOut,sOut,sExt);
   }
   else if (i0==0)
   {
      strcpy_addC(sOut,sIn,sNew); strcpy_addC(sOut,sOut,"."); strcpy_addC(sOut,sOut,sExt);
   }
}
void strcpy_addC_keepE(char *sOut, char *sIn, const char *sNew)
{
   int i0=index(sIn,"."); char sPre[100]="", sExt[5]="";
   if (i0>=1)
   {
      strcpy_bn(sPre,sIn,i0); strcpy_an(sExt,sIn,i0-1);
      strcpy_addC(sOut,sPre,sNew); strcpy_addC(sOut,sOut,sExt);
   }
   else if (i0==0)
   {
      strcpy_addC(sOut,sIn,sNew); strcpy_addC(sOut,sOut,"."); strcpy_addC(sOut,sOut,sExt);
   }
}


// string copy, add char before ".", and renew extension at the end
void strcpy_addC_newE(char *sOut,char *sIn,char *sAdd,char *sExt)
{
   int i0=index(sIn,"."); char sPre[100]="";
   if      (i0>=1) { strcpy_bn(sPre,sIn,i0); strcpy_addC(sOut,sPre,sAdd); }
   else if (i0==0) { strcpy_addC(sOut,sIn,sAdd); }
   strcpy_addC(sOut,sOut,"."); strcpy_addC(sOut,sOut,sExt);
}
void strcpy_addC_newE(char *sOut,const char *sIn,char *sAdd,char *sExt)
{
   int i0=index(sIn,"."); char sPre[100]="";
   if      (i0>=1) { strcpy_bn(sPre,sIn,i0); strcpy_addC(sOut,sPre,sAdd); }
   else if (i0==0) { strcpy_addC(sOut,sIn,sAdd); }
   strcpy_addC(sOut,sOut,"."); strcpy_addC(sOut,sOut,sExt);
}

// string update and assign indicator iNA (1=new file, 2=appd file)
void strupd_BCNA(char *sCur, char *sPre, char *sDef, int &iNA)
{
   int i1=strcmp(sCur,sPre), i2=strcmp(sCur,"");
   if      (i1!=0 && i2!=0) { iNA=1; strcpy(sPre,sCur); }
   else if (i1!=0 && i2==0) { iNA=2; strcpy(sCur,sPre); }
   else if (i1==0 && i2!=0) { iNA=2; }
   else { iNA=1; strcpy(sCur,sDef); strcpy(sPre,sDef); }
}


// strcmp ext for non-case-sensitive str match
int strcmp_extNCS(char *fn,const char *ext)
{
   int dif; char eFn[100]="",eUp[100]="";

   strcpy_upas(eFn,fn,"."); strcpy_up(eUp,ext);
   dif=strcmp(eFn,eUp);
   return dif;
}


// STRING CONVERSION TO/FROM OTHER DATA TYPES:
// convert double to string
string dToS(double x)
{
   ostringstream out; out<<x;
   string s1=out.str();
   return s1;
}

// convert double to string in scientific fmt
string dToS_sci(double x,int dec)
{
   ostringstream out; int neg,i1,i2;

   neg=0; if (x<0) { neg=1; x=-x; }
   i1=Ceil(-log10(x)); i2=Floor(log10(x));

   out<<setiosflags(ios::fixed|ios::showpoint);
   if (neg==1) out<<"-";
   if (x<10)
   {
      if (dec==0) out<<(int)(x*pow(10.,i1)); else out<<setprecision(dec)<<x*pow(10.,i1);
      if (i1>=1) out<<"E-"<<i1;
   }
   else       { out<<setprecision(dec)<<x/pow(10.,i2)<<"E"<<i2; }
   string s1=out.str();
   return s1;
}

// convert double to string in scientific or fix(for dmiS<=x<=1) fmt
string dToS_sciFix(double x,double xmiS,int decS,int decF)
{
   ostringstream out; int neg,i1,i2;

   neg=0; if (x<0) { neg=1; x=-x; }
   i1=Ceil(-log10(x)); i2=Floor(log10(x));

   out<<setiosflags(ios::fixed|ios::showpoint);
   if (neg==1) out<<"-";
   if (x>=-xmiS && x<=xmiS)
   {
      if (decS==0) out<<(int)(x*pow(10.,i1)); else out<<setprecision(decS)<<x*pow(10.,i1);
      if (i1>=1) out<<"E-"<<i1;
   }
   else if (x>=-1. && x<=1) out<<setprecision(decF)<<x;
   else { out<<setprecision(decS)<<x/pow(10.,i2)<<"E"<<i2; }
   string s1=out.str();
   return s1;
}


// [4.3] INPUT VARIABLE, ARRAY, LINE

// input after variable name for 1 integer
void ipt_avn(char *sIn,const char *vNm, int &vOut)
{
   char s1[512]=""; strcpy_asup(s1,sIn,vNm); strcpy_rmbk(s1,s1); strcpy_abup(s1,s1,"="," ;");
   if ((int)strlen(s1)>=1) vOut=atol(s1);
}

// input after variable name for 1 double
void ipt_avn(char *sIn,const char *vNm, double &vOut)
{
   char s1[512]="";
   strcpy_asup(s1,sIn,vNm);
   strcpy_rmbk(s1,s1);
   strcpy_adbd(s1,s1,"="," ;");
   if ((int)strlen(s1)>=1) vOut=atof(s1);
}

// input after variable name for 1st string
void ipt_avn_s1(char *sIn,const char *vNm, char *vOut)
{
   char c[512]="";

   strcpy_asup(c,sIn,vNm);
   strcpy_rmbk(c,c); // note: remove blanks so that input string contains no blanks
   if      (index(c,"\"")>=1) strcpy_adbd(c,c,"\"","\"");
   else if (index(c,"'")>=1 ) strcpy_adbd(c,c,"\'","\'");
   if (index(c,"=")==1) { strcpy_bdlm(c,c," /;"); strcpy_adlm(c,c,"="); compressc(c,c,"()"); }
   if ((int)strlen(c)>=1) strcpy(vOut,c);
}

// input after variable name for whole string
void ipt_avn_str(char *sIn,const char *vNm, char *vOut)
{
   char c[512]="";

   strcpy_asup(c,sIn,vNm);
   strcpy_rmbk(c,c); // note: remove blanks so that input string contains no blanks
   if      (index(c,"\"")>=1) strcpy_adbd(c,c,"\"","\"");
   else if (index(c,"'")>=1 ) strcpy_adbd(c,c,"\'","\'");
   //else if (index(c,"(")>=1 ) strcpy_adbd(c,c,"(",")");
   if (index(c,"=")==1) { strcpy_bdlm(c,c,"/;"); strcpy_adlm(c,c,"="); compressc(c,c,"()"); }
   if ((int)strlen(c)>=1) strcpy(vOut,c);
}

// input after variable name for string array (within "()" and seperated by ",")
void ipt_avn_sa(char *sIn,const char *vNm,char **vOut, int &vSz)
{
   int i0,i1,i2, i,j1; char c[512]="",v[512]="";

   strcpy_asup(c,sIn,vNm);
   strcpy_rmbk(c,c); // note: remove multi-blanks so that input string array contains no blanks
   i0=index(c,"="); i1=index(c,"("); i2=index(c,")");
   if (i1==0) i1=index(c,"=");
   if (i2==0) i2=index(c,";");
   if (i0>=1 && i1>=i0 && i1<i0+5 && i2>i1)
   {
      strcpy_ambn(c,c,i1,i2); vSz=dlmCount(c,",")+1;
      for (i=0;i<vSz;i++)
      {
         j1=index(c,",");
         if      (j1>=1) { strcpy_bn(v,c,j1); compressc(vOut[i],v," "); strcpy_an(c,c,j1); }
         else if (j1==0) strcpy_rmbk(vOut[i],c);
      }
   }
   else if (i0>=1) { strcpy_bdlm(c,c," /;"); strcpy_adlm(vOut[0],c,"="); vSz=1; }
}

// input after variable "title" for 1 string
void ipt_att(char *sIn, char *vOut)
{
   char c[512]=""; strcpy_asup(c,sIn,"TITLE");
   if      (index(c,"\"")>=1) strcpy_adbd(c,c,"\"","\"");
   else if (index(c,"'")>=1 ) strcpy_adbd(c,c,"\'","\'");
   else strcpy(c,"");
   strcpy(vOut,c);
}


// input integer before column n
void ipt_bn(char *sIn, int n, int &vOut)
{
   int i,i1,j1,j2,nmLst; char c[20]="",d[10]="";

   strcpy_bn(c,sIn,n);
   i1=(int)strlen(c); j1=0; j2=-1; nmLst=0;
   for (i=i1-1; i>=0; i--)
   {
      if (nmLst==0 && c[i]!=' ') { j2=i; nmLst=1; }
      if (nmLst==1 && c[i]==' ') { j1=i+1; break; }
   }
   if (j2>=0) strncpy(d,c+j1,j2-j1+1);
   if ((int)strlen(d)>=1) vOut=atol(d);
}

// input integer after column n
void ipt_an(char *sIn, int n, int &vOut)
{
   int i,i1,j1,j2,nmFst; char c[1000]="",d[100]="";

   strcpy_an(c,sIn,n);
   i1=(int)strlen(c); j1=-1; j2=i1-1; nmFst=0;
   for (i=0; i<i1; i++)
   {
      if (nmFst==0 && c[i]!=' ') { j1=i; nmFst=1; }
      if (nmFst==1 && c[i]==' ') { j2=i-1; break; }
   }
   if (j1>=0) strncpy(d,c+j1,j2-j1+1);
   if ((int)strlen(d)>=1) vOut=atol(d);
}

// input nth integer (n=x,...)
void ipt_nth_int(char *s,int n,int &vOut)
{
   int i,len,spaPre,nInt;

   len=(int)strlen(s);
   nInt=0; spaPre=1;
   for (i=0; i<len; i++)
   {
      if      (s[i]!=' ' && s[i]!=',' && s[i]!='*' && s[i]!='|' && spaPre==1)
      {
         nInt++; spaPre=0;
      }
      else if ((s[i]==' ' || s[i]==',' || s[i]=='*' || s[i]=='|') && spaPre==0) spaPre=1;
      if (nInt==n)
      {
         ipt_an(s,i, vOut); break;
      }
   }
}
void ipt_nth_int(const char *s,int n,int &vOut)
{
   char c[1000]=""; strcpy(c,s);
   ipt_nth_int(c,n,vOut);
}

// input nth double
void ipt_nth_dbl(char *s,int n,double &vOut)
{
   char s1[LS]="",c[20]=""; int nDlm,b,e,jk;
   int i1;

   strcpy_rmbk_op(s1,s);
   nDlm=dlmCount(s1," ,*|"); b=-1; e=-1;
   if (n>=nDlm+2)
   {
      cerr<<"\nERROR: Can not input the ";
      if (n>=1 && n<=3) cerr<<(n==1?"1st":(n==2?"2nd":"3rd")); else cerr<<n<<"th";
      cerr<<" value.\n\n"; exit(1);
   }
   if (n==1)
   {
      b=0; if (nDlm==0) e=(int)strlen(s1)-1; else e=indexc(s1," ,*|",n)-2;
   }
   else if (n>=2)
   {
      jk=indexc(s1," ,*|",n-1);
      i1=indexc(s1+(jk),".0123456789+-");
      b=jk+indexc(s1+(jk),".0123456789+-");
      if (n<=nDlm) e=indexc(s1," ,*|",n)-2; else e=(int)strlen(s1)-1;
   }
   if (e>=b && b>=0) { strncpy(c,s1+b,e-b+1); c[e-b+1]='\0'; } else c[0]='\0';
   if ((int)strlen(c)>=1) vOut=atof(c);
}
void ipt_nth_dbl(const char *s, int n, double &vOut)
{
   char c[LS]=""; strcpy(c,s); ipt_nth_dbl(c,n,vOut);
}

// input nth word (delimiter as space=' ' or comma=',' or tab='	')
void ipt_nth_word(char *sIn, int n, char *vOut)
{
   int i,j1,j2,sSz,spaPre,nOrd;

   j1=-1; j2=-1; sSz=(int)strlen(sIn); nOrd=0; spaPre=1;
   for (i=0;i<sSz;i++)
   {
      if      ((sIn[i]!=' ' && sIn[i]!=',' && sIn[i]!='	') && spaPre==1) { nOrd++; spaPre=0; }
      else if ((sIn[i]==' ' || sIn[i]==',' || sIn[i]=='	') && spaPre==0) spaPre=1;
      if (nOrd==n && j1==-1) j1=i;
      if (nOrd==n && j1>=0) { if (sIn[i]!=' ' && sIn[i]!=',' && sIn[i]!='	') j2=i; else break; }
   }
   if (j1>=0 && j2>=j1)
   {
      if      (j1==0 && j2<sSz-1)  strcpy_bn(vOut,sIn,j2+2);
      else if (j1==0 && j2==sSz-1) strcpy(vOut,sIn);
      else if (j1>=1 && j2<sSz-1)  strcpy_ambn(vOut,sIn,j1,j2+2);
      else if (j1>=1 && j2==sSz-1) strcpy_an(vOut,sIn,j1);
   }
   else vOut[0]='\0';
}
void ipt_nth_word(const char *sIn,int n,char *vOut)
{
   int sSz; char *s1;

   sSz=(int)strlen(sIn); Array1(s1,sSz+1);
   strcpy(s1,sIn); ipt_nth_word(s1,n,vOut);
   Drray1(s1,sSz+1);
}
void ipt_nth_word(string sIn, int n, char *vOut)
{
   int sSz; char *s1;

   sSz=(int)sIn.length(); Array1(s1,sSz+1);
   strcpy(s1,sIn.c_str()); ipt_nth_word(s1,n,vOut);
   Drray1(s1,sSz+1);
}

// input nth word & upcase
void ipt_nth_wdup(char *sIn,int n, char *vOut)
{
   ipt_nth_word(sIn,n,vOut); strcpy_up(vOut,vOut);
}
void ipt_nth_wdup(string sIn,int n, char *vOut)
{
   ipt_nth_word(sIn,n,vOut); strcpy_up(vOut,vOut);
}


// input anyIntDbl array size (w/o array)
int ipt_sz(char *sIn)
{
   int vSz,n,nDlm,cnt,j; double d1,d2,d3; char s[512]="";
   strcpy_rmbk_op(s,sIn);

   nDlm=dlmCount(s," *|"); vSz=0;
   if (nDlm==0) { vSz++; }
   else if (nDlm>=1)
   {
      for (n=1;n<=nDlm;n++)
      {
         ipt_nth_dbl(s,n,d1); j=indexc(s," *|",n)-1;
         if      (s[j]==' ') { vSz++; }
         else if (s[j]=='*')
         {
            ipt_nth_int(s,n+1,cnt); n++;
            vSz+=cnt;
         }
         else if (s[j]=='|')
         {
            ipt_nth_dbl(s,n+1,d2); ipt_nth_dbl(s,n+2,d3); n+=2;
            cnt=(int)Round((d2-d1)/d3,1.);
            vSz+=(cnt+1);
         }
         if (n==nDlm) { vSz++; }
      }
   }
   return vSz;
}

// input integer array ('*' = nRepValues, '|' = ) & its size
// TBD: to be revised to handle "(1 2 3)*2" as "1 2 3 1 2 3".
void ipt_ia(char *sIn,int *v,int &vSz)
{
   int nDlm,n,cnt,j,k; int i1,i2,i3; char s[512]="";
   strcpy_rmbk_op(s,sIn); strcpy_lw(s,s);

   nDlm=dlmCount(s," ,*|"); vSz=0;
   if (nDlm>=1)
   {
      for (n=1;n<=nDlm;n++)
      {
         ipt_nth_int(s,n,i1); j=indexc(s," ,*|",n)-1;
         if      (s[j]==' ' || s[j]==',') { v[vSz]=(int)i1; vSz++; }
         else if (s[j]=='*') // i1*cnt = multiple integers
         {
             ipt_nth_int(s,n+1,cnt); n++;
             for (k=vSz;k<vSz+cnt;k++) v[k]=i1; vSz+=cnt;
         }
         else if (s[j]=='|') // i1|i2|i3 = iBeg|iEnd|iBy
         {
            ipt_nth_int(s,n+1,i2); ipt_nth_int(s,n+2,i3); n+=2;
            cnt=(int)((i2-i1)/i3);
            for (k=0;k<=cnt;k++) v[vSz+k]=i1+k*i3; vSz+=(cnt+1);
         }
         if (n==nDlm) { ipt_nth_int(s,n+1,i2); v[vSz]=i2; vSz++; }
      }
   }
   else if (nDlm==0)
   {
      ipt_nth_int(s,1,i1); v[0]=i1; vSz++;
   }
}
void ipt_ia(const char *sIn,int *v,int &vSz)
{
   char c[1000]=""; strcpy(c, sIn);
   ipt_ia(c,v,vSz);
}

// input double array
void ipt_da(char *sIn,double *v,int &vSz)
{
   int n,nDlm,cnt,j,k; double d1,d2,d3,d3Exact; char s[LS]="";
   strcpy_rmbk_op(s,sIn);

   nDlm=dlmCount(s," *|"); vSz=0;
   if (nDlm==0)
   {
      v[0]=atof(s); vSz++;
   }
   else if (nDlm>=1)
   {
      for (n=1;n<=nDlm;n++)
      {
         ipt_nth_dbl(s,n,d1);
         j=indexc(s," *|",n)-1;
         if      (s[j]==' ') { v[vSz]=d1; vSz++; }
         else if (s[j]=='*')
         {
            ipt_nth_int(s,n+1,cnt); n++;
            for (k=vSz;k<vSz+cnt;k++) v[k]=d1; vSz+=cnt;
         }
         else if (s[j]=='|')
         {
            ipt_nth_dbl(s,n+1,d2); ipt_nth_dbl(s,n+2,d3); n+=2;
            cnt=(int)Round((d2-d1)/d3,1.); d3Exact=(d2-d1)/(double)cnt;
            for (k=0;k<=cnt;k++) v[vSz+k]=(d1+k*d3Exact); vSz+=(cnt+1);
         }
         if (n==nDlm) { ipt_nth_dbl(s,n+1,d2); v[vSz]=d2; vSz++; }
      }
   }
}
void ipt_da(const char *sIn,double *v,int &vSz)
{
   char c[1000]=""; strcpy(c, sIn);
   ipt_da(c,v,vSz);
}

// input integer array after variable name
void ipt_ia_avn(char *sIn,char *vNm,int *v,int &vSz)
{
   char s[512]=""; strcpy_asup(s,sIn,vNm);
   strcpy_adbd(s,s,"=","abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ;");
   ipt_ia(s,v,vSz);
}

// input double array after variable name
void ipt_da_avn(char *sIn, char *vNm,double *v,int &vSz)
{
   char s1[512]=""; strcpy_asup(s1,sIn,vNm);
   strcpy_adbd(s1,s1,"=","abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ;");
   ipt_da(s1,v,vSz);
}

// input integer array from 1 or multiple lines
void ipt_ia_ml(istream &din,char *vNm,int *v,int sIn)
{
   const int Lsz=512; char Ln[Lsz]; int sOut,nNonEptLn, sDlm,nv,sJk;
   int *vJk; Array1(vJk,sIn);
   sOut=0; nNonEptLn=0;
   while (sOut<sIn)
   {
      din.getline(Ln,Lsz,'\n');
      sDlm=index(Ln,"//"); if (sDlm>=1) strcpy_bstr(Ln,Ln,"//");
      if (indexc(Ln,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")>=1)
         {cerr<<"\nERROR: No. of "<<vNm<<" integer elements is "<<sOut<<", not "<<sIn<<".\n\n"; exit(1);}
      nv=wordCount(Ln);
      if (nv>=1)
      {
         nNonEptLn++;
         ipt_ia(Ln,vJk,sJk);
         passArray_OP(vJk,v,sOut,sOut+sJk-1);
         sOut=sOut+sJk;
      }
      if (nNonEptLn==1 && sOut==1) break;
   }
   if (sOut!=1 && sOut!=sIn)
      { cerr<<"\nERROR: No. of "<<vNm<<" integer elements is "<<sOut<<", not "<<sIn<<".\n\n"; exit(1); }
   Drray1(vJk,sIn);
}

// input double array from 1 or multiple lines
void ipt_da_ml(istream &din,char *vNm,double *v,int sIn)
{
   const int Lsz=512; char Ln[Lsz]; int sOut,nNonEptLn, sDlm,nv,sJk;
   double *vJk; Array1(vJk,sIn);
   sOut=0; nNonEptLn=0;
   while (sOut<sIn)
   {
      din.getline(Ln,Lsz,'\n');
      sDlm=index(Ln,"//"); if (sDlm>=1) strcpy_bstr(Ln,Ln,"//");
      if (indexc(Ln,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")>=1)
         { cerr<<"\nERROR: No. of "<<vNm<<" double elements is "<<sOut<<", not "<<sIn<<".\n\n"; exit(1); }
      nv=wordCount(Ln);
      if (nv>=1)
      {
         nNonEptLn++;
         ipt_da(Ln,vJk,sJk);
         passArray_OP(vJk,v,sOut,sOut+sJk-1);
         sOut=sOut+sJk;
      }
      if (nNonEptLn==1 && sOut==1) break;
   }
   if (sOut!=1 && sOut!=sIn)
   { cerr<<"\nERROR: No. of "<<vNm<<" double elements is "<<sOut<<", not "<<sIn<<".\n\n"; exit(1); }
   Drray1(vJk,sIn);
}


// input akw=(after key word from program lines) for int/dbl array size & check its max value
void ipt_akw_szPli_chkSz(char **Pgm,int Npl,const char *kw,int szMa, int &vSz,int &Pli)
{
/*
   int *ia; Array1(ia,szMa);
   ipt_akw_iaSzPli(Pgm,Npl,kw, ia,vSz,Pli);
   Drray1(ia,szMa);
*/
   int p; char c[512]="",kwP[512]=""; vSz=0; Pli=-1;
   for (p=0;p<Npl;p++)
   {
      strcpy_up(c,Pgm[p]); strcpy_bdlm(kwP,c,"=");
      if (strcmp(kwP,kw)==0)
      {
         strcpy_abup(c,c,"=",";");
         if ((int)strlen(c)>=1)
         {
            vSz=ipt_sz(c);
            if (vSz>szMa) { cerr<<"\nERROR: Array size is out of range [1, "<<szMa<<"]: "
               <<Pgm[Pli]<<"\n\n"; exit(1); }
            Pli=p;
         }
      }
   }
}

// input akw for 1 integer+pgmLnIdx
void ipt_akw_i1Pli(char **Pgm,int Npl,const char *kw, int &v1,int &Pli)
{
   int p; char c[512]="",kwP[512]=""; v1=-1; Pli=-1;
   for (p=0;p<Npl;p++)
   {
      strcpy_up(c,Pgm[p]); strcpy_bdlm(kwP,c,"=");
      if (strcmp(kwP,kw)==0)
      {
         strcpy_abup(c,c,"=",";");
         if ((int)strlen(c)>=1) { v1=atol(c); Pli=p; }
      }
   }
}

// input akw for 1int & prtErrForNotFd
void ipt_akw_i1Pli_errNF(char **Pgm,int Npl,const char *kw, int &v1,int &Pli)
{
   ipt_akw_i1Pli(Pgm,Npl,kw, v1,Pli);
   if (Pli==-1) prt_akw_err(Pgm,kw,Pli);
}

void ipt_akw_i1_Nu(char **Pgm,int Npl,const char *kw,int &v1,int v1Nu) // ipt akw 1int w/ null val forNotFd
{
   int Pli; ipt_akw_i1Pli(Pgm,Npl,kw, v1,Pli);
   if (Pli==-1) v1=v1Nu;
}


// input akw from program lines for v[]=integer array, vSz=array size & Pli=pgmLnIdx
void ipt_akw_iaSzPli(char **Pgm,int Npl,const char *kw,int *v,int &vSz,int &Pli)
{
   int p,i1,i2; char c[512]="",kwP[512]=""; vSz=0; Pli=-1;
   for (p=0;p<Npl;p++)
   {
      strcpy_up(c,Pgm[p]); strcpy_bdlm(kwP,c,"=");
      if (strcmp(kwP,kw)==0)
      {
         i1=index(c,"="); i2=index(c,"(");
         if      (i1>=2 && i2==0) strcpy_abup(c,c,"=",";");
         else if (i1>=2 && i2>=3) strcpy_abup(c,c,"(",")");
         else { cerr<<"\nERROR: Invalid statement: "<<Pgm[Pli]<<"\n\n"; exit(1); }
         if ((int)strlen(c)>=1) { ipt_ia(c,v,vSz); Pli=p; }
      }
   }
}
void ipt_akw_iaSz(char **Pgm,int Npl,const char *kw,int *v,int &vSz)
{
   int Pli; ipt_akw_iaSzPli(Pgm,Npl,kw, v,vSz,Pli);
}

// input after key word from program lines for integer array w/ check of ArraySz+Rng
void ipt_akw_ia_chkSz(char **Pgm,int Npl,const char *kw,int *v,int vSzIn)
{
   int Pli,vSz; ipt_akw_iaSzPli(Pgm,Npl,kw, v,vSz,Pli);
   if      (Pli==-1) { cerr<<"\nERROR: Keyword "<<kw<<" is not found.\n\n"; exit(1); }
   else if (Pli>=0 && vSz!=vSzIn) {
      cerr<<"\nERROR: Invalid array size: "<<Pgm[Pli]<<"\n\n"; exit(1); }
}
void ipt_akw_ia_chkSzRg(char **Pgm,int Npl,const char *kw,int *v,int vSzIn,int vMi,int vMa)
{
   int Pli,vSz; ipt_akw_iaSzPli(Pgm,Npl,kw, v,vSz,Pli);
   if      (Pli==-1) { cerr<<"\nERROR: Keyword "<<kw<<" is not found.\n\n"; exit(1); }
   else if (Pli>=0 && vSz!=vSzIn) {
      cerr<<"\nERROR: Invalid array size: "<<Pgm[Pli]<<"\n\n"; exit(1); }
   else if (Pli>=0 && (Min(v,vSzIn)<vMi || Max(v,vSzIn)>vMa)) {
      cerr<<"\nERROR: Out of range values: "<<Pgm[Pli]<<"\n\n"; exit(1); }
}

// input after key word from program lines for 2 integers
void ipt_akw_i2(char **Pgm,int Npl,const char *kw,int &v1,int &v2)
{
   int *v; Array1(v,2);
   ipt_akw_ia_chkSz(Pgm,Npl,kw, v,2); v1=v[0]; v2=v[1];
   Drray1(v,2);
}
void ipt_akw_i4(char **Pgm,int Npl,const char *kw,int &v1,int &v2,int &v3,int &v4)
{
   int *v; Array1(v,4);
   ipt_akw_ia_chkSz(Pgm,Npl,kw, v,4); v1=v[0]; v2=v[1]; v3=v[2]; v4=v[3];
   Drray1(v,4);
}


// input after key word from program lines for 1 double+pgmLnIdx
void ipt_akw_d1Pli(char **Pgm,int Npl,const char *kw, double &v1,int &Pli)
{
   int p,i1; char c[512]="",w1[512]=""; Pli=-1;
   for (p=0;p<Npl;p++)
   {
      strcpy_up(c,Pgm[p]);
      strcpy_bdlm(w1,c,"=");
      if (strcmp(w1,kw)==0)
      {
         i1=index(c,"=");
         strcpy_abup(c,c,"=",";");
         if ((int)strlen(c)>=1) { v1=atof(c); Pli=p; }
      }
   }
}
void ipt_akw_d1_Nu(char **Pgm,int Npl,const char *kw, double &v1,double v1Nu)
{
   int Pli; ipt_akw_d1Pli(Pgm,Npl,kw, v1,Pli);
   if (Pli==-1) v1=v1Nu;
}

// input after key word from program lines for 1 double w/ check range
void ipt_akw_d1_chkRg(char **Pgm,int Npl,const char *kw,double &v1,double v1Mi,double v1Ma)
{
   int Pli; ipt_akw_d1Pli(Pgm,Npl,kw, v1,Pli);
   if      (Pli==-1) { cerr<<"\nERROR: Keyword "<<kw<<" is not found.\n\n"; exit(1); }
   else if (Pli>=0 && (v1<v1Mi || v1>v1Ma)) {
      cerr<<"\nERROR: Out of range value: "<<Pgm[Pli]<<"\n\n"; exit(1); }
}

// input after keyword+pliPre for v[]=double array, vSz=array size, Pli=pgmLnIdx
void ipt_akwPli_daSzPli(char **Pgm,int Npl,const char *kw,int pliP, double *v,int &vSz,int &Pli)
{
   int p,p0,i1; char c[512]="",w1[512]=""; vSz=0; Pli=-1;
   p0=Max(pliP+1,0);
   for (p=p0;p<Npl;p++)
   {
      strcpy_up(c,Pgm[p]); strcpy_bdlm(w1,c,"=");
      if (strcmp(w1,kw)==0)
      {
         i1=index(c,"="); strcpy_abup(c,c,"=",";");
         if ((int)strlen(c)>=1) { ipt_da(c,v,vSz); Pli=p; break; }
      }
   }
}

// input after kwyword for v[]=double array, vSz=array size, Pli=pgmLnIdx
void ipt_akw_daSzPli(char **Pgm,int Npl,const char *kw, double *v,int &vSz,int &Pli)
{
   int p,i1,i2; char c[512]="",w1[512]=""; vSz=0; Pli=-1;
   for (p=0;p<Npl;p++)
   {
      strcpy_up(c,Pgm[p]); strcpy_bdlm(w1,c,"=");
      if (strcmp(w1,kw)==0)
      {
         i1=index(c,"="); i2=index(c,"(");
         if      (i1>=2 && i2==0) strcpy_abup(c,c,"=",";");
         else if (i1>=2 && i2>=3) strcpy_abup(c,c,"(",")");
         else { cerr<<"\nERROR: Invalid statement: "<<Pgm[Pli]<<"\n\n"; exit(1); }
         if ((int)strlen(c)>=1) { ipt_da(c,v,vSz); Pli=p; }
      }
   }
}
void ipt_akw_daSz(char **Pgm,int Npl,const char *kw,double *v,int &vSz)
{
   int Pli; ipt_akw_daSzPli(Pgm,Npl,kw, v,vSz,Pli);
}
void ipt_akw_daSz_chkSz(char **Pgm,int Npl,const char *kw,double *v,int &vSz,int vSzMi,int vSzMa)
{
   int Pli; ipt_akw_daSzPli(Pgm,Npl,kw, v,vSz,Pli);
   if      (Pli==-1) { cerr<<"\nERROR: Keyword "<<kw<<" is not found.\n\n"; exit(1); }
   else if (Pli>=0)
   {
      if      (vSz==0) { cerr<<"\nERROR: Invalid statement:\n"<<Pgm[Pli]<<"\n\n"; exit(1); }
      else if (vSz<vSzMi || vSz>vSzMa) { cerr<<"\nERROR: Invalid array size: "<<Pgm[Pli]<<"\n\n";
         exit(1); }
   }
}

// input after key word from program lines for 2 doubles
void ipt_akw_d2(char **Pgm,int Npl,const char *kw,double &v1,double &v2)
{
   int vSz; double *v; Array1(v,2);
   ipt_akw_daSz_chkSz(Pgm,Npl,kw, v,vSz,2,2); v1=v[0]; v2=v[1];
   Drray1(v,2);
}
void ipt_akw_d4(char **Pgm,int Npl,const char *kw,double &v1,double &v2,double &v3,double &v4)
{
   int vSz; double *v; Array1(v,4);
   ipt_akw_daSz_chkSz(Pgm,Npl,kw, v,vSz,4,4); v1=v[0]; v2=v[1]; v3=v[2]; v4=v[3];
   Drray1(v,4);
}


// iptAftKeyWd fromPgmLn for 1+ int & 1+ dbl w/ null values
void ipt_akw_i1d1_Nu(char **Pgm,int Npl,const char *kw, int &i1,double &d1, int i1Nu,double d1Nu)
{
   char str1[100],str2[100]; int Pli;

   ipt_akw_s2(Pgm,Npl,kw,str1,str2,Pli);
   if (Pli>=0 && Pli<20 && strlen(str1)>=1 && strlen(str2)>=1)
   {
      i1=atol(str1); d1=atof(str2);
   }
   else
   {
      i1=i1Nu; d1=d1Nu;
   }
}
void ipt_akw_i1d2_Nu(char **Pgm,int Npl,const char *kw, int &i1,double &d1,double &d2,
                     int i1Nu,double d1Nu,double d2Nu)
{
   char str1[100],str2[100],str3[100]; int Pli;

   ipt_akw_s3(Pgm,Npl,kw,str1,str2,str3,Pli);
   if (Pli>=0 && Pli<20 && strlen(str1)>=1 && strlen(str2)>=1 && strlen(str3)>=1)
   {
      i1=atol(str1); d1=atof(str2); d2=atof(str3);
   }
   else
   {
      i1=i1Nu; d1=d1Nu; d2=d2Nu;
   }
}

// iptAftKeyWd fromPgmLn for 1+ dbl & 1+ int w/ null values
void ipt_akw_d1i1_Nu(char **Pgm,int Npl,const char *kw, double &d1,int &i1, double d1Nu,int i1Nu)
{
   char str1[100],str2[100]; int Pli;

   ipt_akw_s2(Pgm,Npl,kw,str1,str2,Pli);
   if (Pli>=0 && Pli<20 && strlen(str1)>=1 && strlen(str2)>=1)
   {
      d1=atof(str1); i1=atol(str2);
   }
   else
   {
      d1=d1Nu; i1=i1Nu;
   }
}


// print error msg after keyword for erroneous keyword statement
void prt_akw_err(char **Pgm,const char *kw,int Pli)
{
   if      (Pli< 0) { cerr<<"\nERROR: Keyword "<<kw<<" is not found.\n\n"; exit(1); }
   else if (Pli>=0) { cerr<<"\nERROR: Invalid statement:\n"<<Pgm[Pli]<<"\n\n"; exit(1); }
}

// print error msg after keyword for erroneous keyword statement+msg
void prt_akw_err_msg(char **Pgm,const char *kw,int Pli,const char *msg)
{
   cerr<<"\nERROR: Invalid statement for "<<msg<<":\n"<<Pgm[Pli]<<"\n\n"; exit(1);
}

// print error msg after keyword for erroneous keyword statement or wrong array size
void prt_akw_errSz(char **Pgm,const char *kw,int Pli,int szData,int szTrue)
{
   if (szData>=1)
   {
      if (szData!=szTrue) { cerr<<"\nERROR: Invalid array size: "<<Pgm[Pli]<<"\n\n"; exit(1); }
      else                { cerr<<"\nERROR: Invalid statement: "<<Pgm[Pli]<<"\n\n"; exit(1); }
   }
}


// check existence of kw from pgm lines
int chk_kw(char **Pgm,int Npl,const char *kw)
{
   int fd,p; char pln[512]="",pkw[512]="";

   fd=0;
   for (p=0;p<Npl;p++)
   {
      strcpy_up(pln,Pgm[p]);
      strcpy_bdlm(pkw,pln,"=");
      if (index(pkw,kw)==1) fd=1;
   }
   return fd;
}


// check aftKeyWdFmPgmLn for each str in strArrayIpt
int chk_akw_sa(char **Pgm,int Npl,const char *kw,const char *saI)
{
   int saiSz,saoSz; char **saO; Array2(saO,10,512,"");
   int fd,p,i,j; char pln[512]="",pkw[512]="",ikw[50]=""; 

   saiSz=wordCount(saI);
   fd=0; saoSz=0;
   for (p=0;p<Npl;p++)
   {
      strcpy_up(pln,Pgm[p]);
      strcpy_bdlm(pkw,pln,"=");
      if (strcmp(pkw,kw)==0)
      {
         ipt_avn_sa(pln,kw,saO,saoSz);
         if (saoSz>=1)
         {
            for (i=0;i<saoSz;i++) for (j=0;j<saiSz;j++)
            {
               ipt_nth_word(saI,j+1, ikw);
               if (strcmp(saO[i],ikw)==0)
               {
                  fd=1; goto end_p;
               }
            }
         }
      }
   }
end_p:
   Drray2(saO,10,512);
   return fd;
}

// check aftKeyWdFmPgmLn for any str in strArrayIpt
int chk_akw_saAny(char **Pgm,int Npl,const char *kw,const char *saI)
{
   int saiSz,saoSz; char **saO; Array2(saO,10,512,"");
   int fd,p,i,j; char pln[512]="",pkw[512]="",ikw[50]=""; 

   saiSz=wordCount(saI);
   fd=0; saoSz=0;
   for (p=0;p<Npl;p++)
   {
      strcpy_up(pln,Pgm[p]);
      strcpy_bdlm(pkw,pln,"=");
      if (strcmp(pkw,kw)==0)
      {
         ipt_avn_sa(pln,kw,saO,saoSz);
         if (saoSz>=1)
         {
            for (i=0;i<saoSz;i++) for (j=0;j<saiSz;j++)
            {
               ipt_nth_word(saI,j+1, ikw);
               if (index(saO[i],ikw)==1) // if (strcmp(saO[i],ikw)==0)
               {
                  fd=1; goto end_p;
               }
            }
         }
      }
   }
end_p:
   Drray2(saO,10,512);
   return fd;
}


// input after key word from program lines for the first 1+ strings
void ipt_akw_s1(char **Pgm,int Npl,const char *kw,char *s1,int &Pli)
{
   int saSz,p; char **sa,c[512]="",w1[512]=""; Array2(sa,10,512,"");

   strcpy(s1,""); saSz=0;
   for (p=0;p<Npl;p++)
   {
      strcpy_up(c,Pgm[p]); strcpy_bdlm(w1,c,"=");
      if (strcmp(w1,kw)==0)
      {
         ipt_avn_sa(c,kw,sa,saSz);
         if (saSz>=1) { strcpy(s1,sa[0]); Pli=p; }
      }
   }
   Drray2(sa,10,512);
}
void ipt_akw_s2(char **Pgm,int Npl,const char *kw,char *s1,char *s2,int &Pli)
{
   int saSz,p; char **sa,c[512]="",w1[512]=""; Array2(sa,10,512,"");

   strcpy(s1,""); strcpy(s2,""); saSz=0;
   for (p=0;p<Npl;p++)
   {
      strcpy_up(c,Pgm[p]); strcpy_bdlm(w1,c,"=");
      if (strcmp(w1,kw)==0)
      {
         ipt_avn_sa(Pgm[p],kw,sa,saSz);
         if (saSz>=1)
         {
            strcpy(s1,sa[0]);
            if (saSz>=2) strcpy(s2,sa[1]);
            Pli=p;
         }
      }
   }
   Drray2(sa,10,512);
}
void ipt_akw_s3(char **Pgm,int Npl,const char *kw,char *s1,char *s2,char *s3,int &Pli)
{
   int saSz,p; char **sa,c[512]="",w1[512]=""; Array2(sa,10,512,"");

   strcpy(s1,""); strcpy(s2,""); strcpy(s3,""); saSz=0;
   for (p=0;p<Npl;p++)
   {
      strcpy_up(c,Pgm[p]); strcpy_bdlm(w1,c,"=");
      if (strcmp(w1,kw)==0)
      {
         ipt_avn_sa(c,kw,sa,saSz);
         if (saSz>=1)
         {
            strcpy(s1,sa[0]);
            if (saSz>=2) strcpy(s2,sa[1]);
            if (saSz>=3) strcpy(s3,sa[2]);
            Pli=p;
         }
      }
   }
   Drray2(sa,10,512);
}
void ipt_akw_s4(char **Pgm,int Npl,const char *kw,char *s1,char *s2,char *s3,char *s4,int &Pli)
{
   int saSz,p; char **sa,c[512]="",w1[512]=""; Array2(sa,10,512,"");

   strcpy(s1,""); strcpy(s2,""); strcpy(s3,""); strcpy(s4,""); saSz=0;
   for (p=0;p<Npl;p++)
   {
      strcpy_up(c,Pgm[p]); strcpy_bdlm(w1,c,"=");
      if (strcmp(w1,kw)==0)
      {
         ipt_avn_sa(c,kw,sa,saSz);
         if (saSz>=1)
         {
            strcpy(s1,sa[0]);
            if (saSz>=2) strcpy(s2,sa[1]);
            if (saSz>=3) strcpy(s3,sa[2]);
            if (saSz>=4) strcpy(s4,sa[3]);
            Pli=p;
         }
      }
   }
   Drray2(sa,10,512);
}

// get data size=(no. of observations)
int getDsz(char *fn,int vhead)
{
   int i1,i2,dsz;
   getRsz(fn, i1,i2); dsz=i1-(vhead==1?1:0)-i2;
   return dsz;
}

// input varArrayFromDataFileFn for 1 int var
void ipt_var_i1(char *fn,int vhead,int vn, int *v)
{
   int i,j; double d1; char Ln[LS],c1[20];
   ifstream din(fn,ios::in); if (!din) errOpenIF_exit(fn);

   if (vhead==1) din.getline(Ln,LS,'\n'); i=0;
   while (din>>c1)
   {
      if (strcmp(c1,"")==0) break;
      for (j=1;j<=vn-2;j++) din>>d1;
      if (vn==1) v[i]=atoi(c1); else din>>v[i]; i++;
      din.getline(Ln,LS,'\n');
   }
}

// input varArrayFromDataFileFn for 1 dbl var
void ipt_var_d1(char *fn,int vhead,int vn, double *v)
{
   int i,j; double d1; char Ln[LS],c1[20];
   ifstream din(fn,ios::in); if (!din) errOpenIF_exit(fn);

   if (vhead==1) din.getline(Ln,LS,'\n'); i=0;
   while (din>>c1)
   {
      if (strcmp(c1,"")==0) break;
      for (j=1;j<=vn-2;j++) din>>d1;
      if (vn==1) v[i]=atof(c1); else din>>v[i]; i++;
      din.getline(Ln,LS,'\n');
   }
}


// [4.4] OPEN+MANAGE FILE

// open file (OK for ofstream &out or ofstream out)
void openOF(ofstream &out, const char *fn)
{
   out.open(fn, ios::out); if (!out) errOpenOF_exit(fn);
   //out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4);
}
void openOF(ofstream &out, char *fn)
{
   out.open(fn, ios::out); if (!out) errOpenOF_exit(fn);
   //out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4);
}

void openOFNewExt(ofstream &out, const char *fn, char *sExt)
{
   int i1=(int)strlen(fn), i2=(int)strlen(sExt);
   char *fnJk; Array1(fnJk, i1+1);
   strncpy(fnJk, fn, i1); strncpy(fnJk+i1-i2, sExt, i2); fnJk[i1]='\0';
   out.open(fnJk, ios::out); if (!out) errOpenOF_exit(fnJk);
}

void openOFAddExt(ofstream &out, const char *fn, char *sExt)
{
   int i1=(int)strlen(fn), i2=(int)strlen(sExt);
   char *fnJk; Array1(fnJk, i1+i2+1);
   strncpy(fnJk, fn, i1); strncpy(fnJk+i1, sExt, i2); fnJk[i1+i2]='\0';
   out.open(fnJk, ios::out); if (!out) errOpenOF_exit(fnJk);
}

void openOFNewExtInt(ofstream &out, const char *fn, int iExt)
{
   int i1=(int)strlen(fn);
   char *sIn; Array1(sIn, i1+1);
   strncpy(sIn, fn, i1); sIn[i1]='\0';
   char *sOut; Array1(sOut, i1+1);
   strcpy_newiExt(sOut,sIn,iExt);
   out.open(sOut, ios::out); if (!out) errOpenOF_exit(sOut);
}

void openOFAddExtInt(ofstream &out, const char *fn, int iExt)
{
   int i1=(int)strlen(fn), i2=4;
   char *sIn; Array1(sIn, i1+1);
   strncpy(sIn, fn, i1); sIn[i1]='\0';
   char *sOut; Array1(sOut, i1+i2+1);
   strcpy_addiExt(sOut, sIn, iExt);
   out.open(sOut, ios::out); if (!out) errOpenOF_exit(sOut);
}

void openOF_addSuf_keepE(ofstream &out,const char *fn,char *sSuf)
{
   int i1=(int)strlen(fn), i2=(int)strlen(sSuf);
   char *fnJk; Array1(fnJk, i1+i2+1);
   strcpy_addC_keepE(fnJk,fn,sSuf);
   out.open(fnJk, ios::out); if (!out) errOpenOF_exit(fnJk);
}

void apndOF(ofstream &out, const char *fn)
{
   out.open(fn, ios::app); if (!out) errOpenOF_exit(fn);
   out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4);
}
void apndOF(ofstream &out, char *fn)
{
   out.open(fn, ios::app); if (!out) errOpenOF_exit(fn);
   out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4);
}

void apndOFNewExt(ofstream &out, const char *fn, char *sExt)
{
   int i1=(int)strlen(fn), i2=(int)strlen(sExt);
   char *fnJk; Array1(fnJk, i1+1);
   strncpy(fnJk, fn, i1); strncpy(fnJk+i1-i2, sExt, i2); fnJk[i1]='\0';
   out.open(fnJk, ios::app); if (!out) errOpenOF_exit(fnJk);
   out<<setiosflags(ios::fixed|ios::showpoint)<<setprecision(4);
}

// open or append file co mode of output file (o1a2) 1=open, 2=apnd
void openApndOF(ofstream &out, const char *fn,int o1a2)
{
   if (o1a2==1) openOF(out, fn); else if (o1a2==2) apndOF(out,fn);
}
void openApndOF(ofstream &out, char *fn,int o1a2)
{
   if (o1a2==1) openOF(out, fn); else if (o1a2==2) apndOF(out,fn);
}

// open file with added characters to filename (Filename.Ext -> Filename||sNew||.Ext)
void openApndOF_addC_keepE(ofstream &out, const char *fn,char *sNew,int o1a2)
{
   char fnNew[100]="";
   strcpy_addC_keepE(fnNew, fn,sNew);
   openApndOF(out, fnNew, o1a2);
}


// append lines n1l-n1u in file f1 followed by lines n2l-n2u in file f2
void apnd2Files(char *f1,char *f2,char *fout,int n1l,int n1u,int n2l,int n2u)
{
   char Ln[LS]; char **f12; int n,n1,n2,len,len_max;

   ifstream cin1(f1, ios::in); if (!cin1) errOpenIF_exit(f1);
   ifstream cin2(f2, ios::in); if (!cin2) errOpenIF_exit(f2);
   len_max=0; n1=0; n2=0;
   while (cin1.getline(Ln,LS,'\n')) { n1++; len=(int)strlen(Ln); if (len>len_max) len_max=len; }
   while (cin2.getline(Ln,LS,'\n')) { n2++; len=(int)strlen(Ln); if (len>len_max) len_max=len; }

   Array2(f12,n1+n2,len_max+1);
   n=0;
   // note: "cin1.clear(); cin1.seekg(0);" is needed only if the data is being inputted again
   cin1.clear(); cin1.seekg(0); while (cin1.getline(Ln,LS,'\n')) { strcpy(f12[n],Ln); n++; }
   cin2.clear(); cin2.seekg(0); while (cin2.getline(Ln,LS,'\n')) { strcpy(f12[n],Ln); n++; }

   ofstream out; openOF(out,fout);
   for (n=0; n<n1+n2; n++)
   {
      if (n1u==0 || (n>=n1l-1 && n<=n1u-1) || (n>=n1+n2l-1 && n<=n1+n2u-1)) out<<f12[n]<<endl;
   }
   Drray2(f12,n1+n2,len_max+1);
}

// append file f1 and file f2
void apnd2Files(char *f1,char *f2,char *fout)
{
   apnd2Files(f1,f2,fout,0,0,0,0);
}


// [4.5] CALCULATION FUNCTIONS

// [4.5a] SIMPLE FUNCTIONS AND MATRIX COMPUTATION
// get sign of x
int sign(double x) { return (x>0.?1:(x<0.?-1:0)); }

// get square root of an integer
double sqrt(int i) { return sqrt((double)i); }

// logit(p) = log(p/(1-p)):
double logit(double p) { return log(p/(1-p)); }

// logistic or invlogit(x) = 1/(1+exp(-x)):
double logistic(double x) { return 1./(1.+exp(-x)); }



// update double array for const sum
void sumArrayConst(double *a,int s1,double csum)
{
   int i; double d1;
   for (i=0;i<s1;i++) if (a[i]<0) { cerr<<"\nERROR: sumArrayConst() failed"
      <<" for array value out of range.\n\n"; exit(1); }
   d1=Sum(a,s1);
   if (d1<csum-1e-50 || d1>csum+1e-50) { d1=csum/d1; for (i=0;i<s1;i++) a[i]*=d1; }
}

// get a complementary sub-matrix associated w/ element a[i][j]
void MatSub(double **a, int aSz, int i, int j, double **aSub)
{
   int k, l, kFull, lFull;

   for(k=0; k<aSz-1; k++)
   {
      kFull = (k<i?k:k+1);
      for(l=0; l<aSz-1; l++)
      {
         lFull = (l<j?l:l+1);
         aSub[k][l] = a[kFull][lFull];
      }
   }
}

// get determinant value of matrix a[aSz]
double MatDet(double **a, int aSz)
{
   int j;
   double det=-1;

   if      (aSz==1) det = a[0][0];
   else if (aSz==2) det = a[0][0]*a[1][1]-a[0][1]*a[1][0];
   else if (aSz==3) det = a[0][0]*( a[1][1]*a[2][2] - a[1][2]*a[2][1] )
                        - a[0][1]*( a[1][0]*a[2][2] - a[1][2]*a[2][0] )
                        + a[0][2]*( a[1][0]*a[2][1] - a[1][1]*a[2][0] );
   else if (aSz>=4)
   {
      det = 0.0;
      double **aSub; Array2(aSub,aSz-1,aSz-1);
      for (j=0; j<aSz; j++)
      {
         MatSub(a, aSz, 0, j, aSub);
         det += a[0][j] * MatDet(aSub, aSz-1) * (1-(j%2)*2);
      }
      Drray2(aSub,aSz-1,aSz-1);
   }

   return det;
}

// get the inverse matrix aInv[aSz][aSz]
void MatInv(double **a, double **aInv, int aSz)
{
   int i, j;
   double det;
   double **aSub; Array2(aSub,aSz-1,aSz-1);

   det = MatDet(a, aSz);
   if (fabs(det)<1e-20) {
      cerr<<"\nERROR: Double value out of range: MatInv() w/ 0 determinant.\n\n"; exit(1); }

   for(i=0; i<aSz; i++)
   {
      for(j=0; j<aSz; j++)
      {
         MatSub(a, aSz, j, i, aSub);
         aInv[i][j] = (1-((i+j)%2)*2)*MatDet(aSub, aSz-1)/det;
      }
   }

   Drray2(aSub,aSz-1,aSz-1);
}


// get no. of consecutive shared intervals surrounding a mk b/ 2 haplotypes
int sharedItv(int *h1, int *h2, int M, int m0)
{
   int m,nSI=0;
   if (h1[m0]==h2[m0] && m0<=M-2)
   {
      for (m=m0;m<M;m++) { if (h1[m]==h2[m] && h1[m+1]==h2[m+1]) nSI++; else break; }
   }
   if (h1[m0]==h2[m0] && m0>=1)
   {
      for (m=m0;m>=1;m--) { if (h1[m]==h2[m] && h1[m-1]==h2[m-1]) nSI++; else break; }
   }
   return nSI;
}



// STATISTICAL TESTS: FISHER'S EXACT TEST, 2 GROUP COMPARISONS, ETC
// get log-factorial of n (i.e., log(n!))
double lfact(int n)
{
   int i; double lf; lf=0; for (i=2; i<=n; i++) lf += log((double)i); return lf;
}

// get hypergeometric prob of a 2 by 2 table (a,b,c,d) for Fisher's exact test
double pFisherEtTab(int a,int b,int c,int d)
{
   int n; double lp,p; n=a+b+c+d;
   lp = lfact(a+b)+lfact(c+d)+lfact(a+c)+lfact(b+d)-lfact(n)-lfact(a)-lfact(b)-lfact(c)-lfact(d);
   p=exp(lp); return p;
}

// get 2-sided p-value of Fisher's exact test
double pFisherEt2Side(int a0,int b0,int c0,int d0)
{
   int n1,n2,n3,n4, nMi,iMi, i,a,b,c,d; double pObs,p2Side,pTab;

   // get nMi=min(n1,n2,n3,n4), iMi=index of nMi
   n1=(a0+b0); n2=(c0+d0); n3=(a0+c0); n4=(b0+d0); nMi=n1; iMi=1;
   if (n2<nMi) { nMi=n2; iMi=2; } if (n3<nMi) { nMi=n3; iMi=3; } if (n4<nMi) { nMi=n4; iMi=4; }

   // sum over all table probs w/ less than or equal to that of observed table
   pObs=pFisherEtTab(a0,b0,c0,d0); p2Side=0;
   if (nMi>=1)
   {
      for (i=0;i<=nMi;i++)
      {
         if      ((iMi==1 && a0<=b0) || (iMi==3 && a0<=c0)) { a=i; b=n1-i; c=n3-i; d=n4-n1+i; }
         else if ((iMi==1 && b0< a0) || (iMi==4 && b0<=d0)) { a=n1-i; b=i; c=n3-n1+i; d=n4-i; }
         else if ((iMi==2 && c0<=d0) || (iMi==3 && c0< a0)) { a=n3-i; b=n1-n3+i; c=i; d=n2-i; }
         else if ((iMi==2 && d0< c0) || (iMi==4 && d0< b0)) { a=n1-n4+i; b=n4-i; c=n2-i; d=i; }
         if (Min(a,b,c,d)>=0)
         {
            pTab=pFisherEtTab(a,b,c,d); if (pTab<=pObs) p2Side+=pTab;
         }
      }
   }
   else if (nMi==0) { p2Side=pObs; }
   return p2Side;
}


// RANK, PERCENTAGE AND PERCENTILE
// get rank of x in xs[n] from largest
int RankLg(double *xs,int n,double x)
{
   int i,nPre,nDp; nPre=0; nDp=0;
   for (i=0;i<n;i++) { if (xs[i]>x) nPre++; else if (xs[i]==x) nDp++; }
   return nPre+(nDp+1)/2;
}

// get rank of x in xs[n] from smallest
int RankSm(double *xs,int n,double x)
{
   int i,nPre,nDp; nPre=0; nDp=0;
   for (i=0;i<n;i++) { if (xs[i]<x) nPre++; else if (xs[i]==x) nDp++; }
   return nPre+(nDp+1)/2;
}


// get percentage of x in array xs[n]
double Percentage(double *xs,int n,double x)
{
   int j,nLt,nDp; double xmi,xma,d1,xLt,xRt,pctg;

   nLt=0; nDp=0; xmi=Min(xs,n); xma=Max(xs,n); xLt=xmi; xRt=xma;
   if (x>=xmi && x<=xma)
   {
      for (j=0;j<n;j++)
      {
         d1=xs[j];
         if       (d1<x-1e-10) { nLt++; xLt=Max(xLt,d1); }
         else if  (d1<x+1e-10) { nDp++; }
         else                  { xRt=Min(xRt,d1); }
      }
      if (nDp>=1) pctg=((double)nLt+(double)nDp/2)/n;
      else        pctg=((double)nLt-0.5+(double)(x-xLt)/(xRt-xLt))/n;
   }
   else if (x<xmi) pctg=0.5/n;
   else if (x>xma) pctg=(n-0.5)/n;

   return pctg;
}

// get percentage array pctg[n] (i.e., quantile of rank value) atEachValueOf xs[]
void pctgVec(double *xs,int n, double *pctg)
{
   int i,j,nLt,nDp;

   for (i=0;i<n;i++)
   {
      nLt=0; nDp=0;
      for (j=0;j<n;j++)
      {
         if      (xs[j]<xs[i]-1e-10) nLt++;
         else if (xs[j]<xs[i]+1e-10) nDp++;
      }
      pctg[i]=((double)nLt+(double)nDp/2)/n;
   }
}

// get correlation coefficient btw xs[] and its percentage pctg[]
double pctgCorr(double *xs,int n)
{
   double r,*pctg; Array1(pctg,n);

   pctgVec(xs,n,pctg);
   r=Corr(xs,pctg,n); Drray1(pctg,n);

   return r;
}


// get percentile of p in array xs[]
double Percentile(double *xs,int n, double p)
{
   int i; double wi,pctl, *b; Array1(b, n);

   if (p<0 || p>1) { cerr<<"\nERROR: Percentile(*,*,"<<p<<") failed.\n\n"; exit(1); }
   passArray(xs,b,n); sortArrayAsc(b,n);
   i=int(p*n); wi=1.-(p-(double)i/n)*n;
   if      (i>=1 && i<=n-1) pctl=wi*b[i-1]+(1-wi)*b[i];
   else if (i==0) pctl=b[i];
   else if (i==n) pctl=b[n-1];

   Drray1(b, n);
   return pctl;
}


// SIMPLE ALGORITHMS

// Newton-Raphson method forFindingTheRootOf f(x)=0 with fd(x)=derivativeOf_f(x)
// note: x0=initialGuessOfRootX, eps=convergenceCriterion
double newton_raphson(double f(double),double fd(double),double x0,double eps=1e-4,int itm=100)
{
   int it; double xpre,x;

   xpre=x0; it=0;
   while (it<itm)
   {
      x=xpre-f(xpre)/fd(xpre);
      if (fabs(x-xpre)<eps) break; else { xpre=x; it++; }
   }

   return x;
}


#endif

// Line(date): 2333(9/10/2009 for Windows & Linux), 2302(9/15),2516(10/3),2733(12/28).
// 2912(5/13/2010),2937(8/13),2954(8/18).
// 3040(4/8/2011),3273(6/13),3365(6/15),2611(6/16),2625(6/19),2698(6/26),2764(6/30),2873(7/4).
// 1972(7/13/2012)
