#include "f1_head.h"

int main(int na, char *arg[])
{
   int nr,nc,ncMi,ncMa,rEmpty, r,n; double x,p,*xs;
   char fn[100]="",fn1[100]="",Ln[LS]; ofstream ou1;
   int decide;
   if (na==1)
   {
	   strncpy(fn, "input.txt", 100);
	   strncpy(fn1, "output.txt", 100);
	  
	   cout << "GUM = 1, PDAD = 2, GPD = 3?\n";
	   cin >> decide;
	   cout << "Running...";

	 
	   //cout<<"Enter input filename for pdad.exe: "; cin.getline(fn,100,'\n'); cout<<endl;
      //cout<<"Output filename for fitted pvalues: "; cin.getline(fn1,100,'\n'); cout<<endl;
   }
   else if (na==2) { strcpy(fn,arg[1]); strcpy_newE(fn1,fn,"pdad"); }
   else if (na>=3) { strcpy(fn,arg[1]); strcpy(fn1,arg[2]); }
   ifstream din(fn,ios::in); 
   if (!din) errOpenIF_exit(fn);
   openOF(ou1,fn1); 
   ou1<<setiosflags(ios::fixed|ios::showpoint);
   getRsz(fn, nr,rEmpty); 
   getCsz(fn, nc,ncMi,ncMa);
   Array1(xs,ncMa);
   cout << '0';
   
   if (decide == 1)
   {
	   for (r = 0; r < nr; r++)
	   { 
		   din >> x; din.getline(Ln, LS, '\n'); 
		   ipt_da(Ln, xs, n);
		   p = p_gum(xs, n, x);
		   cout << "\n";
		   cout << p;
		   ou1 << setprecision(10) << p << endl;
	   }
	   ou1 << "GUMBEL";
   }

   else if (decide == 2)
   {
	   for (r = 0; r < nr; r++)
	   {
		   din >> x; 		 
		   din.getline(Ln, LS, '\n'); 		  
		   ipt_da(Ln, xs, n);		
		   p = p_dad(xs, n, x);
		   cout << "\n";
		   cout << p;
		   ou1 << setprecision(10) << p << endl;
	   }
   }

   else if (decide == 3)
   {
	   for (r = 0; r < nr; r++)
	   {
		   din >> x; din.getline(Ln, LS, '\n'); 
		   ipt_da(Ln, xs, n);
		   p = p_gpd(xs, n, x);
		   cout << "\n";
		   cout << p;
		   ou1 << setprecision(10) << p << endl;
	   }
   }
   return 0;
}
