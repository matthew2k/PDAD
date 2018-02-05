/*
PROGRAM: f9_asso.h
PURPOSE: COMMON FUNCTIONS for association analysis
*/

#ifndef F9_ASSO_H
#define F9_ASSO_H


///////////////////////////////////////////////////////
// FUNCTIONS #9 (asso_f#x_#x) for association analysis
///////////////////////////////////////////////////////

double p_dad(double *xs, int n, double x, char *dig, int L0R1);                   // asso_f1
double p_dad_detail(double *xs, int n, double x, char *dig, int &nd, int *di,     // asso_f1a
	double *w, double &fbst);
void getdp_xsmv_mom(double *xs, int n, double m, double v, int id, double *dp);   // asso_f1b
double getd_dp(double x, double *dp, int id);                                   // asso_f1c
void hooke_dadml(double *sta, int np, int nps, double *par, double &fB);         // asso_f1d
double getp_dp(double x, double *dp, int id);                                   // asso_f1e

double p_gum(double *xs, int n, double x, int L0R1);                             // asso_f2a
double p_gpd(double *xs, int n, double x, int L0R1);                             // asso_f2b

void hooke_minf(double *par, int np, double *sta, int ns, double *dlt, double *pmi,// asso_f3a
	double *pma, double rho, double eps, int itm, double fo(double *, int, double *, int));


///////////////////////////////////////////////////////
// FUNCTIONS #9 (asso_f#x_#x) for association analysis
///////////////////////////////////////////////////////

// asso_f1: get p-value by fitting a weighted sum distribution co xs[]+dig="1,2,3,4,5"
double p_dad(double *xs, int n, double x, char dig[100] = "", int L0R1 = 1)
{
	int ndm, nd, *di;
	double p, fbst, *w;
	cout << "1 ";
	ndm = 10;
	Array1(di, ndm);
	Array1(w, ndm);
	cout << "2 ";
	p = p_dad_detail(xs, n, x, dig, nd, di, w, fbst);
	cout << "3 ";
	if (L0R1 == 0) 
	{
		p = 1 - p; // asso_f1a
	}
	Drray1(di, ndm);
	Drray1(w, ndm);

	return p;
}


// asso_f1a: get right-tail pval & detailed par=di[]+w[] byFitWtSumDen co xs[]+dig="1,2,3,4,5"
double p_dad_detail(double *xs, int n, double x, char *dig, int &nd, int *di,
	double *w, double &fbst)
{
	int ndm, i, j, i1, err; 
	double m, v, pf, mi, ma, md, we, pe, p, **dp2, *ds; 
	char c1[10] = ""; ndm = 10;
	Array2(dp2, ndm, 4);
	Array1(ds, ndm*n);
	cout << "3a ";
	// step 1: estimate distribution specific paramters: nd+di[] & dp2[nd][4]+ds[nd*n]
	if (strcmp(dig, "") != 0) 
		ipt_ia(dig, di, nd); 
	else 
	{ 
		nd = 5; 
		for (j = 0; j<nd; j++) 
			di[j] = j + 1; 
	}
	cout << "4 ";
	err = 0; 
	cout << "5 ";
	for (j = 0; j<nd; j++) 
	{ 
		iToC(di[j], c1); 
		if (indexWd("1,2,3,4,5", c1) == 0) err = 1; 
	}
	cout << "6 ";
	if (err == 1) 
	{ 
		cerr << "\nERROR: Unable to compute p_dad for dig=" << dig << ".\n\n"; 
		exit(1); 
	}
	cout << "7 ";
	m = Mean(xs, n); 
	v = Var(xs, n);
	
	for (j = 0; j<nd; j++)
		getdp_xsmv_mom(xs, n, m, v, di[j], dp2[j]); // asso_f1b
	i1 = nd; 
	nd = 0; // update nd+di[]+dp2[][]+ds[] as needed
	cout << "8 ";
	for (j = 0; j<i1; j++)
	{
		cout << "8a ";
		if (dp2[j][0]>0)
		{
			di[nd] = di[j]; 
			for (i = 0; i<4; i++) 
				dp2[nd][i] = dp2[j][i];
			cout << "8b ";
			for (i = 0; i < n; i++)
				cout << " 8c";
				ds[nd*n + i] = getd_dp(xs[i], dp2[nd], di[nd]); // asso_f1c
			nd++;
		}
	}
	cout << "9 ";
	// step 2: estimate distribution specific weights using HookeAlgo+ML
	hooke_dadml(ds, nd, nd*n, w, fbst); // asso_f1d

	// step 3: update rt-tail p-value as a weighted sum of percentage and fitted p-value
	pf = 0.; 
	for (j = 0; j<nd; j++) 
		pf += (w[j] * getp_dp(x, dp2[j], di[j])); // asso_f1e
	cout << "10 ";
	mi = Min(xs, n); 
	ma = Max(xs, n);
	
	if (x<mi || x>ma) 
		p = 1. - pf;
	else
	{
		pe = Percentage(xs, n, x); md = Median(xs, n); we = (x<md ? (x - mi) / (md - mi) : (ma - x) / (ma - md));
		p = 1. - (we*pe + (1. - we)*pf);
	}
	cout << "11 ";
	Drray2(dp2, ndm, 4); Drray1(ds, ndm*n);
	return p;
}


// asso_f1b: get dp[]=distPar co xs[n]+m+v+id byMtdOfMoment. note: useXs[]ForFastCompWbPar
void getdp_xsmv_mom(double *xs, int n, double m, double v, int id, double *dp)
{
	double a, b, r, cv;

	if (id >= 1 && id <= 4)
	{
		if (m>0 && v>0)
		{
			if (id == 1) // chisq(): dp[0]=df, dp[1]=nc
			{
				if (m>v / 4 && m <= v / 2) { dp[0] = 2 * m - v / 2; dp[1] = v / 2 - m; }
				else { dp[0] = -9.; dp[1] = -9.; }
			}
			else if (id == 2) // gamma(): dp[0]=shape, dp[1]=scale=1/rate
			{
				dp[0] = m*m / v; dp[1] = v / m;
			}
			else if (id == 3) // lnorm(): dp[0]=logNormMu, dp[1]=logNormSigma
			{
				dp[0] = log(m / sqrt((v + m*m) / (m*m))); dp[1] = sqrt(log((v + m*m) / (m*m)));
			}
			else if (id == 4) // weibull(): dp[0]=shape, dp[1]=scale
			{
				r = pctgCorr(xs, n); cv = sqrt(v) / m;
				a = -log(2.) / log(1. - r*cv*sqrt((double)(n + 1) / ((n - 1) * 3))); b = m / exp(lgamma(1. + 1. / a));
				dp[0] = a; dp[1] = b;
			}
		}
		else {
			cerr << "\nERROR: getdp_xsmv_mom() failed for id=" << id << ", m=" << m << ", v=" << v
				<< ".\n\n"; exit(1);
		}
	}
	else if (id == 5) // gumbel(): dp[0]=location, dp[1]=scale
	{
		// gm=Euler-Mascheroni constant=0.57721566490153286060651209008240243104215933593992
		b = sqrt(6 * v) / PI; a = m - 0.5772156649*b; dp[0] = a; dp[1] = b;
	}
	else { cerr << "\nERROR: getdp_xsmv_mom() failed for unknown id=" << id << ".\n\n"; exit(1); }
}


// asso_f1c: get density co x+dp[]+id
double getd_dp(double x, double *dp, int id)
{
	double d;
	if (id == 1)
	{
		d = dchisq(x, dp[0], dp[1]);
		//cout << "chisq ";
	}
	else if (id == 2)
	{
		d = dgamma(x, dp[0], dp[1]);
		//cout << "dgamma ";
	}
	else if (id == 3)
	{
		d = dlnorm(x, dp[0], dp[1]);
		//cout << "dlnorm ";
	}
	else if (id == 4)
	{
		d = dweibull(x, dp[0], dp[1]);
		//cout << "dweibull ";
	}
	else if (id == 5)
	{
		d = dgumbel(x, dp[0], dp[1]);
		//cout << "dgumbel ";
	}
	else 
	{ 
		cerr << "\nERROR: getd_dp() failed for unknown id=" << id << ".\n\n"; exit(1); 
	}
	return Max(d, 1e-50);
}


// asso_f1d: get par[] by Hooke+Jeeves algo (1961, JACM 8:212-229) co wtSumDen+minNll
void hooke_dadml(double *sta, int np, int nps, double *par, double &fB)
{
	double rho, eps; int itm; rho = 0.5; eps = 1e-6; itm = 1000; // constParForHookeAlgo
	int ns, it, i, j, k, jP, sP, jB, sB, doit; // xP=xPre, xB=xBst, sta[nps]=denSeriesForMulDist
	double pmi, pma, psu, stepsz, f, fP, d1, c1, c2, c3, p, *dlt, *da;
	Array1(dlt, np); Array1(da, np);

	pmi = 0.; pma = 1.; psu = 1.; for (j = 0; j<np; j++) dlt[j] = 0.5;
	it = 0; ns = nps / np; stepsz = sqrt(SS(dlt, np)); jP = -9; sP = -9;
	for (j = 0; j<np; j++)
	{
		for (k = 0; k<np; k++) da[k] = (k == j ? pma : pmi);
		f = 0.; for (i = 0; i<ns; i++) {
			d1 = 0.; for (k = 0; k<np; k++) d1 += sta[k*ns + i] * da[k]; f += -log(d1);
		}
		if (j == 0 || f<fP) { passArray(da, par, np); fP = f; }
	}
	sumArrayConst(par, np, psu);
	while (it<itm && stepsz >= eps)
	{
		fB = fP; jB = -9; sB = -9;
		for (j = 0; j<np; j++)
		{
			c1 = dlt[j]; c2 = (psu - par[j] - c1) / (psu - par[j]); c3 = (psu - par[j] + c1) / (psu - par[j]);
			doit = 1; for (k = 0; k<np; k++)
			{
				if (k == j) p = par[k] + c1; else p = par[k] * c2;
				if (p<pmi || p>pma) { doit = 0; break; }
				else da[k] = p;
			}
			if ((j != jP || sP != -1) && doit == 1)
			{
				f = 0.; for (i = 0; i<ns; i++) {
					d1 = 0.; for (k = 0; k<np; k++) d1 += sta[k*ns + i] * da[k]; f += -log(d1);
				}
				if (f<fB) { fB = f; jB = j; sB = 1; }
			}

			doit = 1; for (k = 0; k<np; k++)
			{
				if (k == j) p = par[k] - c1; else p = par[k] * c3;
				if (p<pmi || p>pma) { doit = 0; break; }
				else da[k] = p;
			}
			if ((j != jP || sP != 1) && doit == 1)
			{
				f = 0.; for (i = 0; i<ns; i++) {
					d1 = 0.; for (k = 0; k<np; k++) d1 += sta[k*ns + i] * da[k]; f += -log(d1);
				}
				if (f<fB) { fB = f; jB = j; sB = -1; }
			}
		}
		if (jB >= 0) // update par[] for better move (jB>=0)
		{
			c1 = dlt[jB] * sB; c2 = (psu - par[jB] - c1) / (psu - par[jB]);
			for (j = 0; j<np; j++) { if (j == jB) par[j] += c1; else par[j] *= c2; }
			fP = fB; jP = jB; sP = sB;
		}
		else { for (j = 0; j<np; j++) dlt[j] = dlt[j] * rho; stepsz = sqrt(SS(dlt, np)); jP = -9; sP = -9; }
		it++;
	}
	Drray1(dlt, np); Drray1(da, np);
}


// asso_f1e: get P-value=CDF of a known dist co dp[]+id
double getp_dp(double x, double *dp, int id)
{
	if (id >= 1 && id <= 4 && x <= 0.) return 1e-50;
	if (id == 1) return pchisq(x, dp[0], dp[1]);
	else if (id == 2) return pgamma(x, dp[0], dp[1]);
	else if (id == 3) return plnorm(x, dp[0], dp[1]);
	else if (id == 4) return pweibull(x, dp[0], dp[1]);
	else if (id == 5) return pgumbel(x, dp[0], dp[1]);
	else { cerr << "\nERROR: getp_dp() failed for unknown id=" << id << ".\n\n"; exit(1); }
}


// asso_f2a: get p=CDF by fitting a GumbelDist (2010, Abrams et al, IJHG 9:61)
double p_gum(double *xs, int n, double x, int L0R1 = 1)
{
	double m, v, a, b, p;
	m = Mean(xs, n); v = Var(xs, n); b = sqrt(6 * v) / PI; a = m - 0.5772156649*b;
	p = pgumbel(x, a, b); if (L0R1 == 1) p = 1 - p; return p;
}

// asso_f2b: get p=CDF byFitTailValToGenParetoDist (2009, Knijnenburg, Bioinfo 25:i161-i168)
double p_gpd(double *xs, int n, double x, int L0R1 = 1)
{
	int i, nGe, n1, n2, ne1, ne2, ne; double pemi, pema, texc, x1, m, v, a, b, p, *xe; Array1(xe, n);

	nGe = 0; for (i = 0; i<n; i++) { if (xs[i] >= x) nGe++; }
	n1 = 100; n2 = 5000; pemi = 0.05; pema = 0.5; ne1 = Round(n1*pema); ne2 = Round(n2*pemi);
	if (n >= n1 && n <= n2) ne = ne1 + Round((ne2 - ne1)*(log((double)n / n1)) / (log((double)n2 / n1)));
	else ne = (n <= n1 ? Round(n*pema) : Round(n*pemi));
	passArray_lgst(xs, n, xe, ne);
	texc = (Min(xe, ne) + MaxNth(xs, n, ne + 1)) / 2; for (i = 0; i<ne; i++) xe[i] -= texc; x1 = x - texc;
	m = Mean(xe, ne); v = Var(xe, ne); a = (m*m / v - 1) / 2; b = m*(m*m / v + 1) / 2;

	if (nGe >= 10 || x1<0 || (x1 >= 0 && a>0 && x1 >= b / a)) p = 1 - Percentage(xs, n, x);
	else if (x1 >= 0 && (a <= 0 || (a>0 && x1<b / a)))      p = ((double)ne / n)*(1. - pgpareto(x1, a, b));
	if (L0R1 == 0) p = 1 - p; p = ((p >= 0 && p <= 1) ? Max(p, 1e-50) : -9.);

	Drray1(xe, n);
	return p;
}


// asso_f3a: get par[]=bstParCoMinf by Hooke+Jeeves algo (1961, JACM 8:212-229)
void hooke_minf(double *par, int np, double *sta, int ns, double *dlt, double *pmi,
	double *pma, double rho, double eps, int itm, double fo(double *, int, double *, int))
{
	int it, jP, sP, jB, sB, j, k; double stepsz, fP, fB, d1, f, *da; // xP=xPre, xB=xBst
	Array1(da, np);

	it = 0; stepsz = sqrt(SS(dlt, np)); fP = fo(par, np, sta, ns); jP = -9; sP = -9;
	while (it<itm && stepsz >= eps)
	{
		fB = fP; jB = -9; sB = -9;
		for (j = 0; j<np; j++)
		{
			for (k = 0; k<np; k++) da[k] = par[k];
			d1 = par[j] + dlt[j]; if ((j != jP || sP != -1) && d1 <= pma[j])
			{
				da[j] = d1; f = fo(da, np, sta, ns); if (f<fB) { fB = f; jB = j; sB = 1; }
			}
			d1 = par[j] - dlt[j]; if ((j != jP || sP != 1) && d1 >= pmi[j])
			{
				da[j] = d1; f = fo(da, np, sta, ns); if (f<fB) { fB = f; jB = j; sB = -1; }
			}
		}
		if (jB >= 0) // update par[] for better move (jB>=0)
		{
			par[jB] = par[jB] + dlt[jB] * sB; fP = fB; jP = jB; sP = sB;
		}
		else { for (j = 0; j<np; j++) dlt[j] = dlt[j] * rho; stepsz = sqrt(SS(dlt, np)); jP = -9; sP = -9; }
		it++;
	}
	Drray1(da, np);
}

#endif
// Line(date): 413(5/22/2013), 374(5/30), 241(6/16), 346(8/6), 360(8/15), 279(9/12)
