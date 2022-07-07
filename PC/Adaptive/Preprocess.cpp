#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

extern int Nphi, NBz0, NBz1, NBz2, NBz3, NBz4, NBz5, NBz6, NBz, NNBz, NNBt, Nch;
extern int NVt1, NVt2, NVt3, NVt4, NVt5, NVt6, NVr1, NVr2, NVr3, NVr4;
extern double dNBt, thetach;
extern double R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, revf, phif, phiout_c, hm, phiout_cy, phiout_dy, phiout_vb, phiout_FM, phiout_LK, w, wp;
extern double V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
extern double integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, zB1, zB2, zB3, zB4, zB5, zB6, zB7, wp, e0, theta0, psi0, de0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
extern double Lambda, Loc, Lpj, Lpg, Lbj, Ls, Rbo, Lpc;
extern float FxhD, FyhD, MxhD, MyhD, FxcD, FycD, MxcD, MycD, FxD, FyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, frD, MxD, MyD, MzD, leakageRinD, leakageRoutD, frictionD, friclossD, Hmin;

/* input data */
void Input()
{
	FILE *ep;

	if((ep = fopen("Input_pressure_asiatribo.dat", "r"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}

	fscanf(ep, "%*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*c %lf %*c %lf %*c %lf %*c %lf %*c %lf %*s %lf %*s %lf %*c %lf %*s %lf %*c %lf	\
				%*s %lf %*c %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*c %lf %*s %lf %*s %lf %*c %lf %*s %lf %*c %lf %*s %lf %*c %lf %*s %lf %*s %lf %*s %lf	\
				%*s %lf %*s %lf %*c %lf %*c %lf %*c %lf %*c %lf %*s %lf %*s %d %*c %d %*c %d %*c %d %*c %d %*c %d %*s %d %*s %lf %*c %lf %*c %lf %*c %lf %*c %lf %*c %lf %*s %d %*s %lf %*s %lf %*s %lf %*c %lf %*c %lf %*c %lf %*c %lf %*c %lf %*s %lf %*s %lf %*s %lf %*c %lf %*c %lf %*s %lf %*c %lf %*c %lf %*c %lf %*c %lf", 
		&R, &Rp, &Rpi, &Lpi, &Lc, &Ld, &thetaK, &thetaV1, &thetaV2, &thetaV3, &thetaV4, &thetaV5, &thetaV6, &rB, &rV1, &rV2, &LB1, &LB2, 
		&Cvo, &Cvi, &Cs, &Pin, &Pout, &Ph, &rho, &eta, &K, &Rs1, &Rs2, &Rso, &beta, &beta2, &Vvo, &Vvi, &Avout, &Avin, &Lvo, &Kvo, &Cvin, 
		&revf, &phiout_c, &hm, &phiout_cy, &phiout_dy, &phiout_vb, &w, &NBz0, &NBz1, &NBz2, &NBz3, &NBz4, &NBz5, &NNBt, &zB1, &zB2, &zB3, &zB4, &zB5, &zB6, &Nch, &thetach, 
		&wp, &e0, &theta0, &psi0, &de0, &wx0, &wy0, &Hg, &Mp, &Mb, &Iba, &Ibt, &Loc, &Lpj, &Lpg, &Lbj, &Ls);
	printf("%e %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e	\
			%e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e	\
			%e  %e  %e  %e  %e  %e  %e  %d  %d  %d  %d  %d  %d  %d  %e  %e  %e  %e  %e  %e  %d  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e %e \n",
			R, Rp, Rpi, Lpi, Lc, Ld, thetaK, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6, rB, rV1, rV2, LB1, LB2, 
			Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, 
			revf, phiout_c, hm, phiout_cy, phiout_dy, phiout_vb, w, NBz0, NBz1, NBz2, NBz3, NBz4, NBz5, NNBt, zB1, zB2, zB3, zB4, zB5, zB6, Nch, thetach, 
			wp, e0, theta0, psi0, de0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt, Loc, Lpj, Lpg, Lbj, Ls); 

	fclose(ep);

}
/* ParameteZs & Unit conveZsion & Non-dimensionalization */
void parameter()
{
	Ap = Pi*Rp*Rp;
	Vpi = Pi*Rpi*Rpi*Lpi;
	hp = 10.e-6;
	hs = 10.e-6;
	hv = 10.e-6;	
	
	double rB1 = 11.5e-3, rB2 = 13.75e-3, rB3 = 16.95e-3, rB4 = 18.80e-3, rB5 = 19.70e-3, rB6 = 22.4e-3;
	rB = 0.5*(rB3-rB2);
	thetaK = thetaK*Pi/180;
	thetaB = thetaK/2;
	thetaV1 = thetaV1*Pi/180;
	thetaV2 = thetaV2*Pi/180;
	thetaV3 = thetaV3*Pi/180;
	thetaV4 = thetaV4*Pi/180;
	thetaV5 = thetaV5*Pi/180;
	thetaV6 = thetaV6*Pi/180;

	w = w*2*Pi/60;
	wp = wp*2*Pi/60;
	beta = beta*Pi/180;
	beta2 = beta2*Pi/180;
	phif = revf*2*Pi;
	e0 = e0/hp;
	theta0 = theta0*Pi/180*Rp/hp;
	psi0 = psi0*Pi/180*Rp/hp;
	de0 = de0/hp/w;
	wx0 = wx0*Pi/180/w*Rp/hp;
	wy0 = wy0*Pi/180/w*Rp/hp;
	hm = hm*Pi/180;

	//printf("%d %d %d \n", NVt1, NVt5, NVr4);

	
	// Nondimensionalize Radii in the Barrel 
	Lambda = 6*eta*w*pow(Rp/hp,2);	

	zB6 = Lc - Ld - 2*R*tan(beta);
	zB7 = Lc - Ld;
	zB1 = zB1/Rp;
	zB2 = zB2/Rp;
	zB3 = zB3/Rp;
	zB4 = zB4/Rp;	
	zB5 = zB5/Rp;
	zB6 = zB6/Rp;
	zB7 = zB7/Rp;

	//Lc = Lc/Rp;
	//Ld = Ld/Rp;
	//Lpg = Lpg/Rp;
	//Lpj = Lpj/Rp;
	//R = R/Rp;

	// NBz: number of CV's , NNBz: Number of Node	
	NBz6 = int(2*(zB7-zB6)/(zB6-zB5)*NBz5);
	NBz = NBz0 + NBz1 + NBz2 + NBz3 + NBz4 + NBz5 + NBz6;
	NNBz = NBz + 1;
	dNBt = 2*Pi/NNBt;

	// Fs = Fs/Lambda/Rbo/Rbo;
	Hg = Hg/hp;
	Mb = Mb*w*w/Lambda/Rp/(Rp/hp);
	Mp = Mp*w*w/Lambda/Rp/(Rp/hp);
	Iba = Iba*w*w/Lambda/pow(Rp,3)/(Rp/hp);
	Ibt = Ibt*w*w/Lambda/pow(Rp,3)/(Rp/hp);
}

/* Initial Condition & Configuration */
void initial(float **Ppc, double *ym, double **Hpc, double **Hpc_gro, float **Zpc, float **Tpc, double *Pc, float *PCleakage)
{
	void groove(double **H_gro, float **Tpc);

	int i, j, k;

	V0 = Vpi + (Ld+2*R*tan(beta))*Ap;	// Initial cylinder volume

	// Initial cylinder pressure
	Pc[1] = 0.0e5;
	Pc[2] = Pout*1.;
	Pc[3] = Pout*1.;
	Pc[4] = Pout*1.;
	Pc[5] = Pout*1.;
	Pc[6] = Pin*0;
	Pc[7] = Pin*1;
	Pc[8] = Pin*1;
	Pc[9] = Pin*0;
	Pc[10] = Pout*1.0;
	Pc[11] = Pin*1.0;
	Pc[12] = 0.0000000015/*0.1*w*R*(tan(beta)*sin(2*Pi/9)+tan(beta2)/cos(beta)*cos(2*Pi/9*2))*Ap*/;
	//Pc[13] = Pout*1.0;

	// Initial pressure Distribution
	for(j=1; j<NNBt; j++){
		for(i=2; i<NNBz/2; i++){
			Ppc[i][j] = Pc[1]/Lambda;
		}
	}
	for(j=1; j<=NNBt; j++){
		for(i=NNBz/2; i<NNBz; i++){
			Ppc[i][j] = Pc[1]/20/Lambda;
		}
	}

	// Initial motion variables
	ym[1] = e0/4;
	ym[2] = -e0*sqrt(15.0)/4;
	ym[3] = theta0;
	ym[4] = psi0;
	ym[5] = de0;
	ym[6] = de0;
	ym[7] = wx0;
	ym[8] = wy0;
	double e0 = sqrt(ym[1]*ym[1]+ym[2]*ym[2]);
	double kai0 = atan2(ym[2],ym[1]);
	double gamma0 = sqrt(ym[3]*ym[3]+ym[4]*ym[4])/*acos(cos(ym[3])*cos(ym[4]))*Rp/hp*/;
	double zeta0 = -atan2(ym[3],ym[4])/*-atan2(tan(ym[3]),sin(ym[4]))*/;

	// Initial film thickness
	groove(Hpc_gro, Tpc);
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBz; i++){
			Hpc[i][j] = 1 - e0*sin(Tpc[i][j]+kai0)-Zpc[i][j]*gamma0*sin(Tpc[i][j]+zeta0) + Hpc_gro[i][j];			
		}
	}	

	// Initialize leakages
	for(i=1; i<=N;i++){		
		PCleakage[i] = 0.0;
	}
}

/* Grid & Coordinate */
void grid(float **Zpc, float **Zn, float **Zs, float **Zw, float **Ze, float **T, float **Tn, float **Ts, float **Tw, float **Te, float **dZ, float **dT, float **dN, float **dS, float **dW, float **dE)
{
	int i,j;

	// dr0
	for(j=1; j<=NNBt; j++){
		for(i=2; i<NBz0+1; i++){	
			dN[i][j] = (zB1-0.0)/NBz0;
			dS[i][j] = (zB1-0.0)/NBz0;
			dW[i][j] = 2*Pi/NNBt;
			dE[i][j] = 2*Pi/NNBt;
			dZ[i][j] = (dN[i][j]+dS[i][j])/2.0;
			dT[i][j] = (dW[i][j]+dE[i][j])/2.0;
		}	
		dN[1][j] = (zB1-0.0)/NBz0;
		dS[1][j] = TINY;
		dW[1][j] = 2*Pi/NNBt;
		dE[1][j] = 2*Pi/NNBt;
		dZ[1][j] = (dN[1][j]+dS[1][j])/2.0;
		dT[1][j] = (dW[1][j]+dE[1][j])/2.0;
		dN[NBz0+1][j] = (zB2-zB1)/NBz1;
		dS[NBz0+1][j] = (zB1-0.0)/NBz0;
		dW[NBz0+1][j] = 2*Pi/NNBt;
		dE[NBz0+1][j] = 2*Pi/NNBt;
		dZ[NBz0+1][j] = (dN[NBz0+1][j]+dS[NBz0+1][j])/2.0;
		dT[NBz0+1][j] = (dW[NBz0+1][j]+dE[NBz0+1][j])/2.0;
	}	
	// dr1
	for(j=1; j<=NNBt; j++){
		for(i=1; i<NBz1; i++){	
			dN[NBz0+1+i][j] = (zB2-zB1)/NBz1;
			dS[NBz0+1+i][j] = (zB2-zB1)/NBz1;
			dW[NBz0+1+i][j] = 2*Pi/NNBt;
			dE[NBz0+1+i][j] = 2*Pi/NNBt;
			dZ[NBz0+1+i][j] = (dN[i][j]+dS[i][j])/2.0;
			dT[NBz0+1+i][j] = (dW[i][j]+dE[i][j])/2.0;
		}
		dN[NBz0+NBz1+1][j] = (zB3-zB2)/NBz2;
		dS[NBz0+NBz1+1][j] = (zB2-zB1)/NBz1;
		dW[NBz0+NBz1+1][j] = 2*Pi/NNBt;
		dE[NBz0+NBz1+1][j] = 2*Pi/NNBt;
		dZ[NBz0+NBz1+1][j] = (dN[NBz0+NBz1+1][j]+dS[NBz0+NBz1+1][j])/2.0;
		dT[NBz0+NBz1+1][j] = (dW[NBz0+NBz1+1][j]+dE[NBz0+NBz1+1][j])/2.0;		
	}		
	// dr2
	for(j=1; j<=NNBt; j++){
		for(i=1; i<NBz2; i++){
			dN[NBz0+NBz1+1+i][j] = (zB3-zB2)/NBz2;
			dS[NBz0+NBz1+1+i][j] = (zB3-zB2)/NBz2;
			dW[NBz0+NBz1+1+i][j] = 2*Pi/NNBt;
			dE[NBz0+NBz1+1+i][j] = 2*Pi/NNBt;
			dZ[NBz0+NBz1+1+i][j] = (dN[NBz0+NBz1+1+i][j]+dS[NBz0+NBz1+1+i][j])/2.0;
			dT[NBz0+NBz1+1+i][j] = (dW[NBz0+NBz1+1+i][j]+dE[NBz0+NBz1+1+i][j])/2.0;
		}
		dN[NBz0+NBz1+1+NBz2][j] = (zB4-zB3)/NBz3;
		dS[NBz0+NBz1+1+NBz2][j] = (zB3-zB2)/NBz2;
		dW[NBz0+NBz1+1+NBz2][j] = 2*Pi/NNBt;
		dE[NBz0+NBz1+1+NBz2][j] = 2*Pi/NNBt;
		dZ[NBz0+NBz1+1+NBz2][j] = (dN[NBz0+NBz1+1+NBz2][j]+dS[NBz0+NBz1+1+NBz2][j])/2.0;
		dT[NBz0+NBz1+1+NBz2][j] = (dW[NBz0+NBz1+1+NBz2][j]+dE[NBz0+NBz1+1+NBz2][j])/2.0;
	}
	// dr3
	for(j=1; j<=NNBt; j++){
		for(i=1; i<NBz3; i++){
			dN[NBz0+NBz1+1+NBz2+i][j] = (zB4-zB3)/NBz3;
			dS[NBz0+NBz1+1+NBz2+i][j] = (zB4-zB3)/NBz3;
			dW[NBz0+NBz1+1+NBz2+i][j] = 2*Pi/NNBt;
			dE[NBz0+NBz1+1+NBz2+i][j] = 2*Pi/NNBt;
			dZ[NBz0+NBz1+1+NBz2+i][j] = (dN[NBz0+NBz1+1+NBz2+i][j]+dS[NBz0+NBz1+1+NBz2+i][j])/2.0;
			dT[NBz0+NBz1+1+NBz2+i][j] = (dW[NBz0+NBz1+1+NBz2+i][j]+dE[NBz0+NBz1+1+NBz2+i][j])/2.0;
		}
		dN[NBz0+NBz1+1+NBz2+NBz3][j] = (zB5-zB4)/NBz4;
		dS[NBz0+NBz1+1+NBz2+NBz3][j] = (zB4-zB3)/NBz3;
		dW[NBz0+NBz1+1+NBz2+NBz3][j] = 2*Pi/NNBt;
		dE[NBz0+NBz1+1+NBz2+NBz3][j] = 2*Pi/NNBt;
		dZ[NBz0+NBz1+1+NBz2+NBz3][j] = (dN[NBz0+NBz1+1+NBz2+NBz3][j]+dS[NBz0+NBz1+1+NBz2+NBz3][j])/2.0;
		dT[NBz0+NBz1+1+NBz2+NBz3][j] = (dW[NBz0+NBz1+1+NBz2+NBz3][j]+dE[NBz0+NBz1+1+NBz2+NBz3][j])/2.0;
	}
	// dr4
	for(j=1; j<=NNBt; j++){
		for(i=1; i<NBz4; i++){
			dN[NBz0+NBz1+1+NBz2+NBz3+i][j] = (zB5-zB4)/NBz4;
			dS[NBz0+NBz1+1+NBz2+NBz3+i][j] = (zB5-zB4)/NBz4;
			dW[NBz0+NBz1+1+NBz2+NBz3+i][j] = 2*Pi/NNBt;
			dE[NBz0+NBz1+1+NBz2+NBz3+i][j] = 2*Pi/NNBt;
			dZ[NBz0+NBz1+1+NBz2+NBz3+i][j] = (dN[NBz0+NBz1+1+NBz2+NBz3+i][j]+dS[NBz0+NBz1+1+NBz2+NBz3+i][j])/2.0;
			dT[NBz0+NBz1+1+NBz2+NBz3+i][j] = (dW[NBz0+NBz1+1+NBz2+NBz3+i][j]+dE[NBz0+NBz1+1+NBz2+NBz3+i][j])/2.0;
		}
		dN[NBz0+NBz1+1+NBz2+NBz3+NBz4][j] = (zB6-zB5)/NBz5;
		dS[NBz0+NBz1+1+NBz2+NBz3+NBz4][j] = (zB5-zB4)/NBz4;
		dW[NBz0+NBz1+1+NBz2+NBz3+NBz4][j] = 2*Pi/NNBt;
		dE[NBz0+NBz1+1+NBz2+NBz3+NBz4][j] = 2*Pi/NNBt;
		dZ[NBz0+NBz1+1+NBz2+NBz3+NBz4][j] = (dN[NBz0+NBz1+1+NBz2+NBz3+NBz4][j]+dS[NBz0+NBz1+1+NBz2+NBz3+NBz4][j])/2.0;
		dT[NBz0+NBz1+1+NBz2+NBz3+NBz4][j] = (dW[NBz0+NBz1+1+NBz2+NBz3+NBz4][j]+dE[NBz0+NBz1+1+NBz2+NBz3+NBz4][j])/2.0;
	}
	// dr5
	for(j=1; j<=NNBt; j++){
		for(i=1; i<NBz5; i++){
			dN[NBz0+NBz1+1+NBz2+NBz3+NBz4+i][j] = (zB6-zB5)/NBz5;
			dS[NBz0+NBz1+1+NBz2+NBz3+NBz4+i][j] = (zB6-zB5)/NBz5;
			dW[NBz0+NBz1+1+NBz2+NBz3+NBz4+i][j] = 2*Pi/NNBt;
			dE[NBz0+NBz1+1+NBz2+NBz3+NBz4+i][j] = 2*Pi/NNBt;
			dZ[NBz0+NBz1+1+NBz2+NBz3+NBz4+i][j] = (dN[NBz0+NBz1+1+NBz2+NBz3+NBz4+i][j]+dS[NBz0+NBz1+1+NBz2+NBz3+NBz4+i][j])/2.0;
			dT[NBz0+NBz1+1+NBz2+NBz3+NBz4+i][j] = (dW[NBz0+NBz1+1+NBz2+NBz3+NBz4+i][j]+dE[NBz0+NBz1+1+NBz2+NBz3+NBz4+i][j])/2.0;
		}
		dN[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j] = TINY;
		dS[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j] = (zB6-zB5)/NBz5;
		dW[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j] = 2*Pi/NNBt;
		dE[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j] = 2*Pi/NNBt;
		dZ[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j] = (dN[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j]+dS[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j])/2.0;
		dT[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j] = (dW[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j]+dE[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j])/2.0;
	}
	//// dr6
	//for(j=1; j<=NNBt; j++){
	//	for(i=1; i<NBz6; i++){
	//		dN[NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j] = (zB7-zB6)/NBz6;
	//		dS[NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j] = (zB7-zB6)/NBz6;
	//		dW[NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j] = 2*Pi/NNBt;
	//		dE[NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j] = 2*Pi/NNBt;
	//		dZ[NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j] = (dN[NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j]+dS[NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j])/2.0;
	//		dT[NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j] = (dW[NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j]+dE[NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j])/2.0;
	//	}
	//	dN[NNBz][j] = TINY;
	//	dS[NNBz][j] = (zB7-zB6)/NBz6;
	//	dW[NNBz][j] = 2*Pi/NNBt;
	//	dE[NNBz][j] = 2*Pi/NNBt;
	//	dZ[NNBz][j] = (dN[NNBz][j]+dS[NNBz][j])/2.0;
	//	dT[NNBz][j] = (dW[NNBz][j]+dE[NNBz][j])/2.0;
	//}

// Coordinate of each node
	for(j=1; j<=NNBt; j++){
		Zpc[1][j] = -Lpg/Rp;
		Zw[1][j] = Zpc[1][j];            
		Ze[1][j] = Zpc[1][j];           
		Zn[1][j] = Zpc[1][j] + dN[1][j]/2.0; 
		Zs[1][j] = Zpc[1][j] - dS[1][j]/2.0; 
	}
	for(j=1; j<=NNBt; j++){
		for(i=2; i<=1+NBz0+NBz1+NBz2+NBz3+NBz4+NBz5; i++){
			Zpc[i][j] = Zpc[i-1][j] + dS[i][j];	 // Radius of node[i][j]
			Zw[i][j] = Zpc[i][j];            
			Ze[i][j] = Zpc[i][j];           
			Zn[i][j] = Zpc[i][j] + dN[i][j]/2.0; 
			Zs[i][j] = Zpc[i][j] - dS[i][j]/2.0; 
		}
	}
	// Coordinate of each node
	for(i=1; i<=1+NBz0+NBz1+NBz2+NBz3+NBz4+NBz5; i++){
		T[i][1] = 0.0;
		Tw[i][1] = T[i][1] + dW[i][1]/2.0;
		Te[i][1] = T[i][1] - dE[i][1]/2.0;
		Tn[i][1] = T[i][1]; 
		Ts[i][1] = T[i][1]; 
	}
	for(j=2; j<=NNBt; j++){
		for(i=1; i<=1+NBz0+NBz1+NBz2+NBz3+NBz4+NBz5; i++){
			T[i][j] = T[i][j-1] + dW[i][j];	// Angle of node[i][j]
			Tw[i][j] = T[i][j] + dW[i][j]/2.0;
			Te[i][j] = T[i][j] - dE[i][j]/2.0;
			Tn[i][j] = T[i][j]; 
			Ts[i][j] = T[i][j]; 
		}
	}

	//printf("%e \n", Zpc[NNBz][1]);
}

/*	area of each control volume */
void area(float **Acv, float *Acvt, float **Zn, float **Zs, float **dT)
{
	int i, j;

	//	area of each control volume
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=1+NBz0+NBz1+NBz2+NBz3+NBz4+NBz5; i++){
			Acv[i][j] = (Zn[i][j]-Zs[i][j])*dT[i][j]; 
		}
		//Acv[1][j] = 0.5*(Zn[i][j]*Zn[i][j]-Zs[i][j]*Zs[i][j])*dT[i][j]*0.5;   // half-CV
		//Acv[NNBz][j] = 0.5*(Zn[i][j]*Zn[i][j]-Zs[i][j]*Zs[i][j])*dT[i][j]*0.5;  // half-CV
	}	

	//	Total Nondimentsional face area
	*Acvt = 0.0;		
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=1+NBz0+NBz1+NBz2+NBz3+NBz4+NBz5; i++){
			*Acvt = *Acvt + Acv[i][j];
		}
	}
	//printf("%e \n", *Acvt);
}

/* Ring Groove depth */
void groove(double **H_gro, float **Tpc)
{
	int i,j,k,l;
	int NBch, dNBch;

	double **H_g;

	H_g = dmatrix(1,NNBz,1,NNBt);

	// Initialize groove depth
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBz; i++){
			H_gro[i][j] = 0.0;
			H_g[i][j] = 0.0;
		}
	}
	// groove ring
	//for(j=1; j<=NNBt; j++){
	//	for(i=NBz0+3+NBz1; i<NBz0+NBz1+NBz2; i++){
	//		//H_g[i][j] = Hg;
	//		H_gro[i][j] = H_gro[i][j] + H_g[i][j];
	//	}
	//	H_gro[2+NBz0+NBz1][j] = Hg/2.0;
	//	H_gro[NBz0+NBz1+NBz2][j] = Hg/2.0;
	//	//H_gro[2+NBz1][j] = H_gro[2+NBz1][j] + H_g[2+NBz1][j];
	//	//H_gro[NBz1+NBz2][j] = H_gro[NBz1+NBz2][j] + H_g[NBz1+NBz2][j];
	//}

	// Profiling
	//for(j=1; j<=NNBt; j++){
	//	for(i=1; i<=1+NBz0; i++){
	//		H_gro[i][j] = 5.0e-6/hp*(1-sin(Pi/2*(i-1)/NBz0));
	//	}
	//}


	//// Nch Channels into the housing 
	//NBch = int(thetach/(2*Pi/NNBt)+0.5);
	//dNBch = int(NNBt/Nch+0.5);
	//for(i=1; i<=Nch; i++){
	//	for(k=NBz1+NBz2+NBz3+NBz4; k<=NNBz; k++){
	//		for(j=NNBt/12+1+dNBch*(i-1); j<=NNBt/12-1+NBch/2+dNBch*(i-1); j++){			
	//			H_gro[k][j] = Hg;
	//		}
	//		H_gro[k][NNBt/12+NBch/2+dNBch*(i-1)] = Hg/2.0;
	//		for(j=NNBt/12+3+dNBch*(i-1)-NBch/2; j<=NNBt/12+dNBch*(i-1); j++){				
	//			H_gro[k][j] = Hg;
	//		}
	//		H_gro[k][NNBt/12+2+dNBch*(i-1)-NBch/2] = Hg/2.0;
	//	}
	//}

	// pocket groove
	//int Ngr = 24;
	//float Hgr = 6.0e-6/hv0;
	//float thetagr = 5.0*Pi/180;
	//float drG = 1.6e-3/Rbo;
	//int NBgrdR = int(drG/(zB6-zB5)*NBz5+0.5);
	//int NBgrR1 = int((NBz5-NBgrdR)/2+0.5);
	//int NBgrR2 = NBgrR1 + NBgrdR;
	//int NBgrT = int(thetagr/(2*Pi/NNBt)+0.5);
	//int dNBgrT = int(NNBt/Ngr+0.5);
	//for(i=1; i<=Ngr; i++){
	//	for(k=1+NBz1+NBz2+NBz3+NBz4+NBgrR1+2; k<=1+NBz1+NBz2+NBz3+NBz4+NBgrR2-2; k++){
	//		for(j=1+dNBgrT*(i-1)+2; j<=1+NBgrT+dNBgrT*(i-1)-2; j++){
	//			H_gro[k][j] = Hgr;
	//		}
	//	}
	//	for(k=1+NBz1+NBz2+NBz3+NBz4+NBgrR1+1; k<=1+NBz1+NBz2+NBz3+NBz4+NBgrR2-1; k++){
	//		H_gro[k][1+dNBgrT*(i-1)+1] = Hgr/2;
	//		H_gro[k][1+dNBgrT*i-1] = Hgr/2;
	//	}
	//	for(j=1+dNBgrT*(i-1)+1; j<=1+NBgrT+dNBgrT*(i-1)-1; j++){
	//		H_gro[1+NBz1+NBz2+NBz3+NBz4+NBgrR1+1][j] = Hgr/2;
	//		H_gro[1+NBz1+NBz2+NBz3+NBz4+NBgrR2-1][j] = Hgr/2;
	//	}					
	//}


	// Wavy groove - CircumeZential
	//float Wh1 = 0.6e-6/hv0, Wh2 = 0.8e-6/hv0, Wh3 = 0.8e-6/hv0, Wh5 = 0.4e-6/hv0;
	//int Nwa1 = 30, Nwa2 = 38, Nwa3 = 46, Nwa5 = 54;
	//float thetawa1 = 6.0*Pi/180, thetawa2 = 5.0*Pi/180, thetawa3 = 4.0*Pi/180, thetawa5 = 3.0*Pi/180;
	//int NBwaT1 = int(thetawa1/(2*Pi/NNBt)+0.5);
	//int NBwaT2 = int(thetawa2/(2*Pi/NNBt)+0.5);
	//int NBwaT3 = int(thetawa3/(2*Pi/NNBt)+0.5);
	//int NBwaT5 = int(thetawa5/(2*Pi/NNBt)+0.5);
	//int dNBwaT1 = int(NNBt/Nwa1+0.5);
	//int dNBwaT2 = int(NNBt/Nwa2+0.5);
	//int dNBwaT3 = int(NNBt/Nwa3+0.5);
	//int dNBwaT5 = int(NNBt/Nwa5+0.5);

	//for(k=1; k<=Nwa1; k++){
	//	for(i=1; i<=1+NBz1; i++){
	//		for(j=1+dNBwaT1*(k-1); j<=1+NBwaT1+dNBwaT1*(k-1); j++){
	//			H_g[i][j] = Wh1*sin((j-dNBwaT1*(k-1)-1)*dNBt/thetawa1*2*Pi);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//}
	//for(k=1; k<=Nwa2; k++){
	//	for(i=2+NBz1; i<=1+NBz1+NBz2; i++){
	//		for(j=1+dNBwaT2*(k-1); j<=1+NBwaT2+dNBwaT2*(k-1); j++){
	//			H_g[i][j] = Wh2*sin((j-dNBwaT2*(k-1)-1)*dNBt/thetawa2*2*Pi);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//}
	//for(k=1; k<=Nwa3; k++){
	//	for(i=2+NBz1+NBz2; i<=1+NBz1+NBz2+NBz3; i++){
	//		for(j=1+dNBwaT3*(k-1); j<=1+NBwaT3+dNBwaT3*(k-1); j++){
	//			H_g[i][j] = Wh3*sin((j-dNBwaT3*(k-1)-1)*dNBt/thetawa3*2*Pi);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//}
	//for(k=1; k<=Nwa5; k++){
	//	for(i=2+NBz1+NBz2+NBz3+NBz4; i<=NNBz; i++){
	//		for(j=1+dNBwaT5*(k-1); j<=1+NBwaT5+dNBwaT5*(k-1); j++){
	//			H_g[i][j] = Wh5*sin((j-dNBwaT5*(k-1)-1)*dNBt/thetawa5*2*Pi);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//}

	// Wavy groove - Radial
	//float WhR1 = 1.6e-6/hv0, WhR2 = 0.4e-6/hv0, WhR3 = 2.0e-6/hv0, WhR5 = 0.6e-6/hv0;
	//int NwaR1 = 1, NwaR2 = 2, NwaR3 = 1, NwaR5 = 2;
	//float dZwa1=2*(zB2-zB1), dZwa2=0.5*(zB3-zB2), dZwa3=2*(zB4-zB3), dZwa5=0.5*(zB5-zB4);
	//int NBwaR1 = int(dZwa1/((zB2-zB1)/NBz1)+0.5);
	//int NBwaR2 = int(dZwa2/((zB3-zB2)/NBz2)+0.5);
	//int NBwaR3 = int(dZwa3/((zB4-zB3)/NBz3)+0.5);
	//int NBwaR5 = int(dZwa5/((zB6-zB5)/NBz5)+0.5);
	//int dNBwaR1 = int(NBz1/float(NwaR1)+0.5);
	//int dNBwaR2 = int(NBz2/float(NwaR2)+0.5);
	//int dNBwaR3 = int(NBz3/float(NwaR3)+0.5);
	//int dNBwaR5 = int(NBz5/float(NwaR5)+0.5);

	//for(j=1; j<=NNBt; j++){
	//	for(k=1; k<=NwaR1; k++){
	//		for(i=1+dNBwaR1*(k-1); i<=1+NBwaR1+dNBwaR1*(k-1); i++){
	//			H_g[i][j] = WhR1*sin(float(i-1)/NBwaR1*2.0*Pi*NwaR1);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//	for(k=1; k<=NwaR2; k++){
	//		for(i=1+NBz1+dNBwaR2*(k-1); i<=1+NBz1+NBwaR2+dNBwaR2*(k-1); i++){
	//			H_g[i][j] = WhR2*sin(float(i-(1+NBz1))/NBwaR2*2.0*Pi*NwaR2);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//	for(k=1; k<=NwaR3; k++){
	//		for(i=1+NBz1+NBz2+dNBwaR3*(k-1); i<=1+NBz1+NBz2+NBwaR3+dNBwaR3*(k-1); i++){
	//			H_g[i][j] = WhR3*sin(float(i-(1+NBz1+NBz2))/NBwaR3*2.0*Pi*NwaR3);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//	for(k=1; k<=NwaR5; k++){
	//		for(i=1+NBz1+NBz2+NBz3+NBz4+dNBwaR5*(k-1); i<=1+NBz1+NBz2+NBz3+NBz4+NBwaR5+dNBwaR5*(k-1); i++){
	//			H_g[i][j] = WhR5*sin(float(i-(1+NBz1+NBz2+NBz3+NBz4))/NBwaR5*2.0*Pi*NwaR5);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//}

//	// Wavy groove - 3D coning
//	float Wh1 = 1.2e-6/hv0, Wh2 = 1.2e-6/hv0, Wh3 = 1.2e-6/hv0, Wh5 = 0.0e-6/hv0, Rv;
//	int NwaC1 = 108, NwaC2 = 108, NwaC3 = 108, NwaC5 = 108;
//	float thetawa1 = 3.33333*Pi/180, thetawa2 = 3.33333*Pi/180, thetawa3 = 3.33333*Pi/180, thetawa5 = 3.33333*Pi/180;
//	int NBwaT1 = int(thetawa1/(2*Pi/NNBt)+0.5);
//	int NBwaT2 = int(thetawa2/(2*Pi/NNBt)+0.5);
//	int NBwaT3 = int(thetawa3/(2*Pi/NNBt)+0.5);
//	int NBwaT5 = int(thetawa5/(2*Pi/NNBt)+0.5);
//	int dNBwaT1 = int(NNBt/NwaC1+0.5);
//	int dNBwaT2 = int(NNBt/NwaC2+0.5);
//	int dNBwaT3 = int(NNBt/NwaC3+0.5);
//	int dNBwaT5 = int(NNBt/NwaC5+0.5);
//	//float WhR1 = -0.4e-6/hv0, WhR2 = -0.4e-6/hv0, WhR3 = -0.4e-6/hv0, WhR5 = -0.4e-6/hv0;
//	int NwaR1 = 1, NwaR2 = 1, NwaR3 = 1, NwaR5 = 1;
//	float dZwa1=(zB2-zB1)/NwaR1, dZwa2=(zB3-zB2)/NwaR2, dZwa3=(zB4-zB3)/NwaR3, dZwa5=(zB6-zB5)/NwaR5;
//	int NBwaR1 = int(dZwa1/((zB2-zB1)/NBz1)+0.5);
//	int NBwaR2 = int(dZwa2/((zB3-zB2)/NBz2)+0.5);
//	int NBwaR3 = int(dZwa3/((zB4-zB3)/NBz3)+0.5);
//	int NBwaR5 = int(dZwa5/((zB6-zB5)/NBz5)+0.5);
//	int dNBwaR1 = int(NBz1/float(NwaR1)+0.5);
//	int dNBwaR2 = int(NBz2/float(NwaR2)+0.5);
//	int dNBwaR3 = int(NBz3/float(NwaR3)+0.5);
//	int dNBwaR5 = int(NBz5/float(NwaR5)+0.5);
//
//	//Znf = (pow(Wh1*hv0,2)+Rbo*Rbo)/2/(Wh1*hv0)/hv0;
//	//Rv = 100000*zB6*Rbo;
//	//float A = (pow(Rbo,2)-pow(zB1*Rbo,2))/2/Wh1/hv0-Wh1*hv0/2;
//	//float Znf = sqrt(pow(A,2)+pow(Rbo,2))/hv0;
//	float Rc1, Rc3, Rc5, a1, a3, a5, b1, b3, b5, c1, c3, c5;
//	Rc1 = zB1;
//	a1 = Wh1/((zB2*zB2-zB1*zB1)-2*Rc1*(zB2-zB1));
//	b1 = -2*Rc1*a1;
//	c1 = -a1*zB1*zB1-b1*zB1;
//	Rc3 = zB4;
//	a3 = Wh3/((zB3*zB3-zB4*zB4)-2*Rc3*(zB3-zB4));
//	b3 = -2*Rc3*a3;
//	c3 = -a3*zB4*zB4-b3*zB4;
//	Rc5 = zB5;
//	a5 = Wh5/((zB6*zB6-zB5*zB5)-2*Rc5*(zB6-zB5));
//	b5 = -2*Rc5*a5;
//	c5 = -a5*zB5*zB5-b5*zB5;
//
//	for(i=1; i<1+NBz1; i++){
//		for(j=1; j<=NNBt; j++){
//			H_g[i][j] = a1*pow(zB1+(i-1)*((zB2-zB1)/NBz1),2)+b1*(zB1+(i-1)*((zB2-zB1)/NBz1))+c1;
//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 		
//		}
//	}
//
//	for(i=1+NBz1; i<1+NBz1+NBz2; i++){
//		for(j=1; j<=NNBt; j++){
//			H_g[i][j] = Wh2;
//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
//		}
//	}
//
//	for(i=1+NBz1+NBz2; i<=1+NBz1+NBz2+NBz3; i++){
//		for(j=1; j<=NNBt; j++){
//			H_g[i][j] = a3*pow(zB3+(i-1-NBz1-NBz2)*((zB4-zB3)/NBz3),2)+b3*(zB3+(i-1-NBz1-NBz2)*((zB4-zB3)/NBz3))+c3;
//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
//		}
//	}
//
//
//	for(i=1+NBz1+NBz2+NBz3+NBz4; i<=1+NBz1+NBz2+NBz3+NBz4+NBz5; i++){
//		for(j=1; j<=NNBt; j++){
//			H_g[i][j] = a5*pow(zB5+(i-1-NBz1-NBz2-NBz3-NBz4)*((zB6-zB5)/NBz5),2)+b5*(zB5+(i-1-NBz1-NBz2-NBz3-NBz4)*((zB6-zB5)/NBz5))+c5;
//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
//		}
//	}
//
//
//	free_dmatrix(H_g,1,NNBz,1,NNBt);
//
}