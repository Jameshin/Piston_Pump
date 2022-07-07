/* Calculate Cylinder pressure & pressure in the Valve and Dynamics of the Barrel */
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

double Avo(double );
float Vp, **dPdZ_n, **dPdT_c;
void Pgradient(float **Ppc, float **dPdZ_n, float **dPdT_c, float **dN, float **dS, float **dT);
void FrictionPC(double phi, float *Ffric, float *Mfricx, float *Mfricy, float **Zpc, float **Tpc, double **Hpc, float **Acv, float **dPdT_c, float **Ppc);

int Nd;
int Nphi, NBz0, NBz1, NBz2, NBz3, NBz4, NBz5, NBz6, NBz, NNBz, NNBt, Nch, Imin, Jmin;
int NVt1, NVt2, NVt3, NVt4, NVt5, NVt6, NVr1, NVr2, NVr3, NVr4;
double dNBt, thetach;
double R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, revf, phif, phiout_c, hm, phiout_cy, phiout_dy, phiout_vb, phiout_FM, phiout_LK, w;
double V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
double integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, zB1, zB2, zB3, zB4, zB5, zB6, zB7, wp, e0, theta0, psi0, de0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
double Lambda, Loc, Lpj, Lpg, Lbj, Ls, Rbo, Lpc;
float FxhD, FyhD, MxhD, MyhD, FxcD, FycD, MxcD, MycD, FxD, FyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, FzD, MxD, MyD, MzD, leakageRinD, leakageRoutD, frictionD, friclossD, Hmin;

void main()
{
	// 실린더 압력 관련 함수들
	void integBpre();
	void Init(double *);
	void Derivs(double , double *, double *,  float *);
	void RK4c(double *, double *, int , double *, double , double *, double *, double *, void (*Derivs)(double, double *, double *, float *), float *);
	void Leakage(double , double *, double *, double *, double *, double *, double *, double *);	
	// 유막압력 & 동역학 관련 함수들 
	void Input();
	void parameter();
	void grid(float **R, float **Zn, float **Zs, float **Zw, float **Ze, float **T, float **Tn, float **Ts, float **Tw, float **Te, float **dR, float **dT, float **dN, float **dS, float **dW, float **dE);
	void Upgrid(float **Zpc, float **Zn, float **Zs, float **Zw, float **Ze, float **T, float **Tn, float **Ts, float **Tw, float **Te, float **dR, float **dT, float **dN, float **dS, float **dW, float **dE, float **Acv, float *Acvt);
	void area(float **Acv, float *Acvt, float **Zn, float **Zs, float **dT);	
	void initial(float **Ppc, double *ym, double **Hpc, double **Hpc_gro, float **Zpc, float **Tpc, double *Pc, float *);
	void savedata(int *n, FILE *data0, FILE *data1, FILE *data2, FILE *,  FILE *,  FILE *, FILE *, double phi, double *hout_cy, double *hout_dy, double *hout_vb, double *hout_FM, double *hout_LK, double hout, double *Pc, float **Zpc, float **Tpc, float **Ppc, double *ym, float *, float **, double **);
	void FM(double  phi, double  *Pc, float **Zpc, float **Ppc, float **Zpcn, float **Zpcs, float **Tpc, float **dTpc, float **Acv, float *Fx, float *Fy, float *fr, float *Mx, float *My, float *Mz, double *ym, double **Hpc, float **dPdZ_n);
	void motion_RK4(double phi, double h, double *ym, double e, double kai, double gamma, double zeta, double de, double dkai,double dgamma, double dzeta, double **Hpc, double **Hn, double **Hs, double **Hw, double **He, double **Hpc_gro, 
				float **Zpc, float **Tpc, float **dTpc, float **dZpc, float **Zpcn, float **Zpcs, float **Zpcw, float **Zpce, float **Tpcn, float **Tpcs, float **Tpcw, float **Tpce, 
				float **dNvb, float **dSvb, float **dWvb, float **dEvb, float **aN, float **aS, float **aW, float **aE, float **aC, float **B, float **Ppc, double *Pcn, float **Acv,
				void (*configuration)(double , double *, double *, double *, double *,double *, double *, double *, double *, double *), 
				void (*filmthickness)(double , double *, double **, double **, double **, double **, double **, double **, float **, float **, double , double , double , double ), 
				void (*descretize)(double , double *, double , double , double , double , double , double , double , double , float **, float **, float **, double **, double **, double **, double **, float **, float **, 
									float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **), 
				void (*solver)(double, double *, float **, float **, float **, float **, float **, float **, float **), 
				void (*FM)(double , double *, float **, float **, float **, float **, float **, float **, float **, float *, float *, float *, float *, float *, float *, double *, double **Hpc, float **dPdZ_n) );
	void configuration(double phi, double *ym, double *e, double *kai, double *gamma, double *zeta, double *de, double *dkai, double *dgamma, double *dzeta);
	void filmthickness(double phi, double *ym, double **H, double **Hn, double **Hs, double **Hw, double **He, double **H_gro, float **R, float **T, double , double, double gamma, double zeta);
	void descretize(double phi, double *ym, double e, double kai, double gamma, double zeta, double de, double dkai, double dgamma, double dzeta, float **R, float **dT, float **dR, double **Hn, double **Hs, double **Hw, double **He, 
					float **Zn, float **Zs, float **Zw, float **Ze, float **Tn, float **Ts, float **Tw, float **Te, float **dN, float **dS, float **dW, float **dE, 
					float **aN, float **aS, float **aW, float **aE, float **aC, float **B);
	void solver(double phi, double *Pc, float **P, float **aN, float **aS, float **aW, float **aE, float **aC, float **B);
	void LeakageVB(float *leakage, float *PCleakage, float **Ppc, float **dN, float **dS, float **dT, float **dR, float **Zn, double **Hpc,double **Hn, float **, float **);
	//void Pgradient(float **Ppc, float **dPdZ_n, float **dPdT_c, float **dN, float **dS, float **dT);
	//void FrictionPC(double phi, float *Tb, float *Tv, float **Zpc, double **Hpc, float **Acv, float **dPdT_c, float **Ppc);

	int i, j, n=0;
	double phi=0.0, phiend, phit;
	double h, hout, hout_cy, hout_dy, hout_vb, hout_FM, hout_LK, hdid, hnext, hcal, hmcal;
	double *Pc, *Pcn, *dPc, *Pscal, Qvo; 
	double *ym, e=0.0, kai=0.0, gamma=0.0, zeta=0.0, de=0.0, dkai=0.0, dgamma=0.0, dzeta=0.0; 
	float **Zpc, **Zpcn, **Zpcs, **Zpcw, **Zpce, **Tpc, **Tpcn, **Tpcs, **Tpcw, **Tpce, **dNvb, **dSvb, **dWvb, **dEvb, **dZpc, **dTpc;	
	float **Acv, Acvt;
	double **Hpc, **Hn, **Hs, **Hw, **He, **Hpc_gro; 
	float **Ppc, **aN, **aS, **aW, **aE, **aC, **B;
	double Qp, Qv, Qs, Qpt, Qvt, Qst;
	float *leakage, *PCleakage, Tb, Tv;
	

	FILE *data0, *data1, *data2, *data3, *data4, *data5, *data6;

	if((data0 = fopen("Output_pressure2.dat", "w"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}
	if((data1 = fopen("pressure_VB.dat", "w"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}
	if((data2 = fopen("Dynamics.dat", "w"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}
	if((data3 = fopen("FThickness.dat", "w"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}
	if((data4 = fopen("ForceMoment.dat", "w"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}
	if((data5 = fopen("Leakage.dat", "w"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}
	if((data6 = fopen("eVB.dat", "w"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}
		
		
	Pc = dvector(1,NRK);
	Pcn = dvector(1,NRK);
	dPc = dvector(1,NRK);
	Pscal = dvector(1,NRK);	
	

	/*	PARAMETERS */
	Input();		// Read parameteZs 
	parameter();	// Unit conveZsion	
	
	ym = dvector(1,2*DOF);
	Zpc = matrix(1,NNBz,1,NNBt);	// Rad. coord. matrix
	Zpcn = matrix(1,NNBz,1,NNBt);	// North Rad. coord. matrix
	Zpcs = matrix(1,NNBz,1,NNBt);	// South Rad. 
	Zpcw = matrix(1,NNBz,1,NNBt);	// West Rad. 
	Zpce = matrix(1,NNBz,1,NNBt);	// East Rad. 
	Tpc = matrix(1,NNBz,1,NNBt);	// Circumf. coord. matrix
	Tpcn = matrix(1,NNBz,1,NNBt);	// North Circumf. coord. matrix
	Tpcs = matrix(1,NNBz,1,NNBt);	// South Circumf. 
	Tpcw = matrix(1,NNBz,1,NNBt);	// West Circumf. 
	Tpce = matrix(1,NNBz,1,NNBt);	// East Circumf. 
	dNvb = matrix(1,NNBz,1,NNBt);	// North Rad. coord. matrix
	dSvb = matrix(1,NNBz,1,NNBt);	// South Rad. 
	dWvb = matrix(1,NNBz,1,NNBt);	// West Rad. 
	dEvb = matrix(1,NNBz,1,NNBt);	// East Rad. 
	dZpc = matrix(1,NNBz,1,NNBt);	// delta Rad. CV 
	dTpc = matrix(1,NNBz,1,NNBt);	// delta Circumf. CV 
	Acv = matrix(1,NNBz,1,NNBt);
	Ppc = matrix(1,NNBz,1,NNBt);	// V-B pressure matrix
	Hpc = dmatrix(1,NNBz,1,NNBt);	// V-B Film thickness matrix
	Hn = dmatrix(1,NNBz,1,NNBt);
	Hs = dmatrix(1,NNBz,1,NNBt);
	Hw = dmatrix(1,NNBz,1,NNBt);
	He = dmatrix(1,NNBz,1,NNBt);
	Hpc_gro = dmatrix(1,NNBz,1,NNBt);	// V-B Groove matrix
	aN = matrix(1,NNBz,1,NNBt);			// 차분화 계수aN, aS, aW, aE, B
	aS = matrix(1,NNBz,1,NNBt);  
	aW = matrix(1,NNBz,1,NNBt);
	aE = matrix(1,NNBz,1,NNBt);
	aC = matrix(1,NNBz,1,NNBt);
	B = matrix(1,NNBz,1,NNBt);

	PCleakage = vector(1,N);
	leakage = vector(1,NNBz);

	dPdZ_n = matrix(1,NNBz,1,NNBt);
	dPdT_c = matrix(1,NNBz,1,NNBt);

	/*	I.PRE-PROCESS */
	initial(Ppc, ym, Hpc, Hpc_gro, Zpc, Tpc, Pc, PCleakage);		// Initial configuration & pressure		
	grid(Zpc, Zpcn, Zpcs, Zpcw, Zpce, Tpc, Tpcn, Tpcs, Tpcw, Tpce, dZpc, dTpc, dNvb, dSvb, dWvb, dEvb);		// Grid(Mesh) data	
	area(Acv, &Acvt, Zpcn, Zpcs, dTpc);					//	area of control volume & Total area

	/* Solve equations & Print the P_cylinder Results*/
	//fprintf(data0, "*Input File: Input_pressure_asiatribo.dat, *Output File: Output_pressure2.dat \n");
	fprintf(data0, "  Phi		P1		P2		P3	     P4		   P5		  P6		P7		P8		P9		Po		Pi \n");
	fprintf(data0, "%4.1f  ", phi*180/Pi);	
	for(i=1; i<=NRK; i++) fprintf(data0, "%e  ", Pc[i]);
	fprintf(data0,"\n");
	fprintf(data2, "  Phi		h0	     gamma		zeta		theta	 psi		dh0		wx		wy \n");
	fprintf(data4, " Phi	Balance	FxD		MyxD	MyyD	FyD		MxxD		MxyD		FlatxD		FlatyD		MlatxD		MlatyD		FfricD		 MfricxD		MfricyD		 FzD	MxD		MyD		MzD \n");
	fprintf(data5, "  Phi			C1			C2			C3		     C4			 C5			  C6			C7			C8			C9			Rin		Rout	TotalLK		FricTorque		LeakageLoss		FricLoss	TotalLoss \n");
	
	// Test
	//fprintf(data1, "NCYCLE=0\n VARIABLES=R,T,pressure(Pa)\n ZONE T=P_Film, I=         %d ,J=        %d , F=POINT\n", NNBz, 720);
	//for(j=1; j<NNBt;j+=NNBt/720){
	//	for(i=1; i<=NNBz; i++){
	//		fprintf(data1, "%e %e %e \n", Zpc[i][j]*Rp, Tpc[i][j], Acv[i][j]*Rp*Rp);
	//	}
	///*fprintf(data1, "%e %e %e \n", Zpc[i][NNBt]*Rp, Tpc[i][NNBt], Acv[i][NNBt]*Rp*Rp);*/
	//}

	h = 0.00001*Pi/180;
	hout = 0.0;
	hout_cy = 0.0;
	hout_dy = 0.0;
	hout_vb = 0.0;
	hout_FM = 0.0;
	hout_LK = 0.0;
	Nd = 0;
	while(phi<phif){
		if(phi>=4.25*Pi){
			//Qvo = 0.0;
			exit(1);
			if(hout>=phiout_c){								
				/* Cylinder pressure data */
				//printf("\n%f\n", phi*180/Pi);
				//fprintf(fp,"%4.1f  ", phi*180/Pi);
				//for(i=1; i<=NRK; i++) fprintf(fp, "%e  ", Pc[i]);
				//Qvo = Pc[N+3];										// Outlet flow
				//fprintf(fp,"%e ", Qvo);
				//fprintf(fp,"%e ", Avo(phi));
				//Leakage(phi, Pc, &Qp, &Qv, &Qs, &Qpt, &Qvt, &Qst);	// Leakage in a Cylinder & Total leakage								// Zesultant force action point
				//fprintf(fp,"%e %e %e %e %e %e %e ", Qp, Qv, Qs, Qpt, Qvt, Qst, Qpt+Qvt+Qst);
				//fprintf(fp,"%e %e ", eBVx, eBVy);		
				
				phit = phi;
				while(phit > 2*Pi) phit = phit - 2*Pi;

				Lpc = Lc/Rp - Ld/Rp - R/Rp*tan(beta)*(1+cos(phit));
				Vp = -R/Rp*(tan(beta)*sin(phit));
				//Ap = -R/Rp*(tan(beta)*cos(phit));
				//printf("##### %e %e \n", phi, hout);
				for(i=1; i<=N+2; i++) Pcn[i] = Pc[i]/Lambda;		// Ze-define nondimensional parameteZs in cylinder pressure calculation		
			/* Updating */
				Upgrid(Zpc, Zpcn, Zpcs, Zpcw, Zpce, Tpc, Tpcn, Tpcs, Tpcw, Tpce, dZpc, dTpc, dNvb, dSvb, dWvb, dEvb, Acv, &Acvt);		// Grid(Mesh) data	
				for(j=1; j<=NNBt; j++){
					Ppc[NNBz][j] = Ph/Lambda;
					Ppc[1][j] = Pcn[1];
				}
			/*	II.FORCE & MOMENT */						
				hcal = 0.0;
				hmcal = hm;
				while(hcal<hout){
					if(hout-hcal < hm*180/Pi) hmcal = (hout-hcal)*Pi/180;
			/*	III.BODY DYNAMICS */
					printf("\n --------- %f %f \n",phit*180/Pi, hout);
					motion_RK4(phit, hmcal, ym, e, kai, gamma, zeta, de,  dkai, dgamma, dzeta, Hpc, Hn, Hs, Hw, He, Hpc_gro, 
								Zpc, Tpc, dTpc, dZpc, Zpcn, Zpcs, Zpcw, Zpce, Tpcn, Tpcs, Tpcw, Tpce, 
								dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Ppc, Pcn, Acv, 
								configuration, filmthickness, descretize, solver, FM);
					phit += hmcal;		// Motion angle step
					hcal += hmcal*180/Pi;
				}
				//printf("################################################################");
			/*	IV.PRESSURE DISTRIBUTION */	
				//configuration(phit, ym, &gamma, &zeta, &dgamma, &dzeta);
				//filmthickness(phit, ym, Hpc, Hn, Hs, Hw, He, Hpc_gro, Zpc, Tpc, gamma, zeta);
				//descretize(phit, ym, gamma, zeta, dgamma, dzeta, Zpc, dTpc, dZpc, Hn, Hs, Hw, He, Zpcn, Zpcs, Zpcw, Zpce, 
				//			Tpcn, Tpcs, Tpcw, Tpce, dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B);
				//solver(phit, Pcn, Ppc, aN, aS, aW, aE, aC, B);
			/*  V. Leakage & Friction */
				//Pgradient(Ppc, dPdZ_n, dPdT_c, dNvb, dSvb, dTpc);
				//LeakageVB(leakage, PCleakage, Ppc, dNvb, dSvb, dTpc, dZpc, Zpcn, Hpc, Hn, dPdZ_n, dPdT_c);
				//FrictionPC(phit, &Tb, &Tv, Zpc, Hpc, Acv, dPdT_c, Ppc);
			/*  VI. SAVE Cylinder pressure, Film pressure, and Motion data */
				savedata(&n, data0, data1, data2, data3, data4, data5, data6, phi-2*Pi, &hout_cy, &hout_dy, &hout_vb, &hout_FM, &hout_LK, hout, Pc, Zpc, Tpc, Ppc, ym, PCleakage, Acv, Hpc);	
				hout = 0.0;
				Nd++;
			}
			hout += h*180/Pi;
		}
		Derivs(phi, Pc, dPc, PCleakage);
		for(i=1; i<=NRK; i++) Pscal[i] = fabs(Pc[i]) + fabs(h*dPc[i]) + TINY;
		RK4c(Pc, dPc, NRK, &phi, h, Pscal, &hdid, &hnext, Derivs, PCleakage);
		for(i=1; i<=N+2; i++) if(Pc[i] < Pca/Lambda) Pc[i] = Pca/Lambda;
		printf("%e \n", phi*180/Pi);
		h = hnext;

		hout_cy += hout;		
		// Cylinder pressure
		if(hout_cy >= phiout_cy){
			fprintf(data0,"%4.1f  ", phi*180/Pi);
			for(i=1; i<=NRK; i++) fprintf(data0, "%e  ", Pc[i]);
			fprintf(data0,"\n");
			hout_cy = 0.0;
		}
	}	

	free_dvector(Pc,1,NRK);
	free_dvector(Pcn,1,NRK);
	free_dvector(dPc,1,NRK);
	free_dvector(Pscal,1,NRK);

	free_dvector(ym,1,2*DOF);
	free_matrix(Zpc,1,NNBz,1,NNBt);  // Rad. coord. matrix
	free_matrix(Zpcn,1,NNBz,1,NNBt);  // North Rad. coord. matrix
	free_matrix(Zpcs,1,NNBz,1,NNBt);  // South Rad. 
	free_matrix(Zpcw,1,NNBz,1,NNBt);  // West Rad. 
	free_matrix(Zpce,1,NNBz,1,NNBt);  // East Rad. 
	free_matrix(Tpc,1,NNBz,1,NNBt);  // Circumf. coord. matrix
	free_matrix(Tpcn,1,NNBz,1,NNBt);  // North Circumf. coord. matrix
	free_matrix(Tpcs,1,NNBz,1,NNBt);  // South Circumf. 
	free_matrix(Tpcw,1,NNBz,1,NNBt);  // West Circumf. 
	free_matrix(Tpce,1,NNBz,1,NNBt);  // East Circumf. 
	free_matrix(dNvb,1,NNBz,1,NNBt);  // North Rad. coord. matrix
	free_matrix(dSvb,1,NNBz,1,NNBt);  // South Rad. 
	free_matrix(dWvb,1,NNBz,1,NNBt);  // West Rad. 
	free_matrix(dEvb,1,NNBz,1,NNBt);  // East Rad. 
	free_matrix(dZpc,1,NNBz,1,NNBt);  // delta Rad. CV 
	free_matrix(dTpc,1,NNBz,1,NNBt);  // delta Circumf. CV 
	free_matrix(Acv,1,NNBz,1,NNBt);
	free_matrix(Ppc,1,NNBz,1,NNBt);  // V-B pressure matrix
	free_dmatrix(Hpc,1,NNBz,1,NNBt);  // V-B Film thickness matrix
	free_dmatrix(Hn,1,NNBz,1,NNBt);
	free_dmatrix(Hs,1,NNBz,1,NNBt);
	free_dmatrix(Hw,1,NNBz,1,NNBt);
	free_dmatrix(He,1,NNBz,1,NNBt);
	free_dmatrix(Hpc_gro,1,NNBz,1,NNBt);  // V-B Groove matrix
	free_matrix(aN,1,NNBz,1,NNBt);  // 차분화 계수aN, aS, aW, aE, B
	free_matrix(aS,1,NNBz,1,NNBt);  
	free_matrix(aW,1,NNBz,1,NNBt);
	free_matrix(aE,1,NNBz,1,NNBt);
	free_matrix(aC,1,NNBz,1,NNBt);
	free_matrix(B,1,NNBz,1,NNBt);
	
	free_vector(leakage,1,NNBz);
	free_vector(PCleakage,1,N);

	free_matrix(dPdZ_n,1,NNBz,1,NNBt);
	free_matrix(dPdT_c,1,NNBz,1,NNBt);

	fclose(data0);
	fclose(data1);
	fclose(data2);
	fclose(data3);
	fclose(data4);
	fclose(data5);
	fclose(data6);
}
/* Save output data */
void savedata(int *n, FILE *data0, FILE *data1, FILE *data2, FILE *data3, FILE *data4, FILE *data5, FILE *data6, double phi, double *hout_cy, double *hout_dy, double *hout_vb, double *hout_FM, double *hout_LK, double hout, double  *Pc, float **Zpc, float **Tpc, float **Ppc, double *ym, float *PCleakage, float **Acv, double **Hpc)
{
	int i, j, k;
	double phit;
	float eVBx, eVBy, Pe, eBVx, eBVy;

	//void Point_resul(double , double *, float *, float *);
	
	*hout_cy += hout;
	*hout_dy += hout;
	*hout_vb += hout;

	phit = phi;
	while(phit > 2*Pi) phit = phit - 2*Pi;

	//printf("################################################################");
	// Cylinder pressure
	if(*hout_cy >= phiout_cy){
		fprintf(data0,"%4.1f  ", phi*180/Pi);
		for(i=1; i<=NRK; i++) fprintf(data0, "%e  ", Pc[i]);
		fprintf(data0,"\n");
		*hout_cy = 0.0;
	}

	 //pressure_VB & Film Thickness 
	if(*hout_vb >= phiout_vb){
		fprintf(data1, "NCYCLE=%d\n VARIABLES=R,T,pressure(Pa)\n ZONE T=P_Film, I=         %d ,J=        %d , F=POINT\n", *n, NNBz, NNBt);
		fprintf(data3, "NCYCLE=%d\n VARIABLES=R,T,FThickness(um)\n ZONE T=P_Film, I=         %d ,J=        %d , F=POINT\n", *n, NNBz, NNBt);
		//fprintf(data3, "NCYCLE=%d\n VARIABLES=eVBx, eVBy\n ZONE T=P_Film, I=         1 ,J=        1 , F=POINT\n", *n);
		for(j=1; j<=NNBt; j++){
			//for(i=1; i<=NNBz; i++){ 
			//	//while(Tpc[i][j]+phit >= 2*Pi) Tpc[i][j] = Tpc[i][j] - 2*Pi;
			//	fprintf(data1, "%e %e %e \n", Zpc[i][j]*Rp, Pi/2-(Tpc[i][j]+phit), Ppc[i][j]*Lambda);	
			//	//fprintf(data3, "%e %e %e \n", Zpc[i][j]*Rp, Pi/2-(Tpc[i][j]+phit), (ym[1]*hp-Zpc[i][j]*Rp*acos(cos(ym[2])*cos(ym[3]))*cos(Tpc[i][j]-Pi/2-atan2(tan(ym[2]),sin(ym[3]))+phit))*1000000);
			//}
			for(i=1; i<1+NBz1; i++/*=NBz1/NBz1*/){ 
				fprintf(data1, "%e %e %e \n", Zpc[i][j]*Rp, Tpc[i][j], Ppc[i][j]*Lambda);
				fprintf(data3, "%e %e %e \n", Zpc[i][j]*Rp, Tpc[i][j], Hpc[i][j]*hp);				
			}
			for(i=1+NBz1; i<1+NBz1+NBz2; i++){ 
				fprintf(data1, "%e %e %e \n", Zpc[i][j]*Rp, Tpc[i][j], Ppc[i][j]*Lambda);				
				fprintf(data3, "%e %e %e \n", Zpc[i][j]*Rp, Tpc[i][j], Hpc[i][j]*hp);	
			}
			for(i=1+NBz1+NBz2; i<1+NBz1+NBz2+NBz3; i++/*=NBz3/NBz3*/){ 
				fprintf(data1, "%e %e %e \n", Zpc[i][j]*Rp, Tpc[i][j], Ppc[i][j]*Lambda);
				fprintf(data3, "%e %e %e \n", Zpc[i][j]*Rp, Tpc[i][j], Hpc[i][j]*hp);	
			}
			for(i=1+NBz1+NBz2+NBz3; i<1+NBz1+NBz2+NBz3+NBz4; i++/*=NBz4/4*/){ 
				fprintf(data1, "%e %e %e \n", Zpc[i][j]*Rp, Tpc[i][j], Ppc[i][j]*Lambda);
				fprintf(data3, "%e %e %e \n", Zpc[i][j]*Rp, Tpc[i][j], Hpc[i][j]*hp);	
			}
			for(i=1+NBz1+NBz2+NBz3+NBz4; i<=NNBz; i++/*=NBz5/NBz5*/){ 
				fprintf(data1, "%e %e %e \n", Zpc[i][j]*Rp, Tpc[i][j], Ppc[i][j]*Lambda);
				fprintf(data3, "%e %e %e \n", Zpc[i][j]*Rp, Tpc[i][j], Hpc[i][j]*hp);	
			}
		}
		for(i=1; i<1+NBz1; i++){ 
			fprintf(data1, "%e %e %e \n", Zpc[i][1]*Rp, Tpc[i][1], Ppc[i][1]*Lambda);
			fprintf(data3, "%e %e %e \n", Zpc[i][1]*Rp, Tpc[i][1], Hpc[i][1]*hp);	
		}
		for(i=1+NBz1; i<1+NBz1+NBz2; i++){ 
			fprintf(data1, "%e %e %e \n", Zpc[i][1]*Rp, Tpc[i][1], Ppc[i][1]*Lambda);
			fprintf(data3, "%e %e %e \n", Zpc[i][1]*Rp, Tpc[i][1], Hpc[i][1]*hp);	
		}
		for(i=1+NBz1+NBz2; i<1+NBz1+NBz2+NBz3; i++/*=NBz3/NBz3*/){ 
			fprintf(data1, "%e %e %e \n", Zpc[i][1]*Rp, Tpc[i][1], Ppc[i][1]*Lambda);
			fprintf(data3, "%e %e %e \n", Zpc[i][1]*Rp, Tpc[i][1], Hpc[i][1]*hp);	
		}
		for(i=1+NBz1+NBz2+NBz3; i<1+NBz1+NBz2+NBz3+NBz4; i++/*=NBz4/4*/){ 
			fprintf(data1, "%e %e %e \n", Zpc[i][1]*Rp, Tpc[i][1], Ppc[i][1]*Lambda);
			fprintf(data3, "%e %e %e \n", Zpc[i][1]*Rp, Tpc[i][1], Hpc[i][1]*hp);	
		}
		for(i=1+NBz1+NBz2+NBz3+NBz4; i<=NNBz; i++/*=NBz5/NBz5*/){ 
			fprintf(data1, "%e %e %e \n", Zpc[i][1]*Rp, Tpc[i][1], Ppc[i][1]*Lambda);
			fprintf(data3, "%e %e %e \n", Zpc[i][1]*Rp, Tpc[i][1], Hpc[i][1]*hp);	
		}
		fprintf(data1,"\n");
		fprintf(data3,"\n");

		*hout_vb = 0.0;
		*n = *n + 1;
	}		
	
	// Barrel Motion
	if(*hout_dy >= phiout_dy){				
		//printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
		//int Nphit1 = NNBt/4-int(phit/dNBt);
		//int Nphit2 = NNBt*7/12-int(phit/dNBt);
		//int Nphit3 = NNBt*11/12-int(phit/dNBt);
		//while(Nphit1 <= 0) Nphit1 = Nphit1 + NNBt;
		//while(Nphit2 <= 0) Nphit2 = Nphit2 + NNBt;
		//while(Nphit3 <= 0) Nphit3 = Nphit3 + NNBt;
		//while((Hpc[NNBz][Nphit1]*hp)*1e6 > 500) Hpc[NNBz][Nphit1] = Hpc[NNBz][Nphit1] - 1000e-6/hp;
		//while((Hpc[NNBz][Nphit2]*hp)*1e6 > 500) Hpc[NNBz][Nphit2] = Hpc[NNBz][Nphit2] - 1000e-6/hp;
		//while((Hpc[NNBz][Nphit3]*hp)*1e6 > 500) Hpc[NNBz][Nphit3] = Hpc[NNBz][Nphit3] - 1000e-6/hp;
		fprintf(data2, "\nphi = %e \n %4.2e %4.2e %4.2e %4.2e %6.5e %4.3e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %e %d %d \n\n", phi*180/Pi, sqrt(ym[1]*ym[1]+ym[2]*ym[2])*hp, atan2(ym[2],ym[1])*180/Pi, (acos(cos(ym[3])*cos(ym[4]))*Rp/hp)*hp/Rp*180/Pi, (-atan2(tan(ym[3]),sin(ym[4])))*180/Pi, ym[1]*hp, ym[2]*hp, ym[3]*180/Pi, ym[4]*180/Pi, ym[5]*hp*w, ym[6]*hp*w, ym[7]*w*180/Pi, ym[8]*w*180/Pi, (/*ym[1]-(*gamma)*/Hmin)*hp, Imin, Jmin);
		*hout_dy = 0.0;

		//eVBx = 0.0, eVBy = 0.0, Pe =0.0;
		//Point_resul(phi, Pc, &eBVx, &eBVy);	
		//for(j=1; j<=NNBt; j+=NNBt/540){
		//	for(i=1; i<=NNBz; i++){ 
		//		eVBx = eVBx + Ppc[i][j]*Lambda*Zpc[i][j]*Rp*cos(Pi/2-(Tpc[i][j]+phit))*Acv[i][j];
		//		eVBy = eVBy + Ppc[i][j]*Lambda*Zpc[i][j]*Rp*sin(Pi/2-(Tpc[i][j]+phit))*Acv[i][j];
		//		Pe = Pe + Ppc[i][j]*Lambda*Acv[i][j];
		//	}
		//}	
		//eVBx = eVBx/Pe;
		//eVBy = eVBy/Pe;
		//fprintf(data6, "%e %e %e %e \n", eVBx, eVBy, eBVx, eBVy);
	}

	// Force & Moment on the Barrel
	*hout_FM += hout;
	phiout_FM = 0.1;
	// Barrel Motion
	if(*hout_FM >= phiout_FM){				
		//printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
		fprintf(data4, "%4.1f %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",	phi*180/Pi,	FyD/FxD, FxD, MxcD, MycD, FyD, MxhD, MxhD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, FzD, MxD, MyD, MzD);
		*hout_FM = 0.0;
	}

	// Leakage & Friction loss
	//*hout_LK += hout;
	//phiout_LK = 1;
	//if(*hout_LK >= phiout_LK){	
	//	fprintf(data5, "%4.1f  ", phi*180/Pi);	
	//	for(i=1; i<=N; i++) fprintf(data5, "%e  ", PCleakage[i]);
	//	fprintf(data5, "%e  %e  %e  %e  %e  %e  %e  %e %e \n", -leakageRinD, leakageRoutD, -leakageRinD+leakageRoutD, frictionD, (Pout-Pin)*(-leakageRinD+leakageRoutD), friclossD, friclossD + (Pout-Pin)*(-leakageRinD+leakageRoutD), frictionD*w, frictionD*w + (Pout-Pin)*(-leakageRinD+leakageRoutD));
	//	*hout_LK = 0.0;
	//}

	//// Test - Boundary condition
	//fprintf(data1, "NCYCLE=%d\n VARIABLES=R,T,pressure(Pa)\n ZONE T=P_Film, I=         %d ,J=        %d , F=POINT\n", *n, NNBz, NNBt+1);
	//for(j=1; j<NNBt; j++){
	//	for(i=1; i<=NNBz; i++){ 				
	//		fprintf(data1, "%e %e %e \n", Zpc[i][j]*Rp, Tpc[i][j], Acv[i][j]*Rp*Rp);
	//	}
	//}
	//for(i=1; i<=NNBz; i++){ 
	//	fprintf(data1, "%e %e %e \n", Zpc[i][NNBt]*Rp, Tpc[i][NNBt], Acv[i][NNBt]*Rp*Rp);
	//}
	//fprintf(data1,"\n");

}