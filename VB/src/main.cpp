/* Calculate Cylinder pressure & Pressure in the Valve and Dynamics of the Barrel */
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

double Avo(double );
double Avi(double );

int Nd;
int Nphi, NBr1, NBr2, NBr3, NBr4, NBr5, NBr, NNBr, NNBt, Nch, Imin, Jmin;
int NVt1, NVt2, NVt3, NVt4, NVt5, NVt6, NVr1, NVr2, NVr3, NVr4;
float dNBt, thetach;
float R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, revf, phif, phiout_c, hm, phiout_cy, phiout_dy, phiout_vb, phiout_FM, phiout_LK, w;
float V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
float integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, rB1, rB2, rB3, rB4, rB5, rB6, Fs, hv0, theta0, psi0, dhv0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
float Lambda, Loc, Lpj, Lpg, Lbj, Ls, Rbo;
float FpressD, MpressxD, MpressyD, FliftD, MliftxD, MliftyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, FzD, MxD, MyD, MzD, leakageRinD, leakageRoutD, frictorqueD, frictorqueD2, friclossD, Hmin;

int main()
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
	void grid(float **R, float **Rn, float **Rs, float **Rw, float **Re, float **T, float **Tn, float **Ts, float **Tw, float **Te, float **dR, float **dT, float **dN, float **dS, float **dW, float **dE);
	void area(float **Acv, float *Acvt, float **Rn, float **Rs, float **dT);	
	void initial(float **Pvb, double *ym, double **Hvb, double **Hvb_gro, float **Rvb, float **Tvb, double *Pc, float *);
	void savedata(int *n, FILE *data0, FILE *data1, FILE *data2, FILE *,  FILE *,  FILE *, FILE *, double phi, double *hout_cy, double *hout_dy, double *hout_vb, double *hout_FM, double *hout_LK, double hout, double *Pc, float **Rvb, float **Tvb, float **Pvb, double *ym, float *, float **, double **);
	void FM(double  phi, double  *Pc, float **Rvb, float **Pvb, float **Rvbn, float **Rvbs, float **Tvb, float **dTvb, float **Acv, float *Fx, float *Fy, float *Fz, float *Mx, float *My, float *Mz, double *ym);
	void motion_RK4(double phi, double h, double *ym, double gamma, double zeta, double dgamma, double dzeta, double **Hvb, double **Hn, double **Hs, double **Hw, double **He, double **Hvb_gro, 
				float **Rvb, float **Tvb, float **dTvb, float **dRvb, float **Rvbn, float **Rvbs, float **Rvbw, float **Rvbe, float **Tvbn, float **Tvbs, float **Tvbw, float **Tvbe, 
				float **dNvb, float **dSvb, float **dWvb, float **dEvb, float **aN, float **aS, float **aW, float **aE, float **aC, float **B, float **Pvb, double *Pcn, float **Acv,
				void (*configuration)(double , double *, double *, double *, double *, double *), 
				void (*filmthickness)(double , double *, double **, double **, double **, double **, double **, double **, float **, float **, double , double ), 
				void (*descretize)(double , double *, double , double , double , double , float **, float **, float **, double **, double **, double **, double **, float **, float **, 
									float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **), 
				void (*solver)(double, double *, float **, float **, float **, float **, float **, float **, float **), 
				void (*FM)(double , double *, float **, float **, float **, float **, float **, float **, float **, float *, float *, float *, float *, float *, float *, double *) );
	void configuration(double phi, double *ym, double *gamma, double *zeta, double *dgamma, double *dzeta);
	void filmthickness(double phi, double *ym, double **H, double **Hn, double **Hs, double **Hw, double **He, double **H_gro, float **R, float **T, double gamma, double zeta);
	void descretize(double phi, double *ym, double gamma, double zeta, double dgamma, double dzeta, float **R, float **dT, float **dR, double **Hn, double **Hs, double **Hw, double **He, 
					float **Rn, float **Rs, float **Rw, float **Re, float **Tn, float **Ts, float **Tw, float **Te, float **dN, float **dS, float **dW, float **dE, 
					float **aN, float **aS, float **aW, float **aE, float **aC, float **B);
	void solver(double phi, double *Pc, float **P, float **aN, float **aS, float **aW, float **aE, float **aC, float **B);
	void LeakageVB(float *leakage, float *Baleakage, float **Pvb, float **dN, float **dS, float **dT, float **dR, float **Rn, double **Hvb,double **Hn, float **, float **);
	void Pgradient(float **Pvb, float **dPdR_n, float **dPdT_c, float **dN, float **dS, float **dT);
	void FrictionVB(double phi, float *Tb, float *Tv, float **Rvb, double **Hvb, float **Acv, float **dPdT_c, float **Pvb);

	int i, j, n=0;
	double phi=0.0, phiend, phit, *ph;
	double h, hout, hout_cy, hout_dy, hout_vb, hout_FM, hout_LK, hdid, hnext, hcal, hmcal;
	double *Pc, *Pcn, *dPc, *Pscal, Qvo, Qvi; 
	double *ym, gamma=0.0, zeta=0.0, dgamma=0.0, dzeta=0.0; 
	float **Rvb, **Rvbn, **Rvbs, **Rvbw, **Rvbe, **Tvb, **Tvbn, **Tvbs, **Tvbw, **Tvbe, **dNvb, **dSvb, **dWvb, **dEvb, **dRvb, **dTvb;	
	float **Acv, Acvt;
	double **Hvb, **Hn, **Hs, **Hw, **He, **Hvb_gro; 
	float **Pvb, **aN, **aS, **aW, **aE, **aC, **B;
	//double Qp, Qv, Qs, Qpt, Qvt, Qst;
	float *leakage, *Baleakage, Tb, Tv;
	float **dPdR_n, **dPdT_c;

	FILE *data0, *data1, *data2, *data3, *data4, *data5, *data6;

	if((data0 = fopen("Output_pressure2.dat", "w"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}
	if((data1 = fopen("Pressure_VB.dat", "w"))==NULL) {
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
		
	ph = dvector(1,N);
	Pc = dvector(1,NRK);
	Pcn = dvector(1,NRK);
	dPc = dvector(1,NRK);
	Pscal = dvector(1,NRK);	
	

	/*	PARAMETERS */
	Input();		// Read parameters 
	parameter();	// Unit conversion	
	
	
	ym = dvector(1,2*DOF);
	Rvb = matrix(1,NNBr,1,NNBt);	// Rad. coord. matrix
	Rvbn = matrix(1,NNBr,1,NNBt);	// North Rad. coord. matrix
	Rvbs = matrix(1,NNBr,1,NNBt);	// South Rad. 
	Rvbw = matrix(1,NNBr,1,NNBt);	// West Rad. 
	Rvbe = matrix(1,NNBr,1,NNBt);	// East Rad. 
	Tvb = matrix(1,NNBr,1,NNBt);	// Circumf. coord. matrix
	Tvbn = matrix(1,NNBr,1,NNBt);	// North Circumf. coord. matrix
	Tvbs = matrix(1,NNBr,1,NNBt);	// South Circumf. 
	Tvbw = matrix(1,NNBr,1,NNBt);	// West Circumf. 
	Tvbe = matrix(1,NNBr,1,NNBt);	// East Circumf. 
	dNvb = matrix(1,NNBr,1,NNBt);	// North Rad. coord. matrix
	dSvb = matrix(1,NNBr,1,NNBt);	// South Rad. 
	dWvb = matrix(1,NNBr,1,NNBt);	// West Rad. 
	dEvb = matrix(1,NNBr,1,NNBt);	// East Rad. 
	dRvb = matrix(1,NNBr,1,NNBt);	// delta Rad. CV 
	dTvb = matrix(1,NNBr,1,NNBt);	// delta Circumf. CV 
	Acv = matrix(1,NNBr,1,NNBt);
	Pvb = matrix(1,NNBr,1,NNBt);	// V-B Pressure matrix
	Hvb = dmatrix(1,NNBr,1,NNBt);	// V-B Film thickness matrix
	Hn = dmatrix(1,NNBr,1,NNBt);
	Hs = dmatrix(1,NNBr,1,NNBt);
	Hw = dmatrix(1,NNBr,1,NNBt);
	He = dmatrix(1,NNBr,1,NNBt);
	Hvb_gro = dmatrix(1,NNBr,1,NNBt);	// V-B Groove matrix
	aN = matrix(1,NNBr,1,NNBt);			// 차분화 계수aN, aS, aW, aE, B
	aS = matrix(1,NNBr,1,NNBt);  
	aW = matrix(1,NNBr,1,NNBt);
	aE = matrix(1,NNBr,1,NNBt);
	aC = matrix(1,NNBr,1,NNBt);
	B = matrix(1,NNBr,1,NNBt);

	Baleakage = vector(1,N);
	leakage = vector(1,NNBr);

	dPdR_n = matrix(1,NNBr,1,NNBt);
	dPdT_c = matrix(1,NNBr,1,NNBt);

	/*	I.PRE-PROCESS */
	grid(Rvb, Rvbn, Rvbs, Rvbw, Rvbe, Tvb, Tvbn, Tvbs, Tvbw, Tvbe, dRvb, dTvb, dNvb, dSvb, dWvb, dEvb);		// Grid(Mesh) data	
	area(Acv, &Acvt, Rvbn, Rvbs, dTvb);					//	Area of control volume & Total area
	initial(Pvb, ym, Hvb, Hvb_gro, Rvb, Tvb, Pc, Baleakage);		// Initial configuration & pressure		

	/* Solve equations & Print the P_cylinder results*/
	//fprintf(data0, "*Input File: Input_pressure_asiatribo.dat, *Output File: Output_pressure2.dat \n");
	fprintf(data0, "  Phi		P1		P2		P3	     P4		   P5		  P6		P7		P8		P9		Po		Pi \n");
	fprintf(data0, "%4.1f  ", phi*180/Pi);	
	for(i=1; i<=NRK; i++) fprintf(data0, "%e  ", Pc[i]);
	fprintf(data0,"\n");
	fprintf(data2, "  Phi		h0	     gamma		zeta		theta	 psi		dh0		wx		wy \n");
	fprintf(data4, " Phi	Balance	FpressD		MpressxD	MpressyD	FliftD		MliftxD		MliftyD		FlatxD		FlatyD		MlatxD		MlatyD		FfricD		 MfricxD		MfricyD		 FzD	MxD		MyD		MzD \n");
	fprintf(data5, "  Phi			C1			C2			C3		     C4			 C5			  C6			C7			C8			C9			Rin		Rout	TotalLK		FricTorque		LeakageLoss		FricLoss	TotalLoss \n");
	
	// Test
	//fprintf(data1, "NCYCLE=0\n VARIABLES=R,T,PRESSURE(Pa)\n ZONE T=P_Film, I=         %d ,J=        %d , F=POINT\n", NNBr, 720);
	//for(j=1; j<NNBt;j+=NNBt/720){
	//	for(i=1; i<=NNBr; i++){
	//		fprintf(data1, "%e %e %e \n", Rvb[i][j]*Rbo, Tvb[i][j], Acv[i][j]*Rbo*Rbo);
	//	}
	///*fprintf(data1, "%e %e %e \n", Rvb[i][NNBt]*Rbo, Tvb[i][NNBt], Acv[i][NNBt]*Rbo*Rbo);*/
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
		if(phi>=2*Pi){
			//Qvo = 0.0;			
			if(hout>=phiout_c){			
				phit = phi;
				while(phit > 2*Pi) phit = phit - 2*Pi;
				/* Cylinder pressure data */
				//printf("\n%f\n", phi*180/Pi);
				//fprintf(fp,"%4.1f  ", phi*180/Pi);
				//for(i=1; i<=NRK; i++) fprintf(fp, "%e  ", Pc[i]);
				Qvo = 0.0, Qvi = 0.0;		// Outlet flow
				for(i=1; i<=N; i++){
					ph[i] = phit + 2*Pi/N*(i-1);
					while(ph[i] > 2*Pi) ph[i] = ph[i] - 2*Pi;
					if(ph[i]<=Pi){
						if(Pc[i]>Pc[N+1]) Qvo = Qvo + Cvo*Avo(ph[i])*sqrt(2/rho)*sqrt(Pc[i]-Pc[N+1]);
						else Qvo = Qvo - Cvo*Avo(ph[i])*sqrt(2/rho)*sqrt(Pc[N+1]-Pc[i]);
					}
					else{
						if(Pc[i]>Pc[N+2]) Qvi = Qvi + Cvi*Avi(ph[i])*sqrt(2/rho)*sqrt(Pc[i]-Pc[N+2]);
						else Qvi = Qvi - Cvi*Avi(ph[i])*sqrt(2/rho)*sqrt(Pc[N+2]-Pc[i]);
					}		
				}		
				Pc[N+3] = Qvo;
				Pc[N+4] = Qvi;
				// Outlet flow
				//fprintf(fp,"%e ", Qvo);
				//fprintf(fp,"%e ", Avo(phi));
				//Leakage(phi, Pc, &Qp, &Qv, &Qs, &Qpt, &Qvt, &Qst);	// Leakage in a Cylinder & Total leakage								// Resultant force action point
				//fprintf(fp,"%e %e %e %e %e %e %e ", Qp, Qv, Qs, Qpt, Qvt, Qst, Qpt+Qvt+Qst);
				//fprintf(fp,"%e %e ", eBVx, eBVy);					
				
				//printf("##### %e %e \n", phi, hout);
				for(i=1; i<=N+2; i++) Pcn[i] = Pc[i]/Lambda;		// Re-define nondimensional parameters in cylinder pressure calculation		
			/*	II.FORCE & MOMENT */						
				hcal = 0.0;
				hmcal = hm;
				while(hcal<hout){
					if(hout-hcal < hm*180/Pi) hmcal = (hout-hcal)*Pi/180;
			/*	III.BODY DYNAMICS */
					printf("\n --------- %f \n",phit*180/Pi);
					motion_RK4(phit, hmcal, ym, gamma, zeta, dgamma, dzeta, Hvb, Hn, Hs, Hw, He, Hvb_gro, 
								Rvb, Tvb, dTvb, dRvb, Rvbn, Rvbs, Rvbw, Rvbe, Tvbn, Tvbs, Tvbw, Tvbe, 
								dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Pvb, Pcn, Acv, 
								configuration, filmthickness, descretize, solver, FM);
					phit += hmcal;		// Motion angle step
					hcal += hmcal*180/Pi;
				}
				//printf("################################################################");
			/*	IV.PRESSURE DISTRIBUTION */	
				//configuration(phit, ym, &gamma, &zeta, &dgamma, &dzeta);
				//filmthickness(phit, ym, Hvb, Hn, Hs, Hw, He, Hvb_gro, Rvb, Tvb, gamma, zeta);
				//descretize(phit, ym, gamma, zeta, dgamma, dzeta, Rvb, dTvb, dRvb, Hn, Hs, Hw, He, Rvbn, Rvbs, Rvbw, Rvbe, 
				//			Tvbn, Tvbs, Tvbw, Tvbe, dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B);
				//solver(phit, Pcn, Pvb, aN, aS, aW, aE, aC, B);
			/*  V. Leakage & Friction */
				Pgradient(Pvb, dPdR_n, dPdT_c, dNvb, dSvb, dTvb);
				LeakageVB(leakage, Baleakage, Pvb, dNvb, dSvb, dTvb, dRvb, Rvbn, Hvb, Hn, dPdR_n, dPdT_c);
				FrictionVB(phit, &Tb, &Tv, Rvb, Hvb, Acv, dPdT_c, Pvb);
			/*  VI. SAVE Cylinder Pressure, Film Pressure, and Motion data */
				savedata(&n, data0, data1, data2, data3, data4, data5, data6, phi-2*Pi, &hout_cy, &hout_dy, &hout_vb, &hout_FM, &hout_LK, hout, Pc, Rvb, Tvb, Pvb, ym, Baleakage, Acv, Hvb);	
				hout = 0.0;
				Nd++;
			}
			hout += h*180/Pi;
		}
		Derivs(phi, Pc, dPc, Baleakage);
		for(i=1; i<=NRK; i++) Pscal[i] = fabs(Pc[i]) + fabs(h*dPc[i]) + TINY;
		RK4c(Pc, dPc, NRK, &phi, h, Pscal, &hdid, &hnext, Derivs, Baleakage);
		for(i=1; i<=N+2; i++) if(Pc[i] < Pca) Pc[i] = Pca;
		printf("%e \n", phi*180/Pi);
		h = hnext;
	}	

	free_dvector(ph,1,N);
	free_dvector(Pc,1,NRK);
	free_dvector(Pcn,1,NRK);
	free_dvector(dPc,1,NRK);
	free_dvector(Pscal,1,NRK);

	free_dvector(ym,1,2*DOF);
	free_matrix(Rvb,1,NNBr,1,NNBt);  // Rad. coord. matrix
	free_matrix(Rvbn,1,NNBr,1,NNBt);  // North Rad. coord. matrix
	free_matrix(Rvbs,1,NNBr,1,NNBt);  // South Rad. 
	free_matrix(Rvbw,1,NNBr,1,NNBt);  // West Rad. 
	free_matrix(Rvbe,1,NNBr,1,NNBt);  // East Rad. 
	free_matrix(Tvb,1,NNBr,1,NNBt);  // Circumf. coord. matrix
	free_matrix(Tvbn,1,NNBr,1,NNBt);  // North Circumf. coord. matrix
	free_matrix(Tvbs,1,NNBr,1,NNBt);  // South Circumf. 
	free_matrix(Tvbw,1,NNBr,1,NNBt);  // West Circumf. 
	free_matrix(Tvbe,1,NNBr,1,NNBt);  // East Circumf. 
	free_matrix(dNvb,1,NNBr,1,NNBt);  // North Rad. coord. matrix
	free_matrix(dSvb,1,NNBr,1,NNBt);  // South Rad. 
	free_matrix(dWvb,1,NNBr,1,NNBt);  // West Rad. 
	free_matrix(dEvb,1,NNBr,1,NNBt);  // East Rad. 
	free_matrix(dRvb,1,NNBr,1,NNBt);  // delta Rad. CV 
	free_matrix(dTvb,1,NNBr,1,NNBt);  // delta Circumf. CV 
	free_matrix(Acv,1,NNBr,1,NNBt);
	free_matrix(Pvb,1,NNBr,1,NNBt);  // V-B Pressure matrix
	free_dmatrix(Hvb,1,NNBr,1,NNBt);  // V-B Film thickness matrix
	free_dmatrix(Hn,1,NNBr,1,NNBt);
	free_dmatrix(Hs,1,NNBr,1,NNBt);
	free_dmatrix(Hw,1,NNBr,1,NNBt);
	free_dmatrix(He,1,NNBr,1,NNBt);
	free_dmatrix(Hvb_gro,1,NNBr,1,NNBt);  // V-B Groove matrix
	free_matrix(aN,1,NNBr,1,NNBt);  // 차분화 계수aN, aS, aW, aE, B
	free_matrix(aS,1,NNBr,1,NNBt);  
	free_matrix(aW,1,NNBr,1,NNBt);
	free_matrix(aE,1,NNBr,1,NNBt);
	free_matrix(aC,1,NNBr,1,NNBt);
	free_matrix(B,1,NNBr,1,NNBt);
	
	free_vector(leakage,1,NNBr);
	free_vector(Baleakage,1,N);

	free_matrix(dPdR_n,1,NNBr,1,NNBt);
	free_matrix(dPdT_c,1,NNBr,1,NNBt);

	fclose(data0);
	fclose(data1);
	fclose(data2);
	fclose(data3);
	fclose(data4);
	fclose(data5);
	fclose(data6);
 
  return 0;
}
/* Save output data */
void savedata(int *n, FILE *data0, FILE *data1, FILE *data2, FILE *data3, FILE *data4, FILE *data5, FILE *data6, double phi, double *hout_cy, double *hout_dy, double *hout_vb, double *hout_FM, double *hout_LK, double hout, double  *Pc, float **Rvb, float **Tvb, float **Pvb, double *ym, float *Baleakage, float **Acv, double **Hvb)
{
	int i, j, k;
	double phit;
	float eVBx, eVBy, Pe, eBVx, eBVy;

	void Point_resul(double , double *, float *, float *);
	
	*hout_cy += hout;
	*hout_dy += hout;
	*hout_vb += hout;

	phit = phi;
	while(phit > 2*Pi) phit = phit - 2*Pi;

	//printf("################################################################");
	// Cylinder Pressure
	if(*hout_cy >= phiout_cy){
		fprintf(data0,"%4.1f  ", phi*180/Pi);
		for(i=1; i<=NRK; i++) fprintf(data0, "%e  ", Pc[i]);
		fprintf(data0,"\n");
		*hout_cy = 0.0;
	}

	// Pressure_VB & Film Thickness 
	if(*hout_vb >= phiout_vb){
		fprintf(data1, "NCYCLE=%d\n VARIABLES=R,T,PRESSURE(Pa)\n ZONE T=P_Film, I=         %d ,J=        %d , F=POINT\n", *n, 60, 361);
		fprintf(data3, "NCYCLE=%d\n VARIABLES=R,T,FThickness(um)\n ZONE T=P_Film, I=         %d ,J=        %d , F=POINT\n", *n, 60, 361);
		//fprintf(data3, "NCYCLE=%d\n VARIABLES=eVBx, eVBy\n ZONE T=P_Film, I=         1 ,J=        1 , F=POINT\n", *n);
		for(j=1; j<=NNBt; j+=NNBt/360){
			//for(i=1; i<=NNBr; i++){ 
			//	//while(Tvb[i][j]+phit >= 2*Pi) Tvb[i][j] = Tvb[i][j] - 2*Pi;
			//	fprintf(data1, "%e %e %e \n", Rvb[i][j]*Rbo, Pi/2-(Tvb[i][j]+phit), Pvb[i][j]*Lambda);	
			//	//fprintf(data3, "%e %e %e \n", Rvb[i][j]*Rbo, Pi/2-(Tvb[i][j]+phit), (ym[1]*hv0-Rvb[i][j]*Rbo*acos(cos(ym[2])*cos(ym[3]))*cos(Tvb[i][j]-Pi/2-atan2(tan(ym[2]),sin(ym[3]))+phit))*1000000);
			//}
			for(i=1; i<1+NBr1; i++/*=NBr1/NBr1*/){ 
				fprintf(data1, "%e %e %e \n", Rvb[i][j]*Rbo, Pi/2-(Tvb[i][j]+phit), Pvb[i][j]*Lambda);
				fprintf(data3, "%e %e %e \n", Rvb[i][j]*Rbo, Pi/2-(Tvb[i][j]+phit), Hvb[i][j]*hv0);				
			}
			for(i=1+NBr1; i<1+NBr1+NBr2; i+=NBr2/9){ 
				fprintf(data1, "%e %e %e \n", Rvb[i][j]*Rbo, Pi/2-(Tvb[i][j]+phit), Pvb[i][j]*Lambda);				
				fprintf(data3, "%e %e %e \n", Rvb[i][j]*Rbo, Pi/2-(Tvb[i][j]+phit), Hvb[i][j]*hv0);	
			}
			for(i=1+NBr1+NBr2; i<1+NBr1+NBr2+NBr3; i++/*=NBr3/NBr3*/){ 
				fprintf(data1, "%e %e %e \n", Rvb[i][j]*Rbo, Pi/2-(Tvb[i][j]+phit), Pvb[i][j]*Lambda);
				fprintf(data3, "%e %e %e \n", Rvb[i][j]*Rbo, Pi/2-(Tvb[i][j]+phit), Hvb[i][j]*hv0);	
			}
			for(i=1+NBr1+NBr2+NBr3; i<1+NBr1+NBr2+NBr3+NBr4; i++/*=NBr4/4*/){ 
				fprintf(data1, "%e %e %e \n", Rvb[i][j]*Rbo, Pi/2-(Tvb[i][j]+phit), Pvb[i][j]*Lambda);
				fprintf(data3, "%e %e %e \n", Rvb[i][j]*Rbo, Pi/2-(Tvb[i][j]+phit), Hvb[i][j]*hv0);	
			}
			for(i=1+NBr1+NBr2+NBr3+NBr4; i<=NNBr; i++/*=NBr5/NBr5*/){ 
				fprintf(data1, "%e %e %e \n", Rvb[i][j]*Rbo, Pi/2-(Tvb[i][j]+phit), Pvb[i][j]*Lambda);
				fprintf(data3, "%e %e %e \n", Rvb[i][j]*Rbo, Pi/2-(Tvb[i][j]+phit), Hvb[i][j]*hv0);	
			}
		}
		for(i=1; i<1+NBr1; i++){ 
			fprintf(data1, "%e %e %e \n", Rvb[i][1]*Rbo, Pi/2-(Tvb[i][1]+phit), Pvb[i][1]*Lambda);
			fprintf(data3, "%e %e %e \n", Rvb[i][1]*Rbo, Pi/2-(Tvb[i][1]+phit), Hvb[i][1]*hv0);	
		}
		for(i=1+NBr1; i<1+NBr1+NBr2; i+=NBr2/9){ 
			fprintf(data1, "%e %e %e \n", Rvb[i][1]*Rbo, Pi/2-(Tvb[i][1]+phit), Pvb[i][1]*Lambda);
			fprintf(data3, "%e %e %e \n", Rvb[i][1]*Rbo, Pi/2-(Tvb[i][1]+phit), Hvb[i][1]*hv0);	
		}
		for(i=1+NBr1+NBr2; i<1+NBr1+NBr2+NBr3; i++/*=NBr3/NBr3*/){ 
			fprintf(data1, "%e %e %e \n", Rvb[i][1]*Rbo, Pi/2-(Tvb[i][1]+phit), Pvb[i][1]*Lambda);
			fprintf(data3, "%e %e %e \n", Rvb[i][1]*Rbo, Pi/2-(Tvb[i][1]+phit), Hvb[i][1]*hv0);	
		}
		for(i=1+NBr1+NBr2+NBr3; i<1+NBr1+NBr2+NBr3+NBr4; i++/*=NBr4/4*/){ 
			fprintf(data1, "%e %e %e \n", Rvb[i][1]*Rbo, Pi/2-(Tvb[i][1]+phit), Pvb[i][1]*Lambda);
			fprintf(data3, "%e %e %e \n", Rvb[i][1]*Rbo, Pi/2-(Tvb[i][1]+phit), Hvb[i][1]*hv0);	
		}
		for(i=1+NBr1+NBr2+NBr3+NBr4; i<=NNBr; i++/*=NBr5/NBr5*/){ 
			fprintf(data1, "%e %e %e \n", Rvb[i][1]*Rbo, Pi/2-(Tvb[i][1]+phit), Pvb[i][1]*Lambda);
			fprintf(data3, "%e %e %e \n", Rvb[i][1]*Rbo, Pi/2-(Tvb[i][1]+phit), Hvb[i][1]*hv0);	
		}
		fprintf(data1,"\n");		
		fprintf(data3,"\n");

		*hout_vb = 0.0;
		*n = *n + 1;
	}		
	
	// Barrel Motion
	if(*hout_dy >= phiout_dy){				
		//printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
		int Nphit1 = NNBt/4-int(phit/dNBt);
		int Nphit2 = NNBt*7/12-int(phit/dNBt);
		int Nphit3 = NNBt*11/12-int(phit/dNBt);
		while(Nphit1 <= 0) Nphit1 = Nphit1 + NNBt;
		while(Nphit2 <= 0) Nphit2 = Nphit2 + NNBt;
		while(Nphit3 <= 0) Nphit3 = Nphit3 + NNBt;
		while((Hvb[NNBr][Nphit1]*hv0)*1e6 > 500) Hvb[NNBr][Nphit1] = Hvb[NNBr][Nphit1] - 1000e-6/hv0;
		while((Hvb[NNBr][Nphit2]*hv0)*1e6 > 500) Hvb[NNBr][Nphit2] = Hvb[NNBr][Nphit2] - 1000e-6/hv0;
		while((Hvb[NNBr][Nphit3]*hv0)*1e6 > 500) Hvb[NNBr][Nphit3] = Hvb[NNBr][Nphit3] - 1000e-6/hv0;
		fprintf(data2, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",	phi*180/Pi, ym[1]*hv0,  acos(cos(ym[2])*cos(ym[3]))*180/Pi, -atan2(tan(ym[2]),sin(ym[3]))*180/Pi, ym[2]*180/Pi, ym[3]*180/Pi, ym[4]*hv0*w, ym[5]*w*180/Pi, ym[6]*w*180/Pi, Hmin*hv0/*ym[4]*hv0-Rbo*ym[5]*hv0/Rbo*/,
			(Hvb[NNBr][Nphit1]*hv0)*1000000, (Hvb[NNBr][Nphit2]*hv0)*1000000, (Hvb[NNBr][Nphit3]*hv0)*1000000, (ym[1]*hv0-Rvb[NNBr][NNBt/4]*Rbo*acos(cos(ym[2])*cos(ym[3]))*cos(Pi/2-phit-Pi/2-atan2(tan(ym[2]),sin(ym[3]))+phit))*1000000, (ym[1]*hv0-Rvb[NNBr][NNBt/4]*Rbo*acos(cos(ym[2])*cos(ym[3]))*cos(Pi*7/6-phit-Pi/2-atan2(tan(ym[2]),sin(ym[3]))+phit))*1000000, (ym[1]*hv0-Rvb[NNBr][NNBt/4]*Rbo*acos(cos(ym[2])*cos(ym[3]))*cos(Pi*11/6-phit-Pi/2-atan2(tan(ym[2]),sin(ym[3]))+phit))*1000000, Rvb[Imin][Jmin]*Rbo, Pi/2-(Tvb[Imin][Jmin]+phit) );
			
		*hout_dy = 0.0;

		eVBx = 0.0, eVBy = 0.0, Pe =0.0;
		Point_resul(phi, Pc, &eBVx, &eBVy);	
		for(j=1; j<=NNBt; j+=NNBt/540){
			for(i=1; i<=NNBr; i++){ 
				eVBx = eVBx + Pvb[i][j]*Lambda*Rvb[i][j]*Rbo*cos(Pi/2-(Tvb[i][j]+phit))*Acv[i][j];
				eVBy = eVBy + Pvb[i][j]*Lambda*Rvb[i][j]*Rbo*sin(Pi/2-(Tvb[i][j]+phit))*Acv[i][j];
				Pe = Pe + Pvb[i][j]*Lambda*Acv[i][j];
			}
		}	
		eVBx = eVBx/Pe;
		eVBy = eVBy/Pe;
		fprintf(data6, "%e %e %e %e \n", eVBx, eVBy, eBVx, eBVy);
	}

	// Force & Moment on the Barrel
	*hout_FM += hout;
	phiout_FM = 0.1;
	// Barrel Motion
	if(*hout_FM >= phiout_FM){				
		//printf("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
		fprintf(data4, "%4.1f %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n",	phi*180/Pi,	-FliftD/FpressD, -FpressD, MpressxD, MpressyD, FliftD, MliftxD, MliftyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, FzD, MxD, MyD, MzD);
		*hout_FM = 0.0;
	}

	// Leakage & Friction loss
	*hout_LK += hout;
	phiout_LK = 1;
	if(*hout_LK >= phiout_LK){	
		fprintf(data5, "%4.1f  ", phi*180/Pi);	
		for(i=1; i<=N; i++) fprintf(data5, "%e  ", Baleakage[i]);
		fprintf(data5, "%e  %e  %e  %e  %e  %e  %e  %e %e \n", -leakageRinD, leakageRoutD, -leakageRinD+leakageRoutD, frictorqueD, (Pout-Pin)*(-leakageRinD+leakageRoutD), friclossD, friclossD + (Pout-Pin)*(-leakageRinD+leakageRoutD), frictorqueD2*w, frictorqueD2*w + (Pout-Pin)*(-leakageRinD+leakageRoutD));
		*hout_LK = 0.0;
	}

	// Test - Boundary condition
	//fprintf(data1, "NCYCLE=%d\n VARIABLES=R,T,PRESSURE(Pa)\n ZONE T=P_Film, I=         %d ,J=        %d , F=POINT\n", *n, NNBr, 721);
	//for(j=1; j<NNBt; j+=NNBt/720){
	//	for(i=1; i<=NNBr; i++){ 				
	//		fprintf(data1, "%e %e %e \n", Rvb[i][j]*Rbo, Tvb[i][j], Acv[i][j]*Rbo*Rbo);
	//	}
	//}
	//for(i=1; i<=NNBr; i++){ 
	//	fprintf(data1, "%e %e %e \n", Rvb[i][NNBt]*Rbo, Tvb[i][NNBt], Acv[i][NNBt]*Rbo*Rbo);
	//}
	//fprintf(data1,"\n");

}