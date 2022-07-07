#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

extern float Vp, **dPdZ_n, **dPdT_c;
extern void Pgradient(float **Ppc, float **dPdZ_n, float **dPdT_c, float **dN, float **dS, float **dT);
extern void FrictionPC(double phi, float *Ffric, float *Mfricx, float *Mfricy, float **Zpc, float **Tpc, double **Hpc, float **Acv, float **dPdZ_n, float **Ppc);
extern int Nphi, NBz1, NBz2, NBz3, NBz4, NBz5, NBz6, NBz, NNBz, NNBt, Nch, Imin, Jmin;
extern int NVt1, NVt2, NVt3, NVt4, NVt5, NVt6, NVr1, NVr2, NVr3, NVr4;
extern double dNBt, thetach;
extern double R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, revf, phif, phiout_c, hm, phiout_cy, phiout_dy, phiout_vb, phiout_FM, phiout_LK, w, wp;
extern double V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
extern double integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, zB1, zB2, zB3, zB4, zB5, zB6, zB7, Fs, e0, theta0, psi0, de0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
extern double Lambda, Loc, Lpj, Lpg, Lbj, Ls, Rbo, Lpc;
extern float FxhD, FyhD, MxhD, MyhD, FxcD, FycD, MxcD, MycD, FxD, FyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, FzD, MxD, MyD, MzD, leakageRinD, leakageRoutD, frictionD, friclossD, Hmin;
extern double gamH, zetaH;
/* Solve equations of motion */
/* Runge-Kutta 4th */
void motion_RK4(double *phim, int mid, double hrem, double htry, double *hn, double *ym, double e, double kai, double gamma, double zeta, double de, double dkai,double dgamma, double dzeta, double **Hpc, double **Hn, double **Hs, double **Hw, double **He, double **Hpc_gro, 
				float **Zpc, float **Tpc, float **dTpc, float **dZpc, float **Zpcn, float **Zpcs, float **Zpcw, float **Zpce, float **Tpcn, float **Tpcs, float **Tpcw, float **Tpce, 
				float **dNvb, float **dSvb, float **dWvb, float **dEvb, float **aN, float **aS, float **aW, float **aE, float **aC, float **B, float **Ppc, double *Pcn, float **Acv, 
				void (*configuration)(double , double *, double *, double *, double *, double *, double *, double *, double *, double *), 
				void (*filmthickness)(double , double *, double **, double **, double **, double **, double **, double **, float **, float **, double , double, double , double ), 
				void (*descretize)(double , double *, double , double , double , double , double , double , double , double , float **, float **, float **, double **, double **, double **, double **, float **, float **, 
									float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **),
				void (*solver)(double, double *, float **, float **, float **, float **, float **, float **, float **), 
				void (*FM)(double , double *, float **, float **, float **, float **, float **, float **, float **, float **, float **, float *, float *, float *, float *, float *, float *, double *, double **, float **) ) 
{
	void Derivsm(double , double *, double *, double e, double kai, double gamma, double zeta, double de, double dkai,double dgamma, double dzeta, double **Hpc, double **Hn, double **Hs, double **Hw, double **He, double **Hpc_gro, 
				float **Zpc, float **Tpc, float **dTpc, float **dZpc, float **Zpcn, float **Zpcs, float **Zpcw, float **Zpce, float **Tpcn, float **Tpcs, float **Tpcw, float **Tpce, 
				float **dNvb, float **dSvb, float **dWvb, float **dEvb, float **aN, float **aS, float **aW, float **aE, float **aC, float **B, float **Ppc, double *Pcn, float **Acv,
				void (*configuration)(double , double *, double *, double *, double *, double *, double *, double *, double *, double *), 
				void (*filmthickness)(double , double *, double **, double **, double **, double **, double **, double **, float **, float **, double , double , double , double ), 
				void (*descretize)(double , double *, double , double , double , double , double , double , double , double , float **, float **, float **, double **, double **, double **, double **, float **, float **, 
									float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **),
				void (*solver)(double, double *, float **, float **, float **, float **, float **, float **, float **), 
				void (*FM)(double , double *, float **, float **, float **, float **, float **, float **, float **, float **, float **, float *, float *, float *, float *, float *, float *, double *, double **, float **) ); 
	//void Derivsm_c(double , double *, double *, double gamma, double zeta, double dgamma, double dzeta, float **Hpc, float **Hn, float **Hs, float **Hw, float **He, float **Hpc_gro, 
	//			float **Zpc, float **Tpc, float **dTpc, float **dZpc, float **Zpcn, float **Zpcs, float **Zpcw, float **Zpce, float **Tpcn, float **Tpcs, float **Tpcw, float **Tpce, 
	//			float **dNvb, float **dSvb, float **dWvb, float **dEvb, float **aN, float **aS, float **aW, float **aE, float **aC, float **B, float **Ppc, double *Pcn, float **Acv,
	//			void (*configuration)(double , double *, double *, double *, double *, double *), 
	//			void (*filmthickness)(double , double *, float **, float **, float **, float **, float **, float **, float **, float **, double , double ), 
	//			void (*descretize)(double , double *, double , double , double , double , float **, float **, float **, float **, float **, float **, float **, float **, float **, 
	//								float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **), 
	//			void (*solver)(double, double *, float **, float **, float **, float **, float **, float **, float **), 
	//			void (*FM)(double , double *, float **, float **, float **, float **, float **, float **, float **, float *, float *, float *, float *, float *, float *, double *));

	int i, j, jtry;	
	double *K1, *K2, *K3, *K4, phisav, *ysav, *yscal, phi;	
	double *err, errmax, h;
	//double *K1, *K2, *K3, *K4;
	double *slope, *ymR;

	K1 = dvector(1,2*DOF);
	K2 = dvector(1,2*DOF);
	//K3 = dvector(1,2*DOF);
	//K4 = dvector(1,2*DOF);
	slope = dvector(1,2*DOF);
	ymR = dvector(1,2*DOF);
	err = dvector(1,2*DOF);
	ysav = dvector(1,2*DOF);
	yscal = dvector(1,2*DOF);
	phisav = (*phim);

	if(mid == 0){
		Derivsm(phisav, ym, K1, e, kai, gamma, zeta, de, dkai, dgamma, dzeta, Hpc, Hn, Hs, Hw, He, Hpc_gro, 
					Zpc, Tpc, dTpc, dZpc, Zpcn, Zpcs, Zpcw, Zpce, Tpcn, Tpcs, Tpcw, Tpce, 
					dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Ppc, Pcn, Acv, 
					configuration, filmthickness, descretize, solver, FM);
		for(i=1; i<=2*DOF; i++){
			yscal[i] = fabs(ym[i]) + fabs(h*K1[i]) + TINY;
			ysav[i] = ym[i];
		}
		printf("===== hm = %e \n", htry*180/Pi);

		h = htry;
		for(jtry=1; jtry<=MAXTRY; jtry++){		
			//Derivsm(phisav, ym, K1, e, kai, gamma, zeta, de, dkai, dgamma, dzeta, Hpc, Hn, Hs, Hw, He, Hpc_gro, 
			//		Zpc, Tpc, dTpc, dZpc, Zpcn, Zpcs, Zpcw, Zpce, Tpcn, Tpcs, Tpcw, Tpce, 
			//		dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Ppc, Pcn, Acv, 
			//		configuration, filmthickness, descretize, solver, FM);
			for(i=1; i<=2*DOF; i++) ymR[i] = ysav[i] + K1[i]*h/2.;
			Derivsm(phisav+h, ymR, K2, e, kai, gamma, zeta, de, dkai, dgamma, dzeta, Hpc, Hn, Hs, Hw, He, Hpc_gro, 
					Zpc, Tpc, dTpc, dZpc, Zpcn, Zpcs, Zpcw, Zpce, Tpcn, Tpcs, Tpcw, Tpce, 
					dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Ppc, Pcn, Acv, 
					configuration, filmthickness, descretize, solver, FM);
			for(i=1; i<=2*DOF; i++){
				ym[i] = ysav[i] + (K1[i]+K2[i])*0.5*h;
				err[i] = (-K1[i]+K2[i])*0.5*h;
				//ymE[i] = ym[i] + K1[i]*h;
				//ymME[i] = ym[i] + (K1[i]+K2[i])*0.5*h;
				
			}	
			*phim = phisav + h;
			if(*phim == phisav) nrerror("stepsize not significant in Euler");
			errmax = 0.0;			//Evaluate accuracy		
			for(i=1; i<=2*DOF; i++) errmax = FMAX(errmax, fabs(err[i]/yscal[i]));
			//for(i=1; i<=n; i++) printf("%e \n", h);
			errmax /= epsm;			//scale Zelative to the ZequiZed tolerance 		
			if(errmax <= 1.0){
				//*hdid = h;
				*hn = (errmax > ERRCON ? SAFETY*h*pow(errmax,PGROW) : GROW*h);
				if(*hn > hrem) *hn = hrem;
				//printf("##################################### %e \n",*hn);
				free_dvector(ysav,1,2*DOF);
				free_dvector(yscal,1,2*DOF);
				free_dvector(err,1,2*DOF);
				free_dvector(K1,1,2*DOF);
				free_dvector(K2,1,2*DOF);
				free_dvector(slope,1,2*DOF);
				free_dvector(ymR,1,2*DOF);
				return;
			} 
			else{				//Truncation error too large, Zeduce stepsize
				*hn = SAFETY*h*pow(errmax,PSHRNK);
				h = (h >= 0.0 ? FMAX(*hn,SHRNK*h) : FMIN(*hn,SHRNK*h));
			}
		}
		nrerror("exceeded MAXTRY in RKstep");
	}
	else{
		Derivsm(phisav, ym, K1, e, kai, gamma, zeta, de, dkai, dgamma, dzeta, Hpc, Hn, Hs, Hw, He, Hpc_gro, 
				Zpc, Tpc, dTpc, dZpc, Zpcn, Zpcs, Zpcw, Zpce, Tpcn, Tpcs, Tpcw, Tpce, 
				dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Ppc, Pcn, Acv, 
				configuration, filmthickness, descretize, solver, FM);
		//for(i=1; i<=2*DOF; i++) ymR[i] = ym[i] + K1[i]*h/2.;
		//Derivsm(phi+h, ymR, K2, e, kai, gamma, zeta, de, dkai, dgamma, dzeta, Hpc, Hn, Hs, Hw, He, Hpc_gro, 
		//		Zpc, Tpc, dTpc, dZpc, Zpcn, Zpcs, Zpcw, Zpce, Tpcn, Tpcs, Tpcw, Tpce, 
		//		dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Ppc, Pcn, Acv, 
		//		configuration, filmthickness, descretize, solver, FM);
		//for(i=1; i<=2*DOF; i++) ymR[i] = ym[i] + K2[i]*h/2.;
		//Derivsm(phi+h/2, ymR, K3, gamma, zeta, dgamma, dzeta, Hpc, Hn, Hs, Hw, He, Hpc_gro, 
		//		Zpc, Tpc, dTpc, dZpc, Zpcn, Zpcs, Zpcw, Zpce, Tpcn, Tpcs, Tpcw, Tpce, 
		//		dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Ppc, Pcn, Acv, 
		//		configuration, filmthickness, descretize, solver, FM);
		//for(i=1; i<=2*DOF; i++) ymR[i] = ym[i] + K3[i]*h;	
		//Derivsm(phi+h, ymR, K4, gamma, zeta, dgamma, dzeta, Hpc, Hn, Hs, Hw, He, Hpc_gro, 
		//		Zpc, Tpc, dTpc, dZpc, Zpcn, Zpcs, Zpcw, Zpce, Tpcn, Tpcs, Tpcw, Tpce, 
		//		dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Ppc, Pcn, Acv, 
		//		configuration, filmthickness, descretize, solver, FM);
		for(i=1; i<=2*DOF; i++){
			//slope[i] = (K1[i]+2.0*(K2[i]+K3[i])+K4[i])/6.0;
			//ym[i] = ym[i] + slope[i]*h;
			ym[i] = ym[i] + K1[i]*h;
		}	

		free_dvector(ysav,1,2*DOF);
		free_dvector(yscal,1,2*DOF);
		free_dvector(err,1,2*DOF);
		free_dvector(K1,1,2*DOF);
		free_dvector(K2,1,2*DOF);
		free_dvector(slope,1,2*DOF);
		free_dvector(ymR,1,2*DOF);
	}

}
/* phi, y, dy */
void Derivsm( double phi, double *ym, double *dy, double e, double kai, double gamma, double zeta, double de, double dkai,double dgamma, double dzeta, double **Hpc, double **Hn, double **Hs, double **Hw, double **He, double **Hpc_gro, 
				float **Zpc, float **Tpc, float **dTpc, float **dZpc, float **Zpcn, float **Zpcs, float **Zpcw, float **Zpce, float **Tpcn, float **Tpcs, float **Tpcw, float **Tpce, 
				float **dNvb, float **dSvb, float **dWvb, float **dEvb, float **aN, float **aS, float **aW, float **aE, float **aC, float **B, float **Ppc, double *Pcn, float **Acv,
				void (*configuration)(double phi, double *ym, double *e, double *kai, double *gamma, double *zeta, double *de, double *dkai, double *dgamma, double *dzeta), 
				void (*filmthickness)(double , double *, double **, double **, double **, double **, double **, double **, float **, float **, double , double, double , double ), 
				void (*descretize)(double , double *, double , double , double , double , double , double , double , double , float **, float **, float **, double **, double **, double **, double **, float **, float **, 
									float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **),
				void (*solver)(double, double *, float **, float **, float **, float **, float **, float **, float **), 
				void (*FM)(double , double *, float **, float **, float **, float **, float **, float **, float **, float **, float **, float *, float *, float *, float *, float *, float *, double *, double **, float **) ) 
{
	float Fx, Fy, Fz, Mx, My, Mz;
	float Ru, Lu, Mu;
	double Ibxx, Ibyy, Ibzz, Ibxy, Ibzx, Ibyz, dIbxx, dIbyy, dIbzz, dIbxy, dIbzx, dIbyz;

	(*configuration)(phi, ym, &e, &kai, &gamma, &zeta, &de, &dkai, &dgamma, &dzeta);
	(*filmthickness)(phi, ym, Hpc, Hn, Hs, Hw, He, Hpc_gro, Zpc, Tpc, e, kai, gamma, zeta);
	(*descretize)(phi, ym, e, kai, gamma, zeta, de, dkai, dgamma, dzeta, Zpc, dTpc, dZpc, Hn, Hs, Hw, He, Zpcn, Zpcs, Zpcw, Zpce, 
				Tpcn, Tpcs, Tpcw, Tpce, dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B);
	(*solver)(phi, Pcn, Ppc, aN, aS, aW, aE, aC, B);
	(*FM)(phi, Pcn, Zpc, Ppc, Zpcn, Zpcs, Tpc, dNvb, dSvb, dTpc, Acv, &Fx, &Fy, &Fz, &Mx, &My, &Mz, ym, Hpc, dPdZ_n);

	Ru = Rbo, Lu = 40.0e-3, Mu = 15e-3; 
	Mu = Mu*w*w/Lambda/Rbo;
	Ru = Ru/Rbo;
	Lu = Lu/Rbo;
	Ibxx = ( Ibt + Mu*(pow(Ru*cos(phi),2)+Lu*Lu) );
	Ibyy = ( Ibt + Mu*(pow(Ru*sin(phi),2)+Lu*Lu) );
	Ibzz = ( Iba + Mu*Ru*Ru );
	Ibxy = Mu*Ru*Ru*sin(phi)*cos(phi);
	Ibzx = Mu*Ru*Lu*sin(phi);
	Ibyz = Mu*Ru*Lu*cos(phi);
	dIbxx = Mu*pow(Ru,2)*sin(2*phi);
	dIbyy = -Mu*pow(Ru,2)*sin(2*phi);
	dIbzz = 0.0;
	dIbxy = Mu*Ru*Ru*cos(2*phi);
	dIbzx = -Mu*Ru*Lu*cos(phi);
	dIbyz = Mu*Ru*Lu*sin(phi);

 
	dy[1] = ym[5];
	dy[2] = ym[6];
	dy[3] = ym[7];
	dy[4] = ym[8];
	dy[5] = Fx/Mb;
	dy[6] = Fy/Mb;
	dy[7] = Mx/Ibt;
	dy[8] = My/Ibt;
	//dy[1] = 0;
	//dy[2] = 0;
	//dy[3] = 0;
	//dy[4] = 0;
	//dy[5] = 0;
	//dy[6] = 0;
	//dy[7] = 0;
	//dy[8] = 0;

	//printf("%3.3e %3.3e %3.3e %3.3e %3.3e %3.3e \n", dy[1], dy[2], dy[3], dy[4], dy[5], dy[6]);
}
/* phi, y, dy incase of boundary lubrication */
//void Derivsm_c( double phi, double *ym, double *dy, double gamma, double zeta, double dgamma, double dzeta, float **Hpc, float **Hn, float **Hs, float **Hw, float **He, float **Hpc_gro, 
//				float **Zpc, float **Tpc, float **dTpc, float **dZpc, float **Zpcn, float **Zpcs, float **Zpcw, float **Zpce, float **Tpcn, float **Tpcs, float **Tpcw, float **Tpce, 
//				float **dNvb, float **dSvb, float **dWvb, float **dEvb, float **aN, float **aS, float **aW, float **aE, float **aC, float **B, float **Ppc, double *Pcn, float **Acv,
//				void (*configuration)(double , double *, double *, double *, double *, double *), 
//				void (*filmthickness)(double , double *, float **, float **, float **, float **, float **, float **, float **, float **, double , double ), 
//				void (*descretize)(double , double *, double , double , double , double , float **, float **, float **, float **, float **, float **, float **, float **, float **, 
//									float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **), 
//				void (*solver)(double, double *, float **, float **, float **, float **, float **, float **, float **), 
//				void (*FM)(double , double *, float **, float **, float **, float **, float **, float **, float **, float *, float *, float *, float *, float *, float *, double *)) 
//{
//	float Fx, Fy, fr, Mx, My, Mz;
//	float Ax, Bx, Kx, Ay, By, Ky;
//
//	(*configuration)(phi, ym, &gamma, &zeta, &dgamma, &dzeta);
//	(*filmthickness)(phi, ym, Hpc, Hn, Hs, Hw, He, Hpc_gro, Zpc, Tpc, gamma, zeta);
//	(*descretize)(phi, ym, gamma, zeta, dgamma, dzeta, Zpc, dTpc, dZpc, Hn, Hs, Hw, He, Zpcn, Zpcs, Zpcw, Zpce, 
//				Tpcn, Tpcs, Tpcw, Tpce, dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B);
//	(*solver)(phi, Pcn, Ppc, aN, aS, aW, aE, aC, B);
//	(*FM)(phi, Pcn, Zpc, Ppc, Zpcn, Zpcs, Tpc, dTpc, Acv, &Fx, &Fy, &fr, &Mx, &My, &Mz, ym);
//
//	Ax = Ibt - sin(zeta)*Mb*sin(ym[2])*cos(ym[3])/gamma*pow(Rp/hp,2)*w;
//	Bx = -sin(zeta)*Mb*cos(ym[2])*sin(ym[3])/gamma*pow(Rp/hp,2)*w;
//	Kx = Mx + sin(zeta)*(Mb*(-dgamma*dgamma*pow(hp/Rp,2)*cos(ym[2])*cos(ym[3])+cos(ym[2])*cos(ym[3])*(ym[5]*ym[5]+ym[6]*ym[6])-2*sin(ym[2])*sin(ym[3]))/gamma*pow(Rp/hp,2)*w-fr)
//		-Ibt*pow(ym[6]*cos(ym[2]),2)*ym[2] - Iba*ym[6]*cos(ym[2])*cos(ym[2])*cos(ym[3]);
//	Ay = Ibt*cos(ym[2]) + cos(zeta)*Mb*cos(ym[2])*sin(ym[3])/gamma*pow(Rp/hp,2)*w;
//	By = cos(zeta)*Mb*sin(ym[2])*cos(ym[3])/gamma*pow(Rp/hp,2)*w;
//	Ky = My - cos(zeta)*(Mb*(-dgamma*dgamma*pow(hp/Rp,2)*cos(ym[2])*cos(ym[3])+cos(ym[2])*cos(ym[3])*(ym[5]*ym[5]+ym[6]*ym[6])-2*sin(ym[2])*sin(ym[3]))/gamma*pow(Rp/hp,2)*w-fr)
//		+Ibt*ym[6]*(sin(ym[2])+cos(ym[2]))*ym[5] + Iba*ym[5]*cos(ym[2])*cos(ym[3]);		
//	dy[1] = 0.0/*ym[4]*/;
//	dy[2] = ym[5];
//	dy[3] = ym[6];
//	dy[5] = (Kx*By/Bx-Ky)/(Ax*By/Bx-Ay);
//	dy[6] = (Kx*Ay/Ax-Ky)/(Bx*Ay/Ax-By);
//	dy[4] = 0.0/*(-dgamma*dgamma*pow(hp/Rp,2)*cos(ym[2])*cos(ym[3])+cos(ym[2])*cos(ym[3])*(ym[5]*ym[5]+ym[6]*ym[6])-2*sin(ym[2])*sin(ym[3])+sin(ym[2])*cos(ym[3])*dy[2]+cos(ym[2])*sin(ym[3])*dy[3])/gamma*pow(Rp/hp,2)*w*/;
//
//	printf("%3.3e %3.3e %3.3e %3.3e %3.3e %3.3e \n", dy[1], dy[2], dy[3], dy[4], dy[5], dy[6]);
//}
void configuration(double phi, double *ym, double *e, double *kai, double *gamma, double *zeta, double *de, double *dkai, double *dgamma, double *dzeta)
{	
	// nondimension
	*e = sqrt(ym[1]*ym[1]+ym[2]*ym[2]);
	*kai = atan2(ym[2],ym[1]);
	*gamma = sqrt(ym[3]*ym[3]+ym[4]*ym[4])/*acos(cos(ym[3])*cos(ym[4]))*Rp/hp*/;
	*zeta = -atan2(ym[3],ym[4])/*-atan2(tan(ym[3]),sin(ym[4]))*/;
	*de = (ym[1]*ym[5]+ym[2]*ym[6])/(*e);
	*dkai = (ym[6]/ym[1]-ym[5]*ym[2]/ym[1]/ym[1])*cos(*kai)*cos(*kai);
	*dgamma = (ym[3]*ym[7]+ym[4]*ym[8])/(*gamma);
	*dzeta = (ym[7]/ym[4]-ym[8]*ym[3]/ym[4]/ym[4])*cos(*zeta)*cos(*zeta);

	//if(*gamma == 0. || *zeta == 0.){ 
	//	printf("Gamma or Zeta = 0.0");
	//	*gamma = TINY;
	//	*dgamma = (sin(ym[2])*cos(ym[3])*ym[5]+cos(ym[2])*sin(ym[3])*ym[6])/(*gamma)*pow(Rp/hp,2);
	//	//exit(1);
	//}
	//if(*zeta < 0) *zeta = *zeta + 2*Pi;
	printf("\nphi = %e \n %4.2e %4.2e %4.2e %4.2e %6.5e %4.3e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %e %d %d \n\n", phi*180/Pi, (*e)*hp, (*kai)*180/Pi, (*gamma)*180/Pi*hp/Rp, (*zeta)*180/Pi, ym[1]*hp, ym[2]*hp, ym[3]*180/Pi*hp/Rp, ym[4]*180/Pi*hp/Rp, ym[5]*hp*w, ym[6]*hp*w, ym[7]*w*180/Pi*hp/Rp, ym[8]*w*180/Pi*hp/Rp, (/*ym[1]-(*gamma)*/Hmin)*hp, Imin, Jmin);
	
}
/* Force & Moment of the Cylinder Barrel */
void FM(double phi, double *Pc, float **Zpc, float **Ppc, float **Zpcn, float **Zpcs, float **Tpc, float **dNvb, float **dSvb, float **dTpc, float **Acv, float *Fx, float *Fy, float *Fz, float *Mx, float *My, float *Mz, double *ym, double **Hpc, float **dPdZ_n)
{
	void FMc(double phi, float **Zpc, double *Pc, float *Fxc, float *Fyc, float *Mxc, float *Myc, double *ym);
	void FMh(double phi, float **R, float **Zn, float **Zs, float **Tpc, float **dT, float **Acv, float **P, float *Fxh, float *Fyh, float *Mxh, float *Myh);
	//void FrictionPC(double phi, float *Ffp, float *Mfxp, float *Mfyp, float **Zpc, double **Hpc, float **Acv, float **dPdZ_n, float **Ppc);
	//void FMfriction(double phi, double *Pc, float *Ffric, float *Mfricx, float *Mfricy);
	void FMlateral(double phi, double *Pc, float *Flatx, float *Flaty, float *Mlatx, float *Mlaty, float *Mlatz);

	float Fxc, Fyc, Mxc, Myc, Fxh, Fyh, Mxh, Myh, Ffric, Mfricx, Mfricy;
	float Flatx=0.0, Flaty=0.0, Mlatx=0.0, Mlaty=0.0, Mlatz, Fpa, Fc=0.0;

	Fpa = (Pc[1]*Lambda*Pi*Rp*Rp - Mp/w/w*Lambda*Rbo*(Rp/hp)*Ap + ((Pc[1]*Lambda-Ph)*hp/Lpc/2.0-eta*Vp/hp)*2*Pi*(Rp+hp)*Lpc)/Lambda/pow(Rbo,2);
		
	FMc(phi, Zpc, Pc, &Fxc, &Fyc, &Mxc, &Myc, ym);
	FMh(phi, Zpc, Zpcn, Zpcs, Tpc, dTpc, Acv, Ppc, &Fxh, &Fyh, &Mxh, &Myh);
	Pgradient(Ppc, dPdZ_n, dPdT_c, dNvb, dSvb, dTpc);
	FrictionPC(phi, &Ffric, &Mfricx, &Mfricy, Zpc, Tpc, Hpc, Acv, dPdZ_n, Ppc);
	//FMfriction(phi, Pc, &Ffric, &Mfricx, &Mfricy);
	FMlateral(phi, Pc, &Flatx, &Flaty, &Mlatx, &Mlaty, &Mlatz);

	if(Hmin < Hpcmin/hp) Fc = 9.0e11*(Hpcmin-Hmin*hp)/(Lambda*pow(Rp,2));

	*Fx = Fxc + Fxh + Flatx - Fc*sin(Tpc[Imin][Jmin]);
	*Fy = Fyc + Fyh + Flaty - Fc*cos(Tpc[Imin][Jmin]);
	*Fz = Pc[1]*Ap/Rp/Rp + Ffric;
	*Mx = Mxc + Mxh /*+ Mlatx */ + Mfricx + Fc*Zpc[Imin][Jmin]*cos(Tpc[Imin][Jmin]);
	*My = Myc + Myh /*+ Mlaty */ + Mfricy - Fc*Zpc[Imin][Jmin]*sin(Tpc[Imin][Jmin]);
	//*Mz = Mfricz;
	
	printf(" Fx	Mx	Flatx	Fxc	Mxc	Fxh Mxh  Mfricx Ffric \n");
	printf("%4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e \n", (*Fx)*Lambda*Rp*Rp, (*Mx)*Lambda*Rp*Rp*Rp,  Flatx*Lambda*Rp*Rp,  Fxc*Lambda*Rp*Rp, Mxc*Lambda*Rp*Rp*Rp, Fxh*Lambda*Rp*Rp, Mxh*Lambda*Rp*Rp*Rp, Mfricx*Lambda*Rp*Rp*Rp, Ffric*Lambda*Rp*Rp );
	printf(" Fy	My	Flaty	Fyc	Myc	Fxh	Myh  Mfricy Fz \n");
	printf("%4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e \n", (*Fy)*Lambda*Rp*Rp, (*My)*Lambda*Rp*Rp*Rp,  Flaty*Lambda*Rp*Rp,  Fyc*Lambda*Rp*Rp, Myc*Lambda*Rp*Rp*Rp, Fyh*Lambda*Rp*Rp, Myh*Lambda*Rp*Rp*Rp, Mfricy*Lambda*Rp*Rp*Rp, (*Fz)*Lambda*Rp*Rp );

	//	Dimensional ParameteZs
	FxcD = Fxc * (Lambda*pow(Rp,2));
	FycD = Fyc * (Lambda*pow(Rp,2));
	MxcD = Mxc * (Lambda*pow(Rp,3));
	MycD = Myc * (Lambda*pow(Rp,3));	
	FxhD = Fxh * (Lambda*pow(Rp,2));
	FyhD = Fyh * (Lambda*pow(Rp,2));
	MxhD = Mxh  * (Lambda*pow(Rp,3));
	MyhD = Myh * (Lambda*pow(Rp,3));
	FlatxD = Flatx*Lambda*Rp*Rp;
	FlatyD = Flaty*Lambda*Rp*Rp;
	MlatxD = Mlatx*Lambda*Rp*Rp*Rp;
	MlatyD = Mlaty*Lambda*Rp*Rp*Rp;
	FfricD = Ffric*Lambda*Rp*Rp;
	MfricxD = Mfricx*Lambda*Rp*Rp*Rp;
	MfricyD = Mfricx*Lambda*Rp*Rp*Rp;
	FxD = (*Fx) * (Lambda*pow(Rp,2));
	FyD = (*Fy) * (Lambda*pow(Rp,2));
	FzD = (*Fz)*Lambda*Rp*Rp;
	MxD = (*Mx)*Lambda*Rp*Rp*Rp;
	MyD= (*My)*Lambda*Rp*Rp*Rp;
	MzD = (*Mz)*Lambda*Rp*Rp*Rp;
}
void FMc(double phi, float **Zpc, double *Pc, float *Fxc, float *Fyc, float *Mxc, float *Myc, double *ym)
{
	int i,j;
	float Abc, FxcD, FycD, MxcD, MycD;	// Cylinder pressure closing area

	Abc = (Ap/pow(Rbo,2) - ((zB3*zB3-zB2*zB2)*(thetaK-2*rB/R)/2+Pi*rB*rB/pow(Rbo,2))); 
	
	//	Nondimensional force by cylinder pressure	
	*Fxc = Pc[1]*Ap/Rp/Rp*sqrt(ym[3]*ym[3]+ym[4]*ym[4])*hp/Rp*cos(-atan2(tan(ym[3]),sin(ym[4])));
	*Fyc = Pc[1]*Ap/Rp/Rp*sqrt(ym[3]*ym[3]+ym[4]*ym[4])*hp/Rp*sin(-atan2(tan(ym[3]),sin(ym[4])));

	//	Dimensional force
	FxcD = (*Fxc) * (Lambda*pow(Rp,2));
	FycD = (*Fyc) * (Lambda*pow(Rp,2));

	//	Nondimensional moment by cylinder pressure	
	*Mxc = -(Lpj-Lpg)/Rp*Pc[1]*Ap/Rp/Rp*tan(beta);
	*Myc = 0.0;

	//	Dimensional moment
	MxcD = (*Mxc) * (Lambda*pow(Rp,3));
	MycD = (*Myc) * (Lambda*pow(Rp,3));
}
void FMh(double phi, float **Zpc, float **Zn, float **Zs, float **Tpc, float **dT, float **Acv, float **P, float *Fxh, float *Fyh, float *Mxh, float *Myh)
{
	int i,j;	

	//float temp=0.0;
	//	for(j=1; j<=NNBt; j+=NNBt/NNBt) 
	//		for(i=1; i<=NNBz; i++) temp = temp + P[i][j]; 
	//			//if(P[i][j] > 1.0) 
	//printf("%e  \n", temp);

	//	Nondimensional lifing force
	*Fxh = 0.0;
	*Fyh = 0.0;
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBz; i++){
			*Fxh = (*Fxh) - P[i][j]*Acv[i][j]*sin(Tpc[i][j]);
			*Fyh = (*Fyh) - P[i][j]*Acv[i][j]*cos(Tpc[i][j]);
			//printf("%d %d %e \n", i, j, *Fy);
		}
	}
	//printf("%e \n", *Fy);
	//	Dimensional lifing force
	FxhD = (*Fxh) * (Lambda*pow(Rp,2));
	FyhD = (*Fyh) * (Lambda*pow(Rp,2));
 
	*Mxh = 0.0;
	*Myh = 0.0;
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBz; i++){
			*Mxh = (*Mxh) + P[i][j]*Acv[i][j]*Zpc[i][j]*cos(Tpc[i][j]);
			*Myh = (*Myh) - P[i][j]*Acv[i][j]*Zpc[i][j]*sin(Tpc[i][j]);
		}
	}

	//	Dimensional lifting moment
	MxhD = (*Mxh) * (Lambda*pow(Rp,3));
	MyhD = (*Myh) * (Lambda*pow(Rp,3));

}

/* Piston friction */
//void FMfriction(double phi, double *Pc, float *Ffric, float *Mfricx, float *Mfricy)
//{
//	int i,j;
//	double *phiD, *Lpc, *Vp;
//
//	phiD = dvector(1,N);
//	Lpc = dvector(1,N);
//	Vp = dvector(1,N);
//
//	// Non-dimensional F-M
//	*Ffric = 0.0;
//	*Mfricx = 0.0;
//	*Mfricy = 0.0;
//	for(i=1;i<=N;i++){
//		phiD[i] = phi + 2*Pi/N*(i-1);
//		if(phiD[i] > 2*Pi) phiD[i] = phiD[i] - 2*Pi;
//		Lpc[i] = Lc - Ld - R*(tan(beta)*(1+cos(phiD[i]))-tan(beta2)/cos(beta)*sin(phiD[i]));
//		Vp[i] = -R*w*(tan(beta)*sin(phiD[i])+tan(beta2)/cos(beta)*cos(phiD[i]));
//		*Ffric = *Ffric + ((Pc[i]*Lambda-Ph)*hp/Lpc[i]/2.0+eta*Vp[i]/hp)*2*Pi*(Rp+hp)*Lpc[i] /Lambda/pow(Rbo,2);
//		*Mfricx = *Mfricx + ((Pc[i]*Lambda-Ph)*hp/Lpc[i]/2.0+eta*Vp[i]/hp)*2*Pi*(Rp+hp)*Lpc[i]*R*cos(phiD[i]) /Lambda/pow(Rbo,3);
//		*Mfricy = *Mfricy - ((Pc[i]*Lambda-Ph)*hp/Lpc[i]/2.0+eta*Vp[i]/hp)*2*Pi*(Rp+hp)*Lpc[i]*R*sin(phiD[i]) /Lambda/pow(Rbo,3);
//	}
//
//	free_dvector(phiD,1,N);	
//	free_dvector(Lpc,1,N);
//	free_dvector(Vp,1,N);
//}
void FMlateral(double phi, double *Pc, float *Flatx, float *Flaty, float *Mlatx, float *Mlaty, float *Mlatz)
{
	int i,j;
	double *phiD, *Lpc, *Vp, *Ap;
	float Fpa, Fwpa, Ffws, z0x, z0y;

	phiD = dvector(1,N);
	Lpc = dvector(1,N);
	Vp = dvector(1,N);
	Ap = dvector(1,N);

	*Flatx = 0.0;
	*Flaty = 0.0;
	*Mlatx = 0.0;
	*Mlaty = 0.0;
	*Mlatz = 0.0;
	//for(i=1;i<=N;i++){
	//	phiD[i] = phi + 2*Pi/N*(i-1);
	//	while(phiD[i] > 2*Pi) phiD[i] = phiD[i] - 2*Pi;
	//	Lpc[i] = Lc - Ld - R*(tan(beta)*(1+cos(phiD[i]))-tan(beta2)/cos(beta)*sin(phiD[i]));
	//	Vp[i] = -R*w*(tan(beta)*sin(phiD[i])+tan(beta2)/cos(beta)*cos(phiD[i]));
	//	Ap[i] = -R*w*w*(tan(beta)*cos(phiD[i])-tan(beta2)/cos(beta)*sin(phiD[i]));
		//Fpa = (Pc[i]*Lambda*Pi*Rp*Rp - Mp/w/w*Lambda*Rbo*(Rp/hp)*Ap[i] + ((Pc[i]*Lambda-Ph)*hp/Lpc[i]/2.0-eta*Vp[i]/hp)*2*Pi*(Rp+hp)*Lpc[i])/Lambda/pow(Rbo,2);
		Fwpa = Mp*R/Rp*(Rp/hp); 
		//Ffws = eta*R*w/hs*Pi*(Rs2*Rs2-Rs1*Rs1)/Lambda/pow(Rbo,2);
		//z0y = (Fpa*tan(beta)+Fwpa*cos(phiD[i])+Ffws*sin(phiD[i]))/12/((Fpa*tan(beta)+Ffws*sin(phiD[i]))*(Lpj-Lpc[i]/2)+Fwpa*cos(phiD[i])*(Lpg-Lpc[i]/2))*pow(Lpc[i],2);
		//z0x = (Fwpa*sin(phiD[i])-Ffws*cos(phiD[i]))/12/((Ffws*cos(phiD[i]))*(Lpj-Lpc[i]/2)-Fwpa*sin(phiD[i])*(Lpg-Lpc[i]/2))*pow(Lpc[i],2);
		*Flatx = Fwpa*sin(phi)/*-Ffws*cos(phiD[i]))*/;					  
		*Flaty = /*(Fpa*tan(beta)+*/Fwpa*cos(phi)/*+Ffws*sin(phiD[i]))*/;
		//*Mlaty = *Mlaty + (  ((z0x+Lpc[i]/2)/3+Loc)*(z0x+Lpc[i]/2)*((Fwpa*sin(phiD[i])-Ffws*cos(phiD[i]))/Lpc[i]+6*((Ffws*cos(phiD[i]))*(Lpj-Lpc[i]/2)-Fwpa*sin(phiD[i])*(Lpg-Lpc[i]/2))/pow(Lpc[i],2))/2
		//					- (Lpc[i]*5/6+z0x/3+Loc)*(z0x-Lpc[i]/2)*((Fwpa*sin(phiD[i])-Ffws*cos(phiD[i]))/Lpc[i]-6*((Ffws*cos(phiD[i]))*(Lpj-Lpc[i]/2)-Fwpa*sin(phiD[i])*(Lpg-Lpc[i]/2))/pow(Lpc[i],2))/2  )/Rbo;
		//*Mlatx = *Mlatx + ( ((z0y+Lpc[i]/2)/3+Loc)*(z0y+Lpc[i]/2)*((Fpa*tan(beta)+Fwpa*cos(phiD[i])+Ffws*sin(phiD[i]))/Lpc[i]+6*((Fpa*tan(beta)+Ffws*sin(phiD[i]))*(Lpj-Lpc[i]/2)+Fwpa*cos(phiD[i])*(Lpg-Lpc[i]/2))/pow(Lpc[i],2))/2
		//					- (Lpc[i]*5/6+z0y/3+Loc)*(z0y-Lpc[i]/2)*((Fpa*tan(beta)+Fwpa*cos(phiD[i])+Ffws*sin(phiD[i]))/Lpc[i]-6*((Fpa*tan(beta)+Ffws*sin(phiD[i]))*(Lpj-Lpc[i]/2)+Fwpa*cos(phiD[i])*(Lpg-Lpc[i]/2))/pow(Lpc[i],2))/2 )/Rbo;		
		//*Mlatz = *Mlatz + (R*cos(phi)*(Fpa*tan(beta)+Fwpa*cos(phiD[i])+Ffws*sin(phiD[i]))- R*sin(phi)*(Fwpa*sin(phiD[i])-Ffws*cos(phiD[i])))/Rbo ;

	//}
	//printf("++++++++++++ %e %e %e %e \n", phi, *Flaty, *Mlatx, *Mlaty);

	free_dvector(phiD,1,N);	
	free_dvector(Lpc,1,N);
	free_dvector(Vp,1,N);
	free_dvector(Ap,1,N);
}