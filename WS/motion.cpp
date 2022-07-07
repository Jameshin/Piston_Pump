#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

extern int Nphi, NBr1, NBr2, NBr3, NBr4, NBr5, NBr, NNBr, NNBt, Nch, Imin, Jmin;
extern int NVt1, NVt2, NVt3, NVt4, NVt5, NVt6, NVr1, NVr2, NVr3, NVr4;
extern double dNBt, thetach;
extern double R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, revf, phif, phiout_c, hm, phiout_cy, phiout_dy, phiout_vb, phiout_FM, phiout_LK, w;
extern double V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
extern double integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, rB1, rB2, rB3, rB4, rB5, rB6, Fs, hv0, theta0, psi0, dhv0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
extern double Lambda, Loc, Lpj, Lpg, Lbj, Ls, Rbo;
extern float FpressD, MpressxD, MpressyD, FliftD, MliftxD, MliftyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, FzD, MxD, MyD, MzD, leakageRinD, leakageRoutD, frictorqueD, friclossD, Hmin;
extern double gamH, zetaH;
/* Solve equations of motion */
/* Runge-Kutta 4th */
void motion_RK4(double phi, double h, double *ym, double gamma, double zeta, double dgamma, double dzeta, double **Hvb, double **Hn, double **Hs, double **Hw, double **He, double **Hvb_gro, 
				float **Rvb, float **Tvb, float **dTvb, float **dRvb, float **Rvbn, float **Rvbs, float **Rvbw, float **Rvbe, float **Tvbn, float **Tvbs, float **Tvbw, float **Tvbe, 
				float **dNvb, float **dSvb, float **dWvb, float **dEvb, float **aN, float **aS, float **aW, float **aE, float **aC, float **B, float **Pvb, double *Pcn, float **Acv, 
				void (*configuration)(double , double *, double *, double *, double *, double *), 
				void (*filmthickness)(double , double *, double **, double **, double **, double **, double **, double **, float **, float **, double , double ), 
				void (*descretize)(double , double *, double , double , double , double , float **, float **, float **, double **, double **, double **, double **, float **, float **, 
									float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **),
				void (*solver)(double, double *, float **, float **, float **, float **, float **, float **, float **), 
				void (*FM)(double , double *, float **, float **, float **, float **, float **, float **, float **, float *, float *, float *, float *, float *, float *, double *) )
{
	void Derivsm(double , double *, double *, double gamma, double zeta, double dgamma, double dzeta, double **Hvb, double **Hn, double **Hs, double **Hw, double **He, double **Hvb_gro, 
				float **Rvb, float **Tvb, float **dTvb, float **dRvb, float **Rvbn, float **Rvbs, float **Rvbw, float **Rvbe, float **Tvbn, float **Tvbs, float **Tvbw, float **Tvbe, 
				float **dNvb, float **dSvb, float **dWvb, float **dEvb, float **aN, float **aS, float **aW, float **aE, float **aC, float **B, float **Pvb, double *Pcn, float **Acv,
				void (*configuration)(double , double *, double *, double *, double *, double *), 
				void (*filmthickness)(double , double *, double **, double **, double **, double **, double **, double **, float **, float **, double , double ), 
				void (*descretize)(double , double *, double , double , double , double , float **, float **, float **, double **, double **, double **, double **, float **, float **, 
									float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **),
				void (*solver)(double, double *, float **, float **, float **, float **, float **, float **, float **), 
				void (*FM)(double , double *, float **, float **, float **, float **, float **, float **, float **, float *, float *, float *, float *, float *, float *, double *));
	//void Derivsm_c(double , double *, double *, double gamma, double zeta, double dgamma, double dzeta, float **Hvb, float **Hn, float **Hs, float **Hw, float **He, float **Hvb_gro, 
	//			float **Rvb, float **Tvb, float **dTvb, float **dRvb, float **Rvbn, float **Rvbs, float **Rvbw, float **Rvbe, float **Tvbn, float **Tvbs, float **Tvbw, float **Tvbe, 
	//			float **dNvb, float **dSvb, float **dWvb, float **dEvb, float **aN, float **aS, float **aW, float **aE, float **aC, float **B, float **Pvb, double *Pcn, float **Acv,
	//			void (*configuration)(double , double *, double *, double *, double *, double *), 
	//			void (*filmthickness)(double , double *, float **, float **, float **, float **, float **, float **, float **, float **, double , double ), 
	//			void (*descretize)(double , double *, double , double , double , double , float **, float **, float **, float **, float **, float **, float **, float **, float **, 
	//								float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **), 
	//			void (*solver)(double, double *, float **, float **, float **, float **, float **, float **, float **), 
	//			void (*FM)(double , double *, float **, float **, float **, float **, float **, float **, float **, float *, float *, float *, float *, float *, float *, double *));

	int i,j;
	double *K1, *K2, *K3, *K4;
	double *slope, *ymR;

	K1 = dvector(1,2*DOF);
	K2 = dvector(1,2*DOF);
	K3 = dvector(1,2*DOF);
	K4 = dvector(1,2*DOF);
	slope = dvector(1,2*DOF);
	ymR = dvector(1,2*DOF);

	//if(phi*180/Pi >45) h = 0.001*Pi/180;
	//if(phi*180/Pi >50) h = 0.0001*Pi/180;

	printf("****** hm = %lf \n", h);
	//gamma = acos(cos(ym[2])*cos(ym[3]))*Rbo/hv0;

	//if((ym[1]-acos(cos(ym[2])*cos(ym[3]))*Rbo/hv0)*hv0 > hvbmin){
		Derivsm(phi, ym, K1, gamma, zeta, dgamma, dzeta, Hvb, Hn, Hs, Hw, He, Hvb_gro, 
				Rvb, Tvb, dTvb, dRvb, Rvbn, Rvbs, Rvbw, Rvbe, Tvbn, Tvbs, Tvbw, Tvbe, 
				dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Pvb, Pcn, Acv, 
				configuration, filmthickness, descretize, solver, FM);
		//for(i=1; i<=2*DOF; i++) ymR[i] = ym[i] + K1[i]*h/2.;
		//Derivsm(phi+h/2, ymR, K2, gamma, zeta, dgamma, dzeta, Hvb, Hn, Hs, Hw, He, Hvb_gro, 
		//		Rvb, Tvb, dTvb, dRvb, Rvbn, Rvbs, Rvbw, Rvbe, Tvbn, Tvbs, Tvbw, Tvbe, 
		//		dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Pvb, Pcn, Acv, 
		//		configuration, filmthickness, descretize, solver, FM);
		//for(i=1; i<=2*DOF; i++) ymR[i] = ym[i] + K2[i]*h/2.;
		//Derivsm(phi+h/2, ymR, K3, gamma, zeta, dgamma, dzeta, Hvb, Hn, Hs, Hw, He, Hvb_gro, 
		//		Rvb, Tvb, dTvb, dRvb, Rvbn, Rvbs, Rvbw, Rvbe, Tvbn, Tvbs, Tvbw, Tvbe, 
		//		dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Pvb, Pcn, Acv, 
		//		configuration, filmthickness, descretize, solver, FM);
		//for(i=1; i<=2*DOF; i++) ymR[i] = ym[i] + K3[i]*h;	
		//Derivsm(phi+h, ymR, K4, gamma, zeta, dgamma, dzeta, Hvb, Hn, Hs, Hw, He, Hvb_gro, 
		//		Rvb, Tvb, dTvb, dRvb, Rvbn, Rvbs, Rvbw, Rvbe, Tvbn, Tvbs, Tvbw, Tvbe, 
		//		dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Pvb, Pcn, Acv, 
		//		configuration, filmthickness, descretize, solver, FM);
		for(i=1; i<=2*DOF; i++){
			//slope[i] = (K1[i]+2.0*(K2[i]+K3[i])+K4[i])/6.0;
			//ym[i] = ym[i] + slope[i]*h;
			ym[i] = ym[i] + K1[i]*h;
		}	
	//}
	//else{		
	//	ym[1] = acos(cos(ym[2])*cos(ym[3]))*Rbo/hv0 + hvbmin/hv0;  //initial conf.			
	//	ym[4] = (sin(ym[2])*cos(ym[3])*ym[5]+cos(ym[2])*sin(ym[3])*ym[6])/acos(cos(ym[2])*cos(ym[3]))*pow(Rbo/hv0,1);
	//	Derivsm_c(phi, ym, K1, gamma, zeta, dgamma, dzeta, Hvb, Hn, Hs, Hw, He, Hvb_gro, 
	//			Rvb, Tvb, dTvb, dRvb, Rvbn, Rvbs, Rvbw, Rvbe, Tvbn, Tvbs, Tvbw, Tvbe, 
	//			dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Pvb, Pcn, Acv, 
	//			configuration, filmthickness, descretize, solver, FM);
	//	for(i=1; i<=2*DOF; i++) ymR[i] = ym[i] + K1[i]*h/2.;
	//	Derivsm_c(phi+h/2, ymR, K2, gamma, zeta, dgamma, dzeta, Hvb, Hn, Hs, Hw, He, Hvb_gro, 
	//			Rvb, Tvb, dTvb, dRvb, Rvbn, Rvbs, Rvbw, Rvbe, Tvbn, Tvbs, Tvbw, Tvbe, 
	//			dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Pvb, Pcn, Acv, 
	//			configuration, filmthickness, descretize, solver, FM);
	//	for(i=1; i<=2*DOF; i++) ymR[i] = ym[i] + K2[i]*h/2.;
	//	Derivsm_c(phi+h/2, ymR, K3, gamma, zeta, dgamma, dzeta, Hvb, Hn, Hs, Hw, He, Hvb_gro, 
	//			Rvb, Tvb, dTvb, dRvb, Rvbn, Rvbs, Rvbw, Rvbe, Tvbn, Tvbs, Tvbw, Tvbe, 
	//			dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Pvb, Pcn, Acv, 
	//			configuration, filmthickness, descretize, solver, FM);
	//	for(i=1; i<=2*DOF; i++) ymR[i] = ym[i] + K3[i]*h;	
	//	Derivsm_c(phi+h, ymR, K4, gamma, zeta, dgamma, dzeta, Hvb, Hn, Hs, Hw, He, Hvb_gro, 
	//			Rvb, Tvb, dTvb, dRvb, Rvbn, Rvbs, Rvbw, Rvbe, Tvbn, Tvbs, Tvbw, Tvbe, 
	//			dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B, Pvb, Pcn, Acv, 
	//			configuration, filmthickness, descretize, solver, FM);
	//	for(i=1; i<=2*DOF; i++){
	//		slope[i] = (K1[i]+2.0*(K2[i]+K3[i])+K4[i])/6.0;
	//		ym[i] = ym[i] + slope[i]*h;
	//	}
	//			
	//	ym[1] = acos(cos(ym[2])*cos(ym[3]))*Rbo/hv0 + hvbmin/hv0;  //initial conf.			
	//	ym[4] = (sin(ym[2])*cos(ym[3])*ym[5]+cos(ym[2])*sin(ym[3])*ym[6])/acos(cos(ym[2])*cos(ym[3]))*pow(Rbo/hv0,1);
	//	//ym[2] = ym[2]*0.99999;
	//	//ym[3] = ym[3]*0.99999;
	//}

	free_dvector(K1,1,2*DOF);
	free_dvector(K2,1,2*DOF);
	free_dvector(K3,1,2*DOF);
	free_dvector(K4,1,2*DOF);
	free_dvector(slope,1,2*DOF);
	free_dvector(ymR,1,2*DOF);
}
/* phi, y, dy */
void Derivsm( double phi, double *ym, double *dy, double gamma, double zeta, double dgamma, double dzeta, double **Hvb, double **Hn, double **Hs, double **Hw, double **He, double **Hvb_gro, 
				float **Rvb, float **Tvb, float **dTvb, float **dRvb, float **Rvbn, float **Rvbs, float **Rvbw, float **Rvbe, float **Tvbn, float **Tvbs, float **Tvbw, float **Tvbe, 
				float **dNvb, float **dSvb, float **dWvb, float **dEvb, float **aN, float **aS, float **aW, float **aE, float **aC, float **B, float **Pvb, double *Pcn, float **Acv,
				void (*configuration)(double , double *, double *, double *, double *, double *), 
				void (*filmthickness)(double , double *, double **, double **, double **, double **, double **, double **, float **, float **, double , double ), 
				void (*descretize)(double , double *, double , double , double , double , float **, float **, float **, double **, double **, double **, double **, float **, float **, 
									float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **),
				void (*solver)(double, double *, float **, float **, float **, float **, float **, float **, float **), 
				void (*FM)(double , double *, float **, float **, float **, float **, float **, float **, float **, float *, float *, float *, float *, float *, float *, double *)) 
{
	float Fx, Fy, Fz, Mx, My, Mz;
	float Ru, Lu, Mu;
	double Ibxx, Ibyy, Ibzz, Ibxy, Ibzx, Ibyz, dIbxx, dIbyy, dIbzz, dIbxy, dIbzx, dIbyz;

	(*configuration)(phi, ym, &gamma, &zeta, &dgamma, &dzeta);
	(*filmthickness)(phi, ym, Hvb, Hn, Hs, Hw, He, Hvb_gro, Rvb, Tvb, gamma, zeta);
	(*descretize)(phi, ym, gamma, zeta, dgamma, dzeta, Rvb, dTvb, dRvb, Hn, Hs, Hw, He, Rvbn, Rvbs, Rvbw, Rvbe, 
				Tvbn, Tvbs, Tvbw, Tvbe, dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B);
	(*solver)(phi, Pcn, Pvb, aN, aS, aW, aE, aC, B);
	(*FM)(phi, Pcn, Rvb, Pvb, Rvbn, Rvbs, Tvb, dTvb, Acv, &Fx, &Fy, &Fz, &Mx, &My, &Mz, ym);

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

 
	dy[1] = ym[4];
	dy[2] = ym[5];
	dy[3] = ym[6]/cos(ym[2]);
	//dy[4] = Fz*cos(ym[2])*cos(ym[3])/(Mb+Mu);
	//dy[5] = ( Ibyy*Mx - Ibxy*My - (Ibzx*Ibyy-Ibxy*Ibyz)*(sin(ym[2])*cos(ym[3])*ym[5]+cos(ym[2])*sin(ym[3])*ym[6]) 
	//		- (dIbxy*Ibyy-dIbyy*Ibxy)*ym[6]*cos(ym[2]) + (dIbzx*Ibyy-dIbyz*Ibxy)*cos(ym[2])*cos(ym[3]) 
	//		- (Ibzx*Ibyy+Ibyz*Ibxy+Ibxy*Ibyy*tan(ym[2])+Ibxx*Ibxy*tan(ym[2]))*ym[5]*ym[6]*cos(ym[2]) 
	//		- Ibyz*Ibyy*pow(ym[6]*cos(ym[2]),2) - Ibzx*Ibxy*pow(ym[5],2) - (Ibyy*Ibyy+Ibxy*Ibxy)*pow(ym[6]*cos(ym[2]),2)*tan(ym[2])
	//		+ Ibzz*Ibyy*ym[6]*cos(ym[2])*cos(ym[2])*cos(ym[3]) + Ibzz*Ibxy*ym[5]*cos(ym[2])*cos(ym[3])
	//		+ (Ibyz*Ibyy+Ibzx*Ibxy)*cos(ym[2])*cos(ym[3])*ym[6]*cos(ym[2])*tan(ym[2]) 
	//		) / (Ibxx*Ibyy-Ibxy*Ibxy);
	//dy[6] = ( ( Ibxx*My - Ibxy*Mx - (Ibyz*Ibxx-Ibxy*Ibzx)*(sin(ym[2])*cos(ym[3])*ym[5]+cos(ym[2])*sin(ym[3])*ym[6]) 
	//		- (dIbxy*Ibxx-dIbxx*Ibxy)*ym[6]*cos(ym[2]) + (dIbyz*Ibxx-dIbzx*Ibxy)*cos(ym[2])*cos(ym[3]) 
	//		+ (Ibxx*Ibxx*tan(ym[2])+Ibxy*Ibxy*tan(ym[2])+Ibyz*Ibxx+Ibzx*Ibxy)*ym[5]*ym[6]*cos(ym[2]) 
	//		+ Ibyz*Ibxy*pow(ym[6]*cos(ym[2]),2) + Ibzx*Ibxx*pow(ym[5],2) + (Ibxy*Ibxx-Ibyy*Ibxy)*pow(ym[6]*cos(ym[2]),2)*tan(ym[2])
	//		- Ibzz*Ibxy*ym[6]*cos(ym[2])*cos(ym[2])*cos(ym[3]) - Ibzz*Ibxx*ym[5]*cos(ym[2])*cos(ym[3])
	//		- (Ibzx*Ibxx+Ibyz*Ibxy)*cos(ym[2])*cos(ym[3])*ym[6]*cos(ym[2])*tan(ym[2]) 
	//		) / (Ibxx*Ibyy-Ibxy*Ibxy) + ym[6]*sin(ym[2])*ym[5]
	//		)/cos(ym[2]);
	dy[4] = Fz*cos(ym[2])*cos(ym[3])/Mb;
	dy[5] = (Mx-Lbj/Rbo*Mb*Lbj/Rbo*gamma*sin(zeta)-Lbj/Rbo*acos(cos(ym[2])*cos(ym[3]))*sin(zeta)*Fz)/Ibt-pow(ym[6],2)*tan(ym[2])+Iba/Ibt*ym[6]*cos(ym[2])*cos(ym[3]);
	dy[6] = (My-Lbj/Rbo*Mb*Lbj/Rbo*gamma*cos(zeta)-Lbj/Rbo*acos(cos(ym[2])*cos(ym[3]))*cos(zeta)*Fz)/Ibt+ym[5]*ym[6]*tan(ym[2])-Iba/Ibt*ym[5]*cos(ym[2])*cos(ym[3]);	

	//printf("%3.3e %3.3e %3.3e %3.3e %3.3e %3.3e \n", dy[1], dy[2], dy[3], dy[4], dy[5], dy[6]);
}
/* phi, y, dy incase of boundary lubrication */
//void Derivsm_c( double phi, double *ym, double *dy, double gamma, double zeta, double dgamma, double dzeta, float **Hvb, float **Hn, float **Hs, float **Hw, float **He, float **Hvb_gro, 
//				float **Rvb, float **Tvb, float **dTvb, float **dRvb, float **Rvbn, float **Rvbs, float **Rvbw, float **Rvbe, float **Tvbn, float **Tvbs, float **Tvbw, float **Tvbe, 
//				float **dNvb, float **dSvb, float **dWvb, float **dEvb, float **aN, float **aS, float **aW, float **aE, float **aC, float **B, float **Pvb, double *Pcn, float **Acv,
//				void (*configuration)(double , double *, double *, double *, double *, double *), 
//				void (*filmthickness)(double , double *, float **, float **, float **, float **, float **, float **, float **, float **, double , double ), 
//				void (*descretize)(double , double *, double , double , double , double , float **, float **, float **, float **, float **, float **, float **, float **, float **, 
//									float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **), 
//				void (*solver)(double, double *, float **, float **, float **, float **, float **, float **, float **), 
//				void (*FM)(double , double *, float **, float **, float **, float **, float **, float **, float **, float *, float *, float *, float *, float *, float *, double *)) 
//{
//	float Fx, Fy, Fz, Mx, My, Mz;
//	float Ax, Bx, Kx, Ay, By, Ky;
//
//	(*configuration)(phi, ym, &gamma, &zeta, &dgamma, &dzeta);
//	(*filmthickness)(phi, ym, Hvb, Hn, Hs, Hw, He, Hvb_gro, Rvb, Tvb, gamma, zeta);
//	(*descretize)(phi, ym, gamma, zeta, dgamma, dzeta, Rvb, dTvb, dRvb, Hn, Hs, Hw, He, Rvbn, Rvbs, Rvbw, Rvbe, 
//				Tvbn, Tvbs, Tvbw, Tvbe, dNvb, dSvb, dWvb, dEvb, aN, aS, aW, aE, aC, B);
//	(*solver)(phi, Pcn, Pvb, aN, aS, aW, aE, aC, B);
//	(*FM)(phi, Pcn, Rvb, Pvb, Rvbn, Rvbs, Tvb, dTvb, Acv, &Fx, &Fy, &Fz, &Mx, &My, &Mz, ym);
//
//	Ax = Ibt - sin(zeta)*Mb*sin(ym[2])*cos(ym[3])/gamma*pow(Rbo/hv0,2)*w;
//	Bx = -sin(zeta)*Mb*cos(ym[2])*sin(ym[3])/gamma*pow(Rbo/hv0,2)*w;
//	Kx = Mx + sin(zeta)*(Mb*(-dgamma*dgamma*pow(hv0/Rbo,2)*cos(ym[2])*cos(ym[3])+cos(ym[2])*cos(ym[3])*(ym[5]*ym[5]+ym[6]*ym[6])-2*sin(ym[2])*sin(ym[3]))/gamma*pow(Rbo/hv0,2)*w-Fz)
//		-Ibt*pow(ym[6]*cos(ym[2]),2)*ym[2] - Iba*ym[6]*cos(ym[2])*cos(ym[2])*cos(ym[3]);
//	Ay = Ibt*cos(ym[2]) + cos(zeta)*Mb*cos(ym[2])*sin(ym[3])/gamma*pow(Rbo/hv0,2)*w;
//	By = cos(zeta)*Mb*sin(ym[2])*cos(ym[3])/gamma*pow(Rbo/hv0,2)*w;
//	Ky = My - cos(zeta)*(Mb*(-dgamma*dgamma*pow(hv0/Rbo,2)*cos(ym[2])*cos(ym[3])+cos(ym[2])*cos(ym[3])*(ym[5]*ym[5]+ym[6]*ym[6])-2*sin(ym[2])*sin(ym[3]))/gamma*pow(Rbo/hv0,2)*w-Fz)
//		+Ibt*ym[6]*(sin(ym[2])+cos(ym[2]))*ym[5] + Iba*ym[5]*cos(ym[2])*cos(ym[3]);		
//	dy[1] = 0.0/*ym[4]*/;
//	dy[2] = ym[5];
//	dy[3] = ym[6];
//	dy[5] = (Kx*By/Bx-Ky)/(Ax*By/Bx-Ay);
//	dy[6] = (Kx*Ay/Ax-Ky)/(Bx*Ay/Ax-By);
//	dy[4] = 0.0/*(-dgamma*dgamma*pow(hv0/Rbo,2)*cos(ym[2])*cos(ym[3])+cos(ym[2])*cos(ym[3])*(ym[5]*ym[5]+ym[6]*ym[6])-2*sin(ym[2])*sin(ym[3])+sin(ym[2])*cos(ym[3])*dy[2]+cos(ym[2])*sin(ym[3])*dy[3])/gamma*pow(Rbo/hv0,2)*w*/;
//
//	printf("%3.3e %3.3e %3.3e %3.3e %3.3e %3.3e \n", dy[1], dy[2], dy[3], dy[4], dy[5], dy[6]);
//}
void configuration(double phi, double *ym, double *gamma, double *zeta, double *dgamma, double *dzeta)
{	
	// nondimension
	*gamma = acos(cos(ym[2])*cos(ym[3]))*Rbo/hv0;
	*zeta = -atan2(tan(ym[2]),sin(ym[3]));
	*dgamma = (sin(ym[2])*cos(ym[3])*ym[5]+cos(ym[2])*sin(ym[3])*ym[6])/sin((*gamma)*hv0/Rbo)*pow(Rbo/hv0,1);
	*dzeta = ((tan(ym[2])*cos(ym[3])*ym[6]-sin(ym[3])/cos(ym[2])/cos(ym[2])*ym[5])*cos(*zeta)*cos(*zeta)/sin(ym[3])/sin(ym[3]));

	if(*gamma == 0. || *zeta == 0.){ 
		printf("Gamma or Zeta = 0.0");
		*gamma = TINY;
		*dgamma = (sin(ym[2])*cos(ym[3])*ym[5]+cos(ym[2])*sin(ym[3])*ym[6])/(*gamma)*pow(Rbo/hv0,2);
		//exit(1);
	}
	//if(*zeta < 0) *zeta = *zeta + 2*Pi;
	printf("\nphi = %e \n%4.2e %4.2e %6.5e %4.3e %4.2e %4.2e %4.2e %4.2e %e %d %d \n\n", phi*180/Pi, *gamma*hv0/Rbo*180/Pi, *zeta*180/Pi, ym[1]*hv0, ym[2]*180/Pi, ym[3]*180/Pi, ym[4]*hv0*w, ym[5]*w*180/Pi, ym[6]*w*180/Pi, (/*ym[1]-(*gamma)*/Hmin)*hv0, Imin, Jmin);
	
}
/* Force & Moment of the Cylinder Barrel */
void FM(double phi, double *Pc, float **Rvb, float **Pvb, float **Rvbn, float **Rvbs, float **Tvb, float **dTvb, float **Acv, float *Fx, float *Fy, float *Fz, float *Mx, float *My, float *Mz, double *ym)
{
	void FMpress(double phi, float **Rvb, double *Pc, float *Fpress, float *Mpressx, float *Mpressy);
	void FMlift(double phi, float **R, float **Rn, float **Rs, float **Tvb, float **dT, float **Acv, float **P, float *Flift, float *Mliftx, float *Mlifty);
	void FMfriction(double phi, double *Pc, float *Ffric, float *Mfricx, float *Mfricy);
	void FMlateral(double phi, double *Pc, float *Flatx, float *Flaty, float *Mlatx, float *Mlaty, float *Mlatz);

	float Fpress, Mpressx, Mpressy, Flift, Mliftx, Mlifty, Ffric, Mfricx, Mfricy;
	float Flatx, Flaty, Mlatx, Mlaty, Mlatz;

	FMpress(phi, Rvb, Pc, &Fpress, &Mpressx, &Mpressy);
	FMlift(phi, Rvb, Rvbn, Rvbs, Tvb, dTvb, Acv, Pvb, &Flift, &Mliftx, &Mlifty);
	FMfriction(phi, Pc, &Ffric, &Mfricx, &Mfricy);
	FMlateral(phi, Pc, &Flatx, &Flaty, &Mlatx, &Mlaty, &Mlatz);

	*Fx = Flatx/**cos(phi) - Flaty*sin(phi)*/;
	*Fy = Flaty/**cos(phi) + Flatx*sin(phi)*/;
	//if(Hmin < hvbmin/hv0) *Fz = Fpress + Flift + Ffric + Kc*(hvbmin/hv0-Hmin);
	/*else*/ *Fz = Fpress + Flift + Ffric;
	//if(Hmin < hvbmin/hv0){
	//	*Mx = (Mpressx + Mliftx + Rvb[Imin][Jmin]*cos(Tvb[Imin][Jmin])*Kc*(hvbmin/hv0-Hmin))*cos(phi) + (Mpressy + Mlifty + Rvb[Imin][Jmin]*cos(Tvb[Imin][Jmin])*Kc*(hvbmin/hv0-Hmin))*sin(phi) + Mlatx + Mfricx + (*Fy)*Ls - Lbj/Rbo*Mb*Lbj/Rbo*acos(cos(ym[2])*cos(ym[3]))*sin(-atan2(tan(ym[2]),sin(ym[3])))*cos(ym[2])*cos(ym[3]) + Lbj*acos(cos(ym[2])*cos(ym[3]))*sin(-atan2(tan(ym[2]),sin(ym[3])))*(*Fz)/Rbo;
	//	*My = (Mpressy + Mlifty + Rvb[Imin][Jmin]*cos(Tvb[Imin][Jmin])*Kc*(hvbmin/hv0-Hmin))*cos(phi) - (Mpressx + Mliftx + Rvb[Imin][Jmin]*cos(Tvb[Imin][Jmin])*Kc*(hvbmin/hv0-Hmin))*sin(phi) + Mlaty + Mfricy - (*Fx)*Ls + Lbj/Rbo*Mb*Lbj/Rbo*acos(cos(ym[2])*cos(ym[3]))*cos(-atan2(tan(ym[2]),sin(ym[3])))*cos(ym[2])*cos(ym[3]) - Lbj*acos(cos(ym[2])*cos(ym[3]))*cos(-atan2(tan(ym[2]),sin(ym[3])))*(*Fz)/Rbo;
	//}
	//else{
		*Mx = (Mpressx + Mliftx)*cos(phi) + (Mpressy + Mlifty)*sin(phi) + Mlatx + Mfricx + (*Fy)*Ls /*- Lbj/Rbo*Mb*Lbj/Rbo*acos(cos(ym[2])*cos(ym[3]))*(Rbo/hv0)*sin(-atan2(tan(ym[2]),sin(ym[3])))*cos(ym[2])*cos(ym[3]) + Lbj*acos(cos(ym[2])*cos(ym[3]))*sin(-atan2(tan(ym[2]),sin(ym[3])))*(*Fz)/Rbo*/;
		*My = (Mpressy + Mlifty)*cos(phi) - (Mpressx + Mliftx)*sin(phi) + Mlaty + Mfricy - (*Fx)*Ls /*+ Lbj/Rbo*Mb*Lbj/Rbo*acos(cos(ym[2])*cos(ym[3]))*(Rbo/hv0)*cos(-atan2(tan(ym[2]),sin(ym[3])))*cos(ym[2])*cos(ym[3]) - Lbj*acos(cos(ym[2])*cos(ym[3]))*cos(-atan2(tan(ym[2]),sin(ym[3])))*(*Fz)/Rbo*/;
	//}
	*Mz = Mlatz;
	
	printf(" Ffric	Flift	Flatx	Mpressx	Mliftx Mlateralx Mfricx	Mx \n");
	printf("%e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e  %4.2e \n", Ffric*Lambda*Rbo*Rbo, Flift*Lambda*Rbo*Rbo, Flatx*Lambda*Rbo*Rbo, (Mpressx*cos(phi) + Mpressy*sin(phi))*Lambda*Rbo*Rbo*Rbo, (Mliftx*cos(phi) + Mlifty*sin(phi))*Lambda*Rbo*Rbo*Rbo, Mlatx*Lambda*Rbo*Rbo*Rbo, Mfricx*Lambda*Rbo*Rbo*Rbo, (*Mx)*Lambda*Rbo*Rbo*Rbo);
	printf("  Fz	Fpress	Flaty	Mpressy	Mlifty Mlateraly	Mfricy	My \n");
	printf("%e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e \n", (*Fz)*Lambda*Rbo*Rbo, Fpress*Lambda*Rbo*Rbo, Flaty*Lambda*Rbo*Rbo, (Mpressy*cos(phi) - Mpressx*sin(phi))*Lambda*Rbo*Rbo*Rbo, (Mlifty*cos(phi) - Mliftx*sin(phi))*Lambda*Rbo*Rbo*Rbo, Mlaty*Lambda*Rbo*Rbo*Rbo, Mfricy*Lambda*Rbo*Rbo*Rbo, (*My)*Lambda*Rbo*Rbo*Rbo);

	//	Dimensional Parameters
	FpressD = Fpress * (Lambda*pow(Rbo,2));
	MpressxD = (Mpressx*cos(phi) + Mpressy*sin(phi)) * (Lambda*pow(Rbo,3));
	MpressyD = (Mpressy*cos(phi) - Mpressx*sin(phi)) * (Lambda*pow(Rbo,3));	
	FliftD = Flift * (Lambda*pow(Rbo,2));
	MliftxD = (Mliftx*cos(phi) + Mlifty*sin(phi))  * (Lambda*pow(Rbo,3));
	MliftyD = (Mlifty*cos(phi) - Mliftx*sin(phi)) * (Lambda*pow(Rbo,3));
	FlatxD = Flatx*Lambda*Rbo*Rbo;
	FlatyD = Flaty*Lambda*Rbo*Rbo;
	MlatxD = Mlatx*Lambda*Rbo*Rbo*Rbo;
	MlatyD = Mlaty*Lambda*Rbo*Rbo*Rbo;
	FfricD = Ffric*Lambda*Rbo*Rbo;
	MfricxD = Mfricx*Lambda*Rbo*Rbo*Rbo;
	MfricyD = Mfricx*Lambda*Rbo*Rbo*Rbo;
	FzD = (*Fz)*Lambda*Rbo*Rbo;
	MxD = (*Mx)*Lambda*Rbo*Rbo*Rbo;
	MyD= (*My)*Lambda*Rbo*Rbo*Rbo;
	MzD = (*Mz)*Lambda*Rbo*Rbo*Rbo;
}
void FMpress(double phi, float **Rvb, double *Pc, float *Fpress, float *Mpressx, float *Mpressy)
{
	int i,j;
	float Abc, FpressD, MpressxD, MpressyD;	// Cylinder pressure closing area

	Abc = (Ap/pow(Rbo,2) - ((rB3*rB3-rB2*rB2)*(thetaK-2*rB/R)/2+Pi*rB*rB/pow(Rbo,2))); 
	
	//	Nondimensional pressing force
	*Fpress = -Fs;	
	for(i=1; i<=N; i++){
		*Fpress = (*Fpress) - Pc[i]*Abc;
	}
	//	Dimensional pressing force
	FpressD = (*Fpress) * (Lambda*pow(Rbo,2));

	//	Nondimensional pressing moment
	*Mpressx = 0.0;
	*Mpressy = 0.0;
	for(i=1; i<=N; i++){
		*Mpressx = (*Mpressx) - Pc[i]*Abc*R*cos(2*Pi/N*(i-1))/Rbo;
		*Mpressy = (*Mpressy) + Pc[i]*Abc*R*sin(2*Pi/N*(i-1))/Rbo;
	}

	//	Dimensional pressing moment
	MpressxD = (*Mpressx) * (Lambda*pow(Rbo,3));
	MpressyD = (*Mpressy) * (Lambda*pow(Rbo,3));
}
void FMlift(double phi, float **Rvb, float **Rn, float **Rs, float **Tvb, float **dT, float **Acv, float **P, float *Flift, float *Mliftx, float *Mlifty)
{
	int i,j;
	float FliftD, MliftxD, MliftyD;		

	//float temp=0.0;
	//	for(j=1; j<=NNBt; j+=NNBt/NNBt) 
	//		for(i=1; i<=NNBr; i++) temp = temp + P[i][j]; 
	//			//if(P[i][j] > 1.0) 
	//printf("%e  \n", temp);

	//	Nondimensional lifing force
	*Flift = 0.0;
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBr; i++){
			*Flift = (*Flift) + P[i][j]*Acv[i][j];			
			//printf("%d %d %e \n", i, j, *Flift);
		}
	}
	//printf("%e \n", *Flift);
	//	Dimensional lifing force
	FliftD = (*Flift) * (Lambda*pow(Rbo,2));
 
	*Mliftx = 0.0;
	*Mlifty = 0.0;
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBr; i++){
			*Mliftx = (*Mliftx) + P[i][j]*Acv[i][j]*Rvb[i][j]*cos(Tvb[i][j]);
			*Mlifty = (*Mlifty) - P[i][j]*Acv[i][j]*Rvb[i][j]*sin(Tvb[i][j]);
		}
	}

	//	Dimensional lifting moment
	MliftxD = (*Mliftx) * (Lambda*pow(Rbo,3));
	MliftyD = (*Mlifty) * (Lambda*pow(Rbo,3));

}

/* Piston friction */
void FMfriction(double phi, double *Pc, float *Ffric, float *Mfricx, float *Mfricy)
{
	int i,j;
	double *phiD, *Lpc, *Vp;

	phiD = dvector(1,N);
	Lpc = dvector(1,N);
	Vp = dvector(1,N);

	// Non-dimensional F-M
	*Ffric = 0.0;
	*Mfricx = 0.0;
	*Mfricy = 0.0;
	for(i=1;i<=N;i++){
		phiD[i] = phi + 2*Pi/N*(i-1);
		if(phiD[i] > 2*Pi) phiD[i] = phiD[i] - 2*Pi;
		Lpc[i] = Lc - Ld - R*(tan(beta)*(1+cos(phiD[i]))-tan(beta2)/cos(beta)*sin(phiD[i]));
		Vp[i] = -R*w*(tan(beta)*sin(phiD[i])+tan(beta2)/cos(beta)*cos(phiD[i]));
		*Ffric = *Ffric + (-(Pc[i]*Lambda-Ph)*hp/Lpc[i]/2.0+eta*Vp[i]/hp)*2*Pi*(Rp+hp)*Lpc[i] /Lambda/pow(Rbo,2);
		*Mfricx = *Mfricx + (-(Pc[i]*Lambda-Ph)*hp/Lpc[i]/2.0+eta*Vp[i]/hp)*2*Pi*(Rp+hp)*Lpc[i]*R*cos(phiD[i]) /Lambda/pow(Rbo,3);
		*Mfricy = *Mfricy - (-(Pc[i]*Lambda-Ph)*hp/Lpc[i]/2.0+eta*Vp[i]/hp)*2*Pi*(Rp+hp)*Lpc[i]*R*sin(phiD[i]) /Lambda/pow(Rbo,3);
	}

	free_dvector(phiD,1,N);	
	free_dvector(Lpc,1,N);
	free_dvector(Vp,1,N);
}
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
	for(i=1;i<=N;i++){
		phiD[i] = phi + 2*Pi/N*(i-1);
		while(phiD[i] > 2*Pi) phiD[i] = phiD[i] - 2*Pi;
		Lpc[i] = Lc - Ld - R*(tan(beta)*(1+cos(phiD[i]))-tan(beta2)/cos(beta)*sin(phiD[i]));
		Vp[i] = -R*w*(tan(beta)*sin(phiD[i])+tan(beta2)/cos(beta)*cos(phiD[i]));
		Ap[i] = -R*w*w*(tan(beta)*cos(phiD[i])-tan(beta2)/cos(beta)*sin(phiD[i]));
		Fpa = (Pc[i]*Lambda*Pi*Rp*Rp - Mp/w/w*Lambda*Rbo*(Rbo/hv0)*Ap[i] + ((Pc[i]*Lambda-Ph)*hp/Lpc[i]/2.0-eta*Vp[i]/hp)*2*Pi*(Rp+hp)*Lpc[i])/Lambda/pow(Rbo,2);
		Fwpa = Mp/w/w*Lambda*(Rbo/hv0)*Rbo*R*w*w/Lambda/pow(Rbo,2);
		Ffws = eta*R*w/hs*Pi*(Rs2*Rs2-Rs1*Rs1)/Lambda/pow(Rbo,2);
		z0y = (Fpa*tan(beta)+Fwpa*cos(phiD[i])+Ffws*sin(phiD[i]))/12/((Fpa*tan(beta)+Ffws*sin(phiD[i]))*(Lpj-Lpc[i]/2)+Fwpa*cos(phiD[i])*(Lpg-Lpc[i]/2))*pow(Lpc[i],2);
		z0x = (Fwpa*sin(phiD[i])-Ffws*cos(phiD[i]))/12/((Ffws*cos(phiD[i]))*(Lpj-Lpc[i]/2)-Fwpa*sin(phiD[i])*(Lpg-Lpc[i]/2))*pow(Lpc[i],2);
		*Flatx = *Flatx + (Fwpa*sin(phiD[i])-Ffws*cos(phiD[i]));					  
		*Flaty = *Flaty + (Fpa*tan(beta)+Fwpa*cos(phiD[i])+Ffws*sin(phiD[i]));
		*Mlaty = *Mlaty + (  ((z0x+Lpc[i]/2)/3+Loc)*(z0x+Lpc[i]/2)*((Fwpa*sin(phiD[i])-Ffws*cos(phiD[i]))/Lpc[i]+6*((Ffws*cos(phiD[i]))*(Lpj-Lpc[i]/2)-Fwpa*sin(phiD[i])*(Lpg-Lpc[i]/2))/pow(Lpc[i],2))/2
							- (Lpc[i]*5/6+z0x/3+Loc)*(z0x-Lpc[i]/2)*((Fwpa*sin(phiD[i])-Ffws*cos(phiD[i]))/Lpc[i]-6*((Ffws*cos(phiD[i]))*(Lpj-Lpc[i]/2)-Fwpa*sin(phiD[i])*(Lpg-Lpc[i]/2))/pow(Lpc[i],2))/2  )/Rbo;
		*Mlatx = *Mlatx + ( ((z0y+Lpc[i]/2)/3+Loc)*(z0y+Lpc[i]/2)*((Fpa*tan(beta)+Fwpa*cos(phiD[i])+Ffws*sin(phiD[i]))/Lpc[i]+6*((Fpa*tan(beta)+Ffws*sin(phiD[i]))*(Lpj-Lpc[i]/2)+Fwpa*cos(phiD[i])*(Lpg-Lpc[i]/2))/pow(Lpc[i],2))/2
							- (Lpc[i]*5/6+z0y/3+Loc)*(z0y-Lpc[i]/2)*((Fpa*tan(beta)+Fwpa*cos(phiD[i])+Ffws*sin(phiD[i]))/Lpc[i]-6*((Fpa*tan(beta)+Ffws*sin(phiD[i]))*(Lpj-Lpc[i]/2)+Fwpa*cos(phiD[i])*(Lpg-Lpc[i]/2))/pow(Lpc[i],2))/2 )/Rbo;		
		*Mlatz = *Mlatz + (R*cos(phi)*(Fpa*tan(beta)+Fwpa*cos(phiD[i])+Ffws*sin(phiD[i]))- R*sin(phi)*(Fwpa*sin(phiD[i])-Ffws*cos(phiD[i])))/Rbo ;

	}
	//printf("++++++++++++ %e %e %e %e \n", phi, *Flaty, *Mlatx, *Mlaty);

	free_dvector(phiD,1,N);	
	free_dvector(Lpc,1,N);
	free_dvector(Vp,1,N);
	free_dvector(Ap,1,N);
}