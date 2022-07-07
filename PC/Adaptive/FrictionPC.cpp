#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

extern float Vp, **dPdZ_n, **dPdT_c;
extern void Pgradient(float **Ppc, float **dPdZ_n, float **dPdT_c, float **dN, float **dS, float **dT);
extern void FrictionPC(double phi, float *Ffric, float *Mfricx, float *Mfricy, float **Zpc, float **Tpc, double **Hpc, float **Acv, float **dPdZ_n, float **Ppc);
extern int Nphi, NBz1, NBz2, NBz3, NBz4, NBz5, NBz6, NBz, NNBz, NNBt, Nch;
extern int NVt1, NVt2, NVt3, NVt4, NVt5, NVt6, NVr1, NVr2, NVr3, NVr4;
extern double dNBt, thetach;
extern double R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, Zevf, phif, phiout_c, hm, phiout_cy, phiout_dy, phiout_vb, phiout_FM, phiout_LK, w;
extern double V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
extern double integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, zB1, zB2, zB3, zB4, zB5, zB6, zB7, Fs, C, theta0, psi0, de, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
extern double Lambda, Loc, Lpj, Lpg, Lbj, Ls, Rbo, Lpc;
extern float FxhD, FyhD, MxhD, MyhD, FxcD, FycD, MxcD, MycD, FxD, FyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, FzD, MxD, MyD, MzD, leakageRinD, leakageRoutD, frictionD, friclossD, Hmin;

void FrictionPC(double phi, float *Ffp, float *Mfxp, float *Mfyp, float **Zpc, float **Tpc, double **Hpc, float **Acv, float **dPdZ_n, float **Ppc)
{
	int i, j, k;
	float **Tb_cv, **Ffp_cv, **Mfxp_cv, **Mfyp_cv;
	float **Acvf;

	//Tb_cv = matrix(1,NNBz,1,NNBt);
	Ffp_cv = matrix(1,NNBz,1,NNBt);
	Mfxp_cv = matrix(1,NNBz,1,NNBt);
	Mfyp_cv = matrix(1,NNBz,1,NNBt);
	Acvf = matrix(1,NNBz,1,NNBt);	

	//for(j=1; j<=NNBt; j++){
	//	for(i=1; i<=NNBz; i++){
	//		Acvf[i][j] = Acv[i][j];
	//	}
	//}

	// Friction torque in a control volume
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBz; i++){			
			//Tb_cv[i][j] = (Zpc[i][j]/Hpc[i][j]-3.0*dPdT_c[i][j])*Acvf[i][j]*Zpc[i][j];
			Ffp_cv[i][j] = (hp/Rp)/6*(Vp/Hpc[i][j]+3.0*Hpc[i][j]*dPdZ_n[i][j])*Acv[i][j];
			Mfxp_cv[i][j] = -(hp/Rp)/6*(Vp/Hpc[i][j]+3.0*Hpc[i][j]*dPdZ_n[i][j])*Acv[i][j]*cos(Tpc[i][j]);
			Mfyp_cv[i][j] = (hp/Rp)/6*(Vp/Hpc[i][j]+3.0*Hpc[i][j]*dPdZ_n[i][j])*Acv[i][j]*sin(Tpc[i][j]);
			if(Ppc[i][j] <= Pca/Lambda){
				Ffp_cv[i][j] = 0.0;
				Mfxp_cv[i][j] = 0.0;
				Mfyp_cv[i][j] = 0.0;
			}
		}
	}

	*Ffp = 0.0;
	//*Tb = 0.0;
	*Mfxp = 0.0;
	*Mfyp = 0.0;  
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBz; i++){
			*Ffp = *Ffp + Ffp_cv[i][j];
			//*Tb = *Tb + Tb_cv[i][j];
			*Mfxp = *Mfxp + Mfxp_cv[i][j];
			*Mfyp = *Mfyp + Mfyp_cv[i][j];
		}
	}	

	frictionD = (*Ffp)*(Lambda*pow(Rp,2));
	//frictionD2 = (*Tv)*(eta*w*pow(Rbo,4)/hv0);
	friclossD = frictionD*w;

	//free_matrix(Tb_cv,1,NNBz,1,NNBt);
	free_matrix(Ffp_cv,1,NNBz,1,NNBt);
	free_matrix(Mfxp_cv,1,NNBz,1,NNBt);
	free_matrix(Mfyp_cv,1,NNBz,1,NNBt);
	free_matrix(Acvf,1,NNBz,1,NNBt);

}