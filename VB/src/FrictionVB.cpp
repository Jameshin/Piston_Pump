#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

extern int Nphi, NBr1, NBr2, NBr3, NBr4, NBr5, NBr, NNBr, NNBt, Nch;
extern int NVt1, NVt2, NVt3, NVt4, NVt5, NVt6, NVr1, NVr2, NVr3, NVr4;
extern float dNBt, thetach;
extern float R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, revf, phif, phiout_c, hm, phiout_cy, phiout_dy, phiout_vb, phiout_FM, phiout_LK, w;
extern float V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
extern float integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, rB1, rB2, rB3, rB4, rB5, rB6, Fs, hv0, theta0, psi0, dhv0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
extern float Lambda, Loc, Lpj, Lpg, Lbj, Ls, Rbo;
extern float FpressD, MpressxD, MpressyD, FliftD, MliftxD, MliftyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, FzD, MxD, MyD, MzD, leakageRinD, leakageRoutD, frictorqueD, friclossD, frictorqueD2;

void FrictionVB(double phi, float *Tb, float *Tv, float **Rvb, double **Hvb, float **Acv, float **dPdT_c, float **Pvb)
{
	int i, j, k;
	float **Tb_cv, **Tv_cv;
	float **Acvf;

	Tb_cv = matrix(1,NNBr,1,NNBt);
	Tv_cv = matrix(1,NNBr,1,NNBt);
	Acvf = matrix(1,NNBr,1,NNBt);	

	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBr; i++){
			Acvf[i][j] = Acv[i][j];
		}
	}

	/* Valve Port */
	for(j=1+NVt1-int(phi/dNBt); j<=1+NVt1+int(rV1/R/dNBt)-int(phi/dNBt); j++){
		for(i=(NVr1+NVr2)/2-int(sqrt(rV1*rV1-pow(R*dNBt*(j-(1+NVt1-int(phi/dNBt)))-rV1,2))*NBr2/(rB3-rB2)/Rbo); 
			i<=(NVr1+NVr2)/2+int(sqrt(rV1*rV1-pow(R*dNBt*(j-(1+NVt1-int(phi/dNBt)))-rV1,2))*NBr2/(rB3-rB2)/Rbo); i++){
			k = j;
			if(j<1) k = j + NNBt;
			Acvf[i][k] =0.0;   // discharge valve port
		}
	}
	for(j=1+NVt1+int(rV1/R/dNBt)-int(phi/dNBt); j<=1+NVt1+NVt2-int(rV1/R/dNBt)-int(phi/dNBt); j++){
		for(i=NVr1; i<=NVr2; i++){
			k = j;
			if(j<1) k = j + NNBt;
			Acvf[i][k] =0.0;   // discharge valve port
		}
	}
	for(j=1+NVt1+NVt2-int(rV1/R/dNBt)-int(phi/dNBt); j<=1+NVt1+NVt2-int(phi/dNBt); j++){
		for(i=(NVr1+NVr2)/2-int(sqrt(rV1*rV1-pow(R*dNBt*(j-(1+NVt1+NVt2-int(rV1/R/dNBt)-int(phi/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); 
			i<=(NVr1+NVr2)/2+int(sqrt(rV1*rV1-pow(R*dNBt*(j-(1+NVt1+NVt2-int(rV1/R/dNBt)-int(phi/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); i++){
			k = j;
			if(j<1) k = j + NNBt;
			Acvf[i][k] =0.0;   // discharge valve port
		}
	}
	for(j=1+NVt1+NVt2+NVt3+NVt4-int(phi/dNBt); j<=1+NVt1+NVt2+NVt3+NVt4+int(rV2/R/dNBt)-int(phi/dNBt); j++){
		for(i=(NVr3+NVr4)/2-int(sqrt(rV2*rV2-pow(R*dNBt*(j-(1+NVt1+NVt2+NVt3+NVt4-int(phi/dNBt)))-rV2,2))*NBr2/(rB3-rB2)/Rbo); 
			i<=(NVr3+NVr4)/2+int(sqrt(rV2*rV2-pow(R*dNBt*(j-(1+NVt1+NVt2+NVt3+NVt4-int(phi/dNBt)))-rV2,2))*NBr2/(rB3-rB2)/Rbo); i++){
			k = j;
			if(j<1) k = j + NNBt;
			Acvf[i][k] =0.0;	// suction valve port
		}
	}
		// Valve rounding
	for(j=1+NVt1+NVt2+NVt3+NVt4+int(rV2/R/dNBt)-int(phi/dNBt); j<=1+NVt1+NVt2+NVt3+NVt4+NVt5-int(rV2/R/dNBt)-int(phi/dNBt); j++){
		for(i=NVr3; i<=NVr4; i++){
			k = j;
			if(j<1) k = j + NNBt;
			Acvf[i][k] =0.0;	// suction valve port
		}
	}
	for(j=1+NVt1+NVt2+NVt3+NVt4+NVt5-int(rV2/R/dNBt)-int(phi/dNBt); j<=1+NVt1+NVt2+NVt3+NVt4+NVt5-int(phi/dNBt); j++){
		for(i=(NVr3+NVr4)/2-int(sqrt(rV2*rV2-pow(R*dNBt*(j-(1+NVt1+NVt2+NVt3+NVt4+NVt5-int(rV2/R/dNBt)-int(phi/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); 
			i<=(NVr3+NVr4)/2+int(sqrt(rV2*rV2-pow(R*dNBt*(j-(1+NVt1+NVt2+NVt3+NVt4+NVt5-int(rV2/R/dNBt)-int(phi/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); i++){
			k = j;
			if(j<1) k = j + NNBt;
			Acvf[i][k] =0.0;	// suction valve port
		}
	}

	
	// Friction torque in a control volume
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBr; i++){			
			Tb_cv[i][j] = (Rvb[i][j]/Hvb[i][j]-3.0*dPdT_c[i][j])*Acvf[i][j]*Rvb[i][j];
			Tv_cv[i][j] = (Rvb[i][j]/Hvb[i][j]+3.0*dPdT_c[i][j])*Acvf[i][j]*Rvb[i][j];
			if(Pvb[i][j] <= Pca/Lambda){
				Tb_cv[i][j] = 0.0;
				Tv_cv[i][j] = 0.0;
			}
		}
	}


	*Tb = 0.0;  
	*Tv = 0.0;  
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBr; i++){
			*Tb = *Tb + Tb_cv[i][j];
			*Tv = *Tv + Tv_cv[i][j];
		}
	}
	

	frictorqueD = (*Tb)*(eta*w*pow(Rbo,4)/hv0);
	frictorqueD2 = (*Tv)*(eta*w*pow(Rbo,4)/hv0);
	friclossD = frictorqueD*w;

	free_matrix(Tb_cv,1,NNBr,1,NNBt);
	free_matrix(Tv_cv,1,NNBr,1,NNBt);
	free_matrix(Acvf,1,NNBr,1,NNBt);

}