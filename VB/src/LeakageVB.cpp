#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

extern int Nphi, NBr1, NBr2, NBr3, NBr4, NBr5, NBr, NNBr, NNBt, Nch;
extern int NVt1, NVt2, NVt3, NVt4, NVt5, NVt6, NVr1, NVr2, NVr3, NVr4;
extern float dNBt, thetach;
extern float R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, revf, phif, phiout_c, phiout_vb, phiout_FM, phiout_LK, w;
extern float V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
extern float integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, rB1, rB2, rB3, rB4, rB5, rB6, Fs, hv0, theta0, psi0, dhv0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
extern float Lambda, Loc, Lpj, Lpg, Lbj, Rbo;
extern float FpressD, MpressxD, MpressyD, FliftD, MliftxD, MliftyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, FzD, MxD, MyD, MzD, leakageRinD, leakageRoutD;

/* Leakage in Cylinder Barrel & Valve plate */
void LeakageVB(float *leakage, float *Baleakage, float **Pvb, float **dN, float **dS, float **dT, float **dR, float **Rn, double **Hvb, double **Hn, float **dPdR_n, float **dPdT_c)
{
	//void Pgradient(float **Pvb, float **dPdR_n, float **dPdT_c, float **dN, float **dS, float **dT);

	int i, j, k, k1, k2, NBKt, dNBKt;
	float /***dPdR_n, **dPdT_c,*/ **pointleakage_n, **pointleakage_c, leakageRin, leakageRout;
	
	//dPdR_n = matrix(1,NNBr,1,NNBt);
	//dPdT_c = matrix(1,NNBr,1,NNBt);
	pointleakage_n = matrix(1,NNBr,1,NNBt);
	pointleakage_c = matrix(1,NNBr,1,NNBt);
	
	//Pgradient(Pvb, dPdR_n, dPdT_c, dN, dS, dT);
	
	// leakage in each control volume	
	// (+) VALUE OF pointleakage MEANS FLOW FROM OUT TO IN

	// RADIAL DIRECTION LEAKAGE AT THE CONTROL VOLUME FACE n
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBr-1; i++){
			pointleakage_n[i][j] = -pow(Hn[i][j],3)*dPdR_n[i][j]*Rn[i][j]*dT[i][j];
			pointleakage_c[i][j] = -pow(Hvb[i][j],3)*dPdT_c[i][j]*dR[i][j];
		}
	}
	// Initialize
	for(i=1; i<=N;i++){		
		Baleakage[i] = 0.0;
	}

	// Barrel Kidney(cylinder) leakage
	NBKt = int(thetaK/dNBt+1.5);		// Barrel Kidney 각도방향 grid 수
	dNBKt = int(NNBt/N+0.5);			// Barrel Kidney간격 각도방향 grid 수
	for(i=1; i<N; i++){
		for(j=1+dNBKt*(i-1); j<=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j++){	
			Baleakage[i] = Baleakage[i] - pointleakage_n[NBr1+1][j] + pointleakage_n[NBr1+1+NBr2][j];
		}			
		for(j=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j<=dNBKt*i; j++){
			Baleakage[i+1] = Baleakage[i+1] - pointleakage_n[NBr1+1][j] + pointleakage_n[NBr1+1+NBr2][j];
		}		
		// Barrel rounding
		for(j=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(i-1); j++){
			k1 = 1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo);
			k2 = 1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo);
			Baleakage[i] = Baleakage[i] - pointleakage_n[k1][j] + pointleakage_n[k2][j] + pointleakage_c[k1][j] + pointleakage_c[k2][j];	
		}
		for(j=1+dNBKt*i-NBKt/2; j<=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j++){
			k1 = 1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo);
			k2 = 1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo);
			Baleakage[i+1] = Baleakage[i+1] - pointleakage_n[k1][j] + pointleakage_n[k2][j] - pointleakage_c[k1][j] - pointleakage_c[k2][j];
		}
	}
	for(j=1+dNBKt*(N-1); j<=1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt); j++){	
		Baleakage[N] = Baleakage[N] - pointleakage_n[NBr1+1][j] + pointleakage_n[NBr1+1+NBr2][j];
	}			
	for(j=1+dNBKt*N-NBKt/2+int(rB/R/dNBt); j<=dNBKt*N; j++){
		Baleakage[1] = Baleakage[1] - pointleakage_n[NBr1+1][j] + pointleakage_n[NBr1+1+NBr2][j];
	}		
	// Barrel rounding
	for(j=1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(N-1); j++){
		k1 = 1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo);
		k2 = 1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo);
		Baleakage[N] = Baleakage[N] - pointleakage_n[k1][j] + pointleakage_n[k2][j] + pointleakage_c[k1][j] + pointleakage_c[k2][j];	
	}
	for(j=1+dNBKt*N-NBKt/2; j<=1+dNBKt*N-NBKt/2+int(rB/R/dNBt); j++){
		k1 = 1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo);
		k2 = 1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo);
		Baleakage[1] = Baleakage[1] - pointleakage_n[k1][j] + pointleakage_n[k2][j] - pointleakage_c[k1][j] - pointleakage_c[k2][j];
	}



	// LEAKAGE[i] MEANS LEAKAGE AT THE CONTROL VOLUME FACE n OF Ith RADIUS
	// Initialize
	for(i=1; i<=NNBr-1;i++){
		leakage[i] = 0.0;
	}
	//	Leakage at the radius Rn (control volume face)
	for(j=1; j<=NNBt; j++){
		for(i=1; i<= NNBr-1; i++){
			leakage[i] = leakage[i] + pointleakage_n[i][j];
		}
	}


	leakageRin = leakage[1];        // DIMENSIONLESS LEAKAGE OF INNER 
	leakageRout = leakage[NNBr-NBr4-NBr5];  // DIMENSIONLESS LEAKAGE OF OUTTER

	leakageRinD = leakageRin*(0.5*hv0*(Rbo*Rbo)*w);
	leakageRoutD = leakageRout*(0.5*hv0*(Rbo*Rbo)*w);
	for(i=1; i<=N; i++) Baleakage[i] = Baleakage[i]*(0.5*hv0*(Rbo*Rbo)*w);

	//free_matrix(dPdR_n,1,NNBr,1,NNBt);
	//free_matrix(dPdT_c,1,NNBr,1,NNBt);
	free_matrix(pointleakage_n,1,NNBr,1,NNBt);
	free_matrix(pointleakage_c,1,NNBr,1,NNBt);

}

/* Pressure Gradient for Leakage Calculation */
void Pgradient(float **Pvb, float **dPdR_n, float **dPdT_c, float **dN, float **dS, float **dT)
{
	int i, j, k, k1, k2, NBKt, dNBKt;

	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBr-1; i++){
			dPdR_n[i][j] = (Pvb[i+1][j]-Pvb[i][j])/dN[i][j];
		}
	}

	//	PRESSURE GRADIENT AT EACH NODE POINT
	//	THERE WAS NO DIFFERENCE BTW GRADIENT AT NODE AND AT CV FACE 
	//	j = 1
	for(i=1; i<=NNBr; i++){
		dPdT_c[i][1] = (Pvb[i][2]-Pvb[i][NNBt]) / (dT[i][1]*2.);
	}
	//	j = 2 ~ NNBt-1
	for(j=2; j<NNBt; j++){
		for(i=1; i<=NNBr; i++){
			dPdT_c[i][j] = (Pvb[i][j+1]-Pvb[i][j-1]) / (dT[i][j]*2.);
		}
	}
	//	j = NNBt
	for(i=1; i<=NNBr; i++){
		dPdT_c[i][NNBt] = (Pvb[i][1]-Pvb[i][NNBt-1]) / (dT[i][1]*2.);
	}

	// Barrel Kidney(cylinder) leakage
	NBKt = int(thetaK/dNBt+1.5);		// Barrel Kidney 각도방향 grid 수
	dNBKt = int(NNBt/N+0.5);			// Barrel Kidney간격 각도방향 grid 수
	for(i=1; i<=N; i++){
		for(j=1+dNBKt*(i-1); j<=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j++){		
			dPdR_n[NBr1+1][j] = (Pvb[NBr1+1][j]-Pvb[NBr1][j])/dS[NBr1+1][j];
			dPdR_n[NBr1+1+NBr2][j] = (Pvb[NBr1+1+NBr2+1][j]-Pvb[NBr1+1+NBr2][j])/dN[NBr1+1+NBr2][j];
		}			
		for(j=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j<=dNBKt*i; j++){
			dPdR_n[NBr1+1][j] = (Pvb[NBr1+1][j]-Pvb[NBr1][j])/dS[NBr1+1][j];
			dPdR_n[NBr1+1+NBr2][j] = (Pvb[NBr1+1+NBr2+1][j]-Pvb[NBr1+1+NBr2][j])/dN[NBr1+1+NBr2][j];
		}		
		// Barrel rounding
		for(j=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(i-1); j++){
			k1 = 1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo);
			k2 = 1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo);
			dPdR_n[k1][j] = (Pvb[k1][j]-Pvb[k1-1][j])/dS[k1][j];
			dPdR_n[k2][j] = (Pvb[k2+1][j]-Pvb[k2][j])/dN[k2][j];
			dPdT_c[k1][j] = (Pvb[k1][j+1]-Pvb[k1][j])/(dT[k1][j]*2.);
			dPdT_c[k2][j] = (Pvb[k2][j+1]-Pvb[k2][j])/(dT[k2][j]*2.);			
		}
		for(j=1+dNBKt*i-NBKt/2; j<=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j++){
			k1 = 1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo);
			k2 = 1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo);
			dPdR_n[k1][j] = (Pvb[k1][j]-Pvb[k1-1][j])/dS[k1][j];
			dPdR_n[k2][j] = (Pvb[k2+1][j]-Pvb[k2][j])/dN[k2][j];
			dPdT_c[k1][j] = (Pvb[k1][j]-Pvb[k1][j-1])/(dT[k1][j]*2.);
			dPdT_c[k2][j] = (Pvb[k2][j]-Pvb[k2][j-1])/(dT[k2][j]*2.);	
		}
	}

}