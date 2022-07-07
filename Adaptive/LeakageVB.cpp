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
extern double R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, Zevf, phif, phiout_c, phiout_vb, phiout_FM, phiout_LK, w, wp;
extern double V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
extern double integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, zB1, zB2, zB3, zB4, zB5, zB6, zB7, Fs, C, theta0, psi0, de, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
extern double Lambda, Loc, Lpj, Lpg, Lbj, Rbo, Lpc;
//extern float FxD, MyxD, MyyD, FyD, MxxD, MxyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, frD, MxD, MyD, MzD, leakageRinD, leakageRoutD;

/* Leakage in Cylinder Barrel & Valve plate */
void LeakageVB(float *leakage, float *PCleakage, float **Ppc, float **dN, float **dS, float **dT, float **dZ, float **Zn, double **Hpc, double **Hn, float **dPdZ_n, float **dPdT_c)
{
	//void Pgradient(float **Ppc, float **dPdZ_n, float **dPdT_c, float **dN, float **dS, float **dT);

	int i, j, k, k1, k2, NBKt, dNBKt;
	float **pointleakage_n, **pointleakage_c, leakageZin, leakageZout, leakageZinD, leakageZoutD;
	
	pointleakage_n = matrix(1,NNBz,1,NNBt);
	pointleakage_c = matrix(1,NNBz,1,NNBt);
	
	//Pgradient(Ppc, dPdZ_n, dPdT_c, dN, dS, dT);
	
	// leakage in each control volume	
	// (+) VALUE OF pointleakage MEANS FLOW FROM OUT TO IN

	// RADIAL DIRECTION LEAKAGE AT THE CONTROL VOLUME FACE n
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBz-1; i++){
			pointleakage_n[i][j] = (hp/Rp)*(-Hn[i][j]*Vp-pow(Hn[i][j],3)*dPdZ_n[i][j])*dT[i][j];
			pointleakage_c[i][j] = (hp/Rp)*(Hpc[i][j]*(w+wp)/w-pow(Hpc[i][j],3)*dPdT_c[i][j])*dZ[i][j];
		}
		pointleakage_n[NNBz][j] = pointleakage_n[NNBz-1][j];
		pointleakage_c[NNBz][j] = pointleakage_c[NNBz-1][j];
	}
	// Initialize
	for(i=1; i<=N;i++){		
		PCleakage[i] = 0.0;
	}

	//// Barrel Kidney(cylinder) leakage
	//NBKt = int(thetaK/dNBt+1.5);		// Barrel Kidney 각도방향 grid 수
	//dNBKt = int(NNBt/N+0.5);			// Barrel Kidney간격 각도방향 grid 수
	//for(i=1; i<N; i++){
	//	for(j=1+dNBKt*(i-1); j<=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j++){	
	//		PCleakage[i] = PCleakage[i] - pointleakage_n[NBz1+1][j] + pointleakage_n[NBz1+1+NBz2][j];
	//	}			
	//	for(j=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j<=dNBKt*i; j++){
	//		PCleakage[i+1] = PCleakage[i+1] - pointleakage_n[NBz1+1][j] + pointleakage_n[NBz1+1+NBz2][j];
	//	}		
	//	// Barrel rounding
	//	for(j=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(i-1); j++){
	//		k1 = 1+NBz1+NBz2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBz2/(zB3-zB2)/Rbo);
	//		k2 = 1+NBz1+NBz2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBz2/(zB3-zB2)/Rbo);
	//		PCleakage[i] = PCleakage[i] - pointleakage_n[k1][j] + pointleakage_n[k2][j] + pointleakage_c[k1][j] + pointleakage_c[k2][j];	
	//	}
	//	for(j=1+dNBKt*i-NBKt/2; j<=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j++){
	//		k1 = 1+NBz1+NBz2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBz2/(zB3-zB2)/Rbo);
	//		k2 = 1+NBz1+NBz2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBz2/(zB3-zB2)/Rbo);
	//		PCleakage[i+1] = PCleakage[i+1] - pointleakage_n[k1][j] + pointleakage_n[k2][j] - pointleakage_c[k1][j] - pointleakage_c[k2][j];
	//	}
	//}
	//for(j=1+dNBKt*(N-1); j<=1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt); j++){	
	//	PCleakage[N] = PCleakage[N] - pointleakage_n[NBz1+1][j] + pointleakage_n[NBz1+1+NBz2][j];
	//}			
	//for(j=1+dNBKt*N-NBKt/2+int(rB/R/dNBt); j<=dNBKt*N; j++){
	//	PCleakage[1] = PCleakage[1] - pointleakage_n[NBz1+1][j] + pointleakage_n[NBz1+1+NBz2][j];
	//}		
	//// Barrel rounding
	//for(j=1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(N-1); j++){
	//	k1 = 1+NBz1+NBz2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBz2/(zB3-zB2)/Rbo);
	//	k2 = 1+NBz1+NBz2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBz2/(zB3-zB2)/Rbo);
	//	PCleakage[N] = PCleakage[N] - pointleakage_n[k1][j] + pointleakage_n[k2][j] + pointleakage_c[k1][j] + pointleakage_c[k2][j];	
	//}
	//for(j=1+dNBKt*N-NBKt/2; j<=1+dNBKt*N-NBKt/2+int(rB/R/dNBt); j++){
	//	k1 = 1+NBz1+NBz2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBz2/(zB3-zB2)/Rbo);
	//	k2 = 1+NBz1+NBz2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBz2/(zB3-zB2)/Rbo);
	//	PCleakage[1] = PCleakage[1] - pointleakage_n[k1][j] + pointleakage_n[k2][j] - pointleakage_c[k1][j] - pointleakage_c[k2][j];
	//}



	// LEAKAGE[i] MEANS LEAKAGE AT THE CONTROL VOLUME FACE n OF Ith RADIUS
	// Initialize
	for(i=1; i<=NNBz;i++){
		leakage[i] = 0.0;
	}
	//	Leakage at the radius Zn (control volume face)
	for(j=1; j<=NNBt; j++){
		for(i=1; i<= NNBz; i++){
			leakage[i] = leakage[i] + pointleakage_n[i][j];
		}
	}


	leakageZin = leakage[1];        // DIMENSIONLESS LEAKAGE OF OUTTER
	leakageZout = leakage[NNBz];  // DIMENSIONLESS LEAKAGE OF INNER 

	leakageZinD = leakageZin*(pow(Rbo,3)*w);
	leakageZoutD = leakageZout*(pow(Rbo,3)*w);
	//for(i=1; i<=N; i++) PCleakage[i] = PCleakage[i]*(0.5*hv0*(Rbo*Rbo)*w);

	free_matrix(pointleakage_n,1,NNBz,1,NNBt);
	free_matrix(pointleakage_c,1,NNBz,1,NNBt);

}

/* pressure Gradient for Leakage Calculation */
void Pgradient(float **Ppc, float **dPdZ_n, float **dPdT_c, float **dN, float **dS, float **dT)
{
	int i, j, k, k1, k2, NBKt, dNBKt;

	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBz-1; i++){
			dPdZ_n[i][j] = (Ppc[i+1][j]-Ppc[i][j])/dN[i][j];
		}
		dPdZ_n[NNBz][j] = dPdZ_n[NNBz-1][j];
	}

	//	pressure GRADIENT AT EACH NODE POINT
	//	THEZe WAS NO DIFFEZeNCE BTW GRADIENT AT NODE AND AT CV FACE 
	//	j = 1
	for(i=1; i<=NNBz; i++){
		dPdT_c[i][1] = (Ppc[i][2]-Ppc[i][NNBt]) / (dT[i][1]*2.);
	}
	//	j = 2 ~ NNBt-1
	for(j=2; j<NNBt; j++){
		for(i=1; i<=NNBz; i++){
			dPdT_c[i][j] = (Ppc[i][j+1]-Ppc[i][j-1]) / (dT[i][j]*2.);
		}
	}
	//	j = NNBt
	for(i=1; i<=NNBz; i++){
		dPdT_c[i][NNBt] = (Ppc[i][1]-Ppc[i][NNBt-1]) / (dT[i][1]*2.);
	}

	//// Barrel Kidney(cylinder) leakage
	//NBKt = int(thetaK/dNBt+1.5);		// Barrel Kidney 각도방향 grid 수
	//dNBKt = int(NNBt/N+0.5);			// Barrel Kidney간격 각도방향 grid 수
	//for(i=1; i<=N; i++){
	//	for(j=1+dNBKt*(i-1); j<=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j++){		
	//		dPdZ_n[NBz1+1][j] = (Ppc[NBz1+1][j]-Ppc[NBz1][j])/dS[NBz1+1][j];
	//		dPdZ_n[NBz1+1+NBz2][j] = (Ppc[NBz1+1+NBz2+1][j]-Ppc[NBz1+1+NBz2][j])/dN[NBz1+1+NBz2][j];
	//	}			
	//	for(j=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j<=dNBKt*i; j++){
	//		dPdZ_n[NBz1+1][j] = (Ppc[NBz1+1][j]-Ppc[NBz1][j])/dS[NBz1+1][j];
	//		dPdZ_n[NBz1+1+NBz2][j] = (Ppc[NBz1+1+NBz2+1][j]-Ppc[NBz1+1+NBz2][j])/dN[NBz1+1+NBz2][j];
	//	}		
	//	// Barrel rounding
	//	for(j=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(i-1); j++){
	//		k1 = 1+NBz1+NBz2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBz2/(zB3-zB2)/Rbo);
	//		k2 = 1+NBz1+NBz2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBz2/(zB3-zB2)/Rbo);
	//		dPdZ_n[k1][j] = (Ppc[k1][j]-Ppc[k1-1][j])/dS[k1][j];
	//		dPdZ_n[k2][j] = (Ppc[k2+1][j]-Ppc[k2][j])/dN[k2][j];
	//		dPdT_c[k1][j] = (Ppc[k1][j+1]-Ppc[k1][j])/(dT[k1][j]*2.);
	//		dPdT_c[k2][j] = (Ppc[k2][j+1]-Ppc[k2][j])/(dT[k2][j]*2.);			
	//	}
	//	for(j=1+dNBKt*i-NBKt/2; j<=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j++){
	//		k1 = 1+NBz1+NBz2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBz2/(zB3-zB2)/Rbo);
	//		k2 = 1+NBz1+NBz2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBz2/(zB3-zB2)/Rbo);
	//		dPdZ_n[k1][j] = (Ppc[k1][j]-Ppc[k1-1][j])/dS[k1][j];
	//		dPdZ_n[k2][j] = (Ppc[k2+1][j]-Ppc[k2][j])/dN[k2][j];
	//		dPdT_c[k1][j] = (Ppc[k1][j]-Ppc[k1][j-1])/(dT[k1][j]*2.);
	//		dPdT_c[k2][j] = (Ppc[k2][j]-Ppc[k2][j-1])/(dT[k2][j]*2.);	
	//	}
	//}

}