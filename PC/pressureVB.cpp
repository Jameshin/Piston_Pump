#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

extern float Vp, **dPdZ_n, **dPdT_c;
extern void Pgradient(float **Ppc, float **dPdZ_n, float **dPdT_c, float **dN, float **dS, float **dT);
extern void FrictionPC(double phi, float *Ffric, float *Mfricx, float *Mfricy, float **Zpc, float **Tpc, double **Hpc, float **Acv, float **dPdZ_n, float **Ppc);
extern int Nphi, NBz0, NBz1, NBz2, NBz3, NBz4, NBz5, NBz6, NBz, NNBz, NNBt, Nch, Imin, Jmin;
extern int NVt1, NVt2, NVt3, NVt4, NVt5, NVt6, NVr1, NVr2, NVr3, NVr4;
extern double dNBt, thetach;
extern double R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, revf, phif, phiout_c, phiout_vb, phiout_FM, phiout_LK, w, wp;
extern double V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
extern double integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, zB1, zB2, zB3, zB4, zB5, zB6, zB7, Fs, e0, theta0, psi0, de0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
extern double Lambda, Loc, Lpj, Lpg, Lbj, Rbo, Lpc;
extern float FxD, MyxD, MyyD, FyD, MxxD, MxyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, frD, MxD, MyD, MzD, leakageRinD, leakageRoutD, Hmin;
double gamH, zetaH;

/* Update Grid & Coordinate */
void Upgrid(float **Zpc, float **Zn, float **Zs, float **Zw, float **Ze, float **T, float **Tn, float **Ts, float **Tw, float **Te, float **dZ, float **dT, float **dN, float **dS, float **dW, float **dE, float **Acv, float *Acvt)
{
	int i,j;

	zB7 = Lpc;
	NBz6 = int(2*(zB7-zB6)/(zB6-zB5)*NBz5);
	NBz = NBz0 + NBz1 + NBz2 + NBz3 + NBz4 + NBz5 + NBz6;
	NNBz = NBz + 1;

	if(NBz6 != 0){
		// dr5
		for(j=1; j<=NNBt; j++){
			dN[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j] = (zB7-zB6)/NBz6;
			if(dN[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j] == 0.0) dN[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j] = TINY;
			dZ[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j] = (dN[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j]+dS[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5][j])/2.0;
		}
		// dr6
		for(j=1; j<=NNBt; j++){
			for(i=1; i<NBz6; i++){
				dN[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j] = (zB7-zB6)/NBz6;
				dS[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j] = (zB7-zB6)/NBz6;
				dW[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j] = 2*Pi/NNBt;
				dE[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j] = 2*Pi/NNBt;
				dZ[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j] = (dN[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j]+dS[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j])/2.0;
				dT[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j] = (dW[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j]+dE[NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5+i][j])/2.0;
			}
			dN[NNBz][j] = TINY;
			dS[NNBz][j] = (zB7-zB6)/NBz6;
			dW[NNBz][j] = 2*Pi/NNBt;
			dE[NNBz][j] = 2*Pi/NNBt;
			dZ[NNBz][j] = (dN[NNBz][j]+dS[NNBz][j])/2.0;
			dT[NNBz][j] = (dW[NNBz][j]+dE[NNBz][j])/2.0;
		}
	}

	// Coordinate of each node
	//for(j=1; j<=NNBt; j++){
	//	Zpc[1][j] = -Lpg/Rp;
	//	Zw[1][j] = Zpc[1][j];            
	//	Ze[1][j] = Zpc[1][j];           
	//	Zn[1][j] = Zpc[1][j] + dN[1][j]/2.0; 
	//	Zs[1][j] = Zpc[1][j] - dS[1][j]/2.0; 
	//}
	for(j=1; j<=NNBt; j++){
		for(i=NBz0+NBz1+1+NBz2+NBz3+NBz4+NBz5; i<=NNBz; i++){
			Zpc[i][j] = Zpc[i-1][j] + dS[i][j];	 // Radius of node[i][j]
			Zw[i][j] = Zpc[i][j];            
			Ze[i][j] = Zpc[i][j];           
			Zn[i][j] = Zpc[i][j] + dN[i][j]/2.0; 
			Zs[i][j] = Zpc[i][j] - dS[i][j]/2.0; 
		}
	}
	// Coordinate of each node
	for(i=NBz0+NBz1+2+NBz2+NBz3+NBz4+NBz5; i<=NNBz; i++){
		T[i][1] = 0.0;
		Tw[i][1] = T[i][j] - dW[i][j]/2.0;
		Te[i][1] = T[i][j] + dE[i][j]/2.0;
		Tn[i][1] = T[i][j]; 
		Ts[i][1] = T[i][j]; 
	}
	for(j=2; j<=NNBt; j++){
		for(i=NBz0+NBz1+2+NBz2+NBz3+NBz4+NBz5; i<=NNBz; i++){
			T[i][j] = T[i][j-1] + dW[i][j];	// Angle of node[i][j]
			Tw[i][j] = T[i][j] - dW[i][j]/2.0;
			Te[i][j] = T[i][j] + dE[i][j]/2.0;
			Tn[i][j] = T[i][j]; 
			Ts[i][j] = T[i][j]; 
		}
	}

	//	area of each control volume
	for(j=1; j<=NNBt; j++){
		for(i=1+NBz0+NBz1+NBz2+NBz3+NBz4+NBz5; i<=NNBz; i++){
			Acv[i][j] = (Zn[i][j]-Zs[i][j])*dT[i][j]; 
		}
		//Acv[1][j] = 0.5*(Zn[i][j]*Zn[i][j]-Zs[i][j]*Zs[i][j])*dT[i][j]*0.5;   // half-CV
		//Acv[NNBz][j] = 0.5*(Zn[i][j]*Zn[i][j]-Zs[i][j]*Zs[i][j])*dT[i][j]*0.5;  // half-CV
	}	

	//	Total Nondimensional face area
	*Acvt = 0.0;		
	for(j=1; j<=NNBt; j++){
		for(i=1+NBz0+NBz1+NBz2+NBz3+NBz4+NBz5; i<=NNBz; i++){
			*Acvt = *Acvt + Acv[i][j];
		}
	}
	//printf("%e \n", Zpc[NNBz][1]);
}

/* Film thickness */
void filmthickness(double phi, double *ym, double **H, double **Hn, double **Hs, double **Hw, double **He, double **H_gro, float **Z, float **T, double e, double kai, double gamma, double zeta)
{
	int i,j;
	//gamH = 0.005*Pi/180*Rp/hp;
	//zetaH = 40.0*Pi/180;

	// Solid contact
	//if(ym[1]-gamma <= Hpcmin/hv0){
	//	printf("Minus Film Thickness < Hmin at Hmin = %f \n", (ym[1]-gamma)*hv0);	
	//	//gamma = ym[1] - Hpcmin/hv0;
	//	for(j=1; j<=NNBt; j++){
	//		for(i=1; i<=NNBz; i++){
	//			H[i][j] = ym[1] - Zpc[i][j]*gamma*cos(T[i][j]-(Pi/2-(phi+zeta)))  + H_gro[i][j];
	//		
	//		}
	//	}		
	//}
	//else{
		//	FILM THICKNESS
		Hmin = H[1][1];
		for(j=1; j<=NNBt; j++){
			for(i=1; i<=NNBz; i++){
				H[i][j] = 1 - e*sin(T[i][j]+kai) - Z[i][j]*gamma*sin(T[i][j]+zeta) + H_gro[i][j];	
				//printf("%e %e %e \n", gamma, zeta, H[i][j]);
				if(Hmin > H[i][j]){
					Hmin = H[i][j];
					Imin = i;
					Jmin = j;
				}
			}
		}
		
	//}

	for(j=2; j<NNBt; j++){
		for(i=2; i<NNBz; i++){
			Hn[i][j] = (H[i][j]+H[i+1][j])/2.0;
			Hs[i][j] = (H[i][j]+H[i-1][j])/2.0;
			Hw[i][j] = (H[i][j]+H[i][j-1])/2.0;
			He[i][j] = (H[i][j]+H[i][j+1])/2.0;
		}
	}
	// j=1
	for(i=2; i<NNBz; i++){
		Hn[i][1] = (H[i][1]+H[i+1][1])/2.0;
		Hs[i][1] = (H[i][1]+H[i-1][1])/2.0;
		Hw[i][1] = (H[i][1]+H[i][NNBt])/2.0;
		He[i][1] = (H[i][1]+H[i][2])/2.0;
	}
	// j=NNBt
	for(i=2; i<NNBz; i++){
		Hn[i][NNBt] = (H[i][NNBt]+H[i+1][NNBt])/2.0;
		Hs[i][NNBt] = (H[i][NNBt]+H[i-1][NNBt])/2.0;
		Hw[i][NNBt] = (H[i][NNBt]+H[i][NNBt-1])/2.0;
		He[i][NNBt] = (H[i][NNBt]+H[i][1])/2.0;
	}
	// i=1
	for(j=2; j<NNBt; j++){
		Hn[1][j] = (H[1][j]+H[2][j])/2.0;
		Hw[1][j] = (H[1][j]+H[1][j-1])/2.0;
		He[1][j] = (H[1][j]+H[1][j+1])/2.0;	
	}
	Hn[1][1] = (H[1][1]+H[2][1])/2.0;
	Hw[1][1] = (H[1][1]+H[1][NNBt])/2.0;
	He[1][1] = (H[1][1]+H[1][2])/2.0;
	Hn[1][NNBt] = (H[1][NNBt]+H[2][NNBt])/2.0;
	Hw[1][NNBt] = (H[1][NNBt]+H[1][NNBt-1])/2.0;
	He[1][NNBt] = (H[1][NNBt]+H[1][1])/2.0;
	// i=NNBz
	for(j=2; j<NNBt; j++){
		Hs[NNBz][j] = (H[NNBz][j]+H[NNBz-1][j])/2.0;
		Hw[NNBz][j] = (H[NNBz][j]+H[NNBz][j-1])/2.0;
		He[NNBz][j] = (H[NNBz][j]+H[NNBz][j+1])/2.0;	
	}
	Hs[NNBz][1] = (H[NNBz][1]+H[NNBz-1][1])/2.0;
	Hw[NNBz][1] = (H[NNBz][1]+H[NNBz][NNBt])/2.0;
	He[NNBz][1] = (H[NNBz][1]+H[NNBz][2])/2.0;
	Hs[NNBz][NNBt] = (H[NNBz][NNBt]+H[NNBz-1][NNBt])/2.0;
	Hw[NNBz][NNBt] = (H[NNBz][NNBt]+H[NNBz][NNBt-1])/2.0;
	He[NNBz][NNBt] = (H[NNBz][NNBt]+H[NNBz][1])/2.0;
}
/* FVM coefficient */
void descretize(double phi, double *ym, double e, double kai, double gamma, double zeta, double de, double dkai, double dgamma, double dzeta, float **Zpc, float **dT, float **dZ, double **Hn, double **Hs, double **Hw, double **He, 
				float **Zn, float **Zs, float **Zw, float **Ze, float **Tn, float **Ts, float **Tw, float **Te, float **dN, float **dS, float **dW, float **dE, 
				float **aN, float **aS, float **aW, float **aE, float **aC, float **B)
{
	int i,j;
	float **Bwe1, **Bwe2, **Bsq1, **Bsq2;
	Bwe1 = matrix(1,NNBz,1,NNBt);
	Bwe2 = matrix(1,NNBz,1,NNBt);
	Bsq1 = matrix(1,NNBz,1,NNBt);
	Bsq2 = matrix(1,NNBz,1,NNBt);

	for(j=1; j<=NNBt; j++){
		for(i=2; i<NNBz; i++){
			// aN, aS, aW, aE
			aN[i][j] = pow(Hn[i][j],3)*Rp*dT[i][j]/dN[i][j];
			aS[i][j] = pow(Hs[i][j],3)*Rp*dT[i][j]/dS[i][j];
			aW[i][j] = pow(Hw[i][j],3)*Rp*dZ[i][j]/dW[i][j];
			aE[i][j] = pow(He[i][j],3)*Rp*dZ[i][j]/dE[i][j];

			aC[i][j] = aN[i][j] + aS[i][j] + aW[i][j] + aE[i][j]; //aC

			//if(dN[i][j] == 0){ /*aC[i][j] = TINY;*/ 			
				//printf("%e %e %e %e %e \n", Hn[i][j], Hs[i][j], dN[i][j], dS[i][j], aC[i][j]);
			//	exit(1);
			//}
			
			// Bwe, Bsq
			Bwe1[i][j] = Rp*(Hw[i][j]-He[i][j])*dZ[i][j];  // wedge
			Bwe2[i][j] = Rp*(Hn[i][j]-Hs[i][j])*dT[i][j]; 
			Bsq1[i][j] = Rp*( (cos(Te[i][j]+kai)-cos(Tw[i][j]+kai))*de-e*(sin(Te[i][j]+kai)-sin(Tw[i][j]+kai))*dkai )*dZ[i][j]; // Normal squeeze
			Bsq2[i][j] = (pow(Zn[i][j],2)-pow(Zs[i][j],2))/2*Rp*((sin(Te[i][j]+zeta)-sin(Tw[i][j]+zeta))*dgamma
										-(cos(Te[i][j]+zeta)-cos(Tw[i][j]+zeta))*gamma*dzeta);  // Titing squeeze

			B[i][j] = (w+wp)/w*Bwe1[i][j] - Vp*Bwe2[i][j] /*+ 2*(Bsq1[i][j] + Bsq2[i][j])*/;  // b

			//printf("%e %e \n", aC[i][j], aC[i][j]);
		}
	}
	//float temp=0.0;
	//for(j=1; j<=NNBt; j+=1) 
	//	for(i=1; i<=NNBz; i++) 
	//		//if(P[i][j] > 1.0) 
	//		temp = temp + Bsq2[i][j]; 	

	free_matrix(Bwe1,1,NNBz,1,NNBt);
	free_matrix(Bwe2,1,NNBz,1,NNBt);
	free_matrix(Bsq1,1,NNBz,1,NNBt);
	free_matrix(Bsq2,1,NNBz,1,NNBt);
}
void solver(double phi, double *Pc, float **P, float **aN, float **aS, float **aW, float **aE, float **aC, float **B)
{
	void equation(float **P, float **Pold, float **aN, float **aS, float **aW, float **aE, float **aC, float **B);
	void boundary(double phi, float **P, double *Pc);
	void error(float **P, float **Pold, float *Psum, float *Esum);

	int i, j, inumber;
	float **Pold, Psum, Esum;

	Pold = matrix(1,NNBz,1,NNBt);

	/* SOR */
	inumber = 0;
	float es = efs;
	//if(phi*180/Pi >45) es = 0.000002;
	//if(phi*180/Pi >55) es = 0.0000001;
	
	do{		 
		for(j=1; j<=NNBt; j++){				// Substitute OLD VALUE
			for(i=2; i<NNBz; i++){
				Pold[i][j] = P[i][j];
			}
		}		
		inumber = inumber + 1;				// Iteration number
		// Re-calculation
		equation(P, Pold, aN, aS, aW, aE, aC, B);
		//boundary(phi, P, Pc);	// Maintain Boundary Condition					
		error(P, Pold, &Psum, &Esum);
		//printf("================== %e  %d \n", Esum/Psum, inumber);
		//if(inumber > 6e3) es = es*1.5;
	}while((Esum/Psum) >= es);	
	printf("================== %e  %d   %e\n", Esum/Psum, inumber, es);
	//es = efs;

	free_matrix(Pold,1,NNBz,1,NNBt);
}
void equation(float **P, float **Pold, float **aN, float **aS, float **aW, float **aE, float **aC, float **B)
{
	int i,j;
	float W = 1.7;

	////	j = 0.75*NNBt ~ NNBt-1
	//for(j=int(NNBt*0.75); j<NNBt; j++){
	//	for(i=2; i<=1+NBz1+NBz2/2; i++){
	//		P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
	//					+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
	//	}
	//	for(i=NNBz-1; i>1+NBz1+NBz2/2; i--){
	//		P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
	//					+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
	//	}
	//}
	////	j = NNBt	
	//for(i=2; i<=1+NBz1+NBz2/2; i++){
	//	P[i][NNBt] = (1-W)*Pold[i][NNBt] + W*(aN[i][NNBt]*P[i+1][NNBt] + aS[i][NNBt]*P[i-1][NNBt] + aW[i][NNBt]*P[i][NNBt-1]   
	//					+ aE[i][NNBt]*P[i][1] - B[i][NNBt]) / aC[i][NNBt];
	//}
	//for(i=NNBz-1; i>1+NBz1+NBz2/2; i--){
	//	P[i][NNBt] = (1-W)*Pold[i][NNBt] + W*(aN[i][NNBt]*P[i+1][NNBt] + aS[i][NNBt]*P[i-1][NNBt] + aW[i][NNBt]*P[i][NNBt-1]   
	//					+ aE[i][NNBt]*P[i][1] - B[i][NNBt]) / aC[i][NNBt];
	//}
	////for(i=2; i<NNBz; i++){
	////	P[i][NNBt] = (1-W)*Pold[i][NNBt] + W*(aN[i][NNBt]*P[i+1][NNBt] + aS[i][NNBt]*P[i-1][NNBt] + aW[i][NNBt]*P[i][NNBt-1]   
	////					+ aE[i][NNBt]*P[i][1] - B[i][NNBt]) / aC[i][NNBt];
	////}
	//
	//// j = 1
	//for(i=2; i<=1+NBz1+NBz2/2; i++){
	//	P[i][1] = (1-W)*Pold[i][1] + W*(aN[i][1]*P[i+1][1] + aS[i][1]*P[i-1][1] + aW[i][1]*P[i][NNBt]   
	//				+ aE[i][1]*P[i][2] - B[i][1]) / aC[i][1];
	//}
	//for(i=NNBz-1; i>1+NBz1+NBz2/2; i--){
	//	P[i][1] = (1-W)*Pold[i][1] + W*(aN[i][1]*P[i+1][1] + aS[i][1]*P[i-1][1] + aW[i][1]*P[i][NNBt]   
	//				+ aE[i][1]*P[i][2] - B[i][1]) / aC[i][1];
	//}
	////for(i=2; i<NNBz; i++){
	////	P[i][1] = (1-W)*Pold[i][1] + W*(aN[i][1]*P[i+1][1] + aS[i][1]*P[i-1][1] + aW[i][1]*P[i][NNBt]   
	////				+ aE[i][1]*P[i][2] - B[i][1]) / aC[i][1];
	////}

	////	j = 2 ~ 0.75*NNBt
	//for(j=2; j<int(NNBt*0.75); j++){
	//	for(i=2; i<=1+NBz1+NBz2/2; i++){
	//		P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
	//					+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
	//	}
	//	for(i=NNBz-1; i>1+NBz1+NBz2/2; i--){
	//		P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
	//					+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
	//	}
	//}	

	////	j = 2~NNBt-1
	////for(j=2; j<NNBt; j++){
	////	for(i=2; i<=1+NBz1+NBz2/2; i++){
	////		P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
	////					+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
	////	}
	////	for(i=NNBz-1; i>1+NBz1+NBz2/2; i--){
	////		P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
	////					+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
	////	}
	////}	

	// j = 2~NNBt-1
	for(j=2; j<NNBt; j++){	
		for(i=2; i<NNBz; i++){
			P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
					+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
			if(P[i][j]<Pca/Lambda) P[i][j] = Pca/Lambda; 
			//printf("%e %e \n", aN[i][j], B[i][j]);
		}		
	}
	// j = 1
	for(i=2; i<NNBz; i++){
		P[i][1] = (1-W)*Pold[i][1] + W*(aN[i][1]*P[i+1][1] + aS[i][1]*P[i-1][1] + aW[i][1]*P[i][NNBt]   
					+ aE[i][1]*P[i][2] - B[i][1]) / aC[i][1];
		if(P[i][j]<Pca/Lambda) P[i][j] = Pca/Lambda; 
	}
	//	j = NNBt	
	for(i=2; i<NNBz; i++){
		P[i][NNBt] = (1-W)*Pold[i][NNBt] + W*(aN[i][NNBt]*P[i+1][NNBt] + aS[i][NNBt]*P[i-1][NNBt] + aW[i][NNBt]*P[i][NNBt-1]   
						+ aE[i][NNBt]*P[i][1] - B[i][NNBt]) / aC[i][NNBt];
		if(P[i][j]<Pca/Lambda) P[i][j] = Pca/Lambda; 
	}
	

	////	pressure should be bigger than 0
	//for(j=1; j<NNBt; j++){
	//	for(i=2; i<NNBz; i++){
	//		if(P[i][j]<=Pca/Lambda) P[i][j] = Pca/Lambda;       // Pca (MPa) = ABSOLUTE pressure
	//	}
	//}
}
/* Boundary Condition */
void boundary(double phi, float **Ppc, double *Pc)
{
	int i,j,k;
	int NBKt, dNBKt;

	///* Valve Port pressure */
	//for(j=1+NVt1-int(phi/dNBt); j<=1+NVt1+int(rV1/R/dNBt)-int(phi/dNBt); j++){
	//	for(i=(NVr1+NVr2)/2-int(sqrt(rV1*rV1-pow(R*dNBt*(j-(1+NVt1-int(phi/dNBt)))-rV1,2))*NBz2/(zB3-zB2)/Rbo); 
	//		i<=(NVr1+NVr2)/2+int(sqrt(rV1*rV1-pow(R*dNBt*(j-(1+NVt1-int(phi/dNBt)))-rV1,2))*NBz2/(zB3-zB2)/Rbo); i++){
	//		k = j;
	//		if(j<1) k = j + NNBt;
	//		Ppc[i][k] = Pc[N+1];   // pressure in discharge valve port
	//	}
	//}
	//for(j=1+NVt1+int(rV1/R/dNBt)-int(phi/dNBt); j<=1+NVt1+NVt2-int(rV1/R/dNBt)-int(phi/dNBt); j++){
	//	for(i=NVr1; i<=NVr2; i++){
	//		k = j;
	//		if(j<1) k = j + NNBt;
	//		Ppc[i][k] = Pc[N+1];   // pressure in discharge valve port
	//	}
	//}
	//for(j=1+NVt1+NVt2-int(rV1/R/dNBt)-int(phi/dNBt); j<=1+NVt1+NVt2-int(phi/dNBt); j++){
	//	for(i=(NVr1+NVr2)/2-int(sqrt(rV1*rV1-pow(R*dNBt*(j-(1+NVt1+NVt2-int(rV1/R/dNBt)-int(phi/dNBt))),2))*NBz2/(zB3-zB2)/Rbo); 
	//		i<=(NVr1+NVr2)/2+int(sqrt(rV1*rV1-pow(R*dNBt*(j-(1+NVt1+NVt2-int(rV1/R/dNBt)-int(phi/dNBt))),2))*NBz2/(zB3-zB2)/Rbo); i++){
	//		k = j;
	//		if(j<1) k = j + NNBt;
	//		Ppc[i][k] = Pc[N+1];   // pressure in discharge valve port
	//	}
	//}
	//for(j=1+NVt1+NVt2+NVt3+NVt4-int(phi/dNBt); j<=1+NVt1+NVt2+NVt3+NVt4+int(rV2/R/dNBt)-int(phi/dNBt); j++){
	//	for(i=(NVr3+NVr4)/2-int(sqrt(rV2*rV2-pow(R*dNBt*(j-(1+NVt1+NVt2+NVt3+NVt4-int(phi/dNBt)))-rV2,2))*NBz2/(zB3-zB2)/Rbo); 
	//		i<=(NVr3+NVr4)/2+int(sqrt(rV2*rV2-pow(R*dNBt*(j-(1+NVt1+NVt2+NVt3+NVt4-int(phi/dNBt)))-rV2,2))*NBz2/(zB3-zB2)/Rbo); i++){
	//		k = j;
	//		if(j<1) k = j + NNBt;
	//		Ppc[i][k] = Pc[N+2];	// pressure in suction valve port
	//	}
	//}
	//	// Valve rounding
	//for(j=1+NVt1+NVt2+NVt3+NVt4+int(rV2/R/dNBt)-int(phi/dNBt); j<=1+NVt1+NVt2+NVt3+NVt4+NVt5-int(rV2/R/dNBt)-int(phi/dNBt); j++){
	//	for(i=NVr3; i<=NVr4; i++){
	//		k = j;
	//		if(j<1) k = j + NNBt;
	//		Ppc[i][k] = Pc[N+2];	// pressure in suction valve port
	//	}
	//}
	//for(j=1+NVt1+NVt2+NVt3+NVt4+NVt5-int(rV2/R/dNBt)-int(phi/dNBt); j<=1+NVt1+NVt2+NVt3+NVt4+NVt5-int(phi/dNBt); j++){
	//	for(i=(NVr3+NVr4)/2-int(sqrt(rV2*rV2-pow(R*dNBt*(j-(1+NVt1+NVt2+NVt3+NVt4+NVt5-int(rV2/R/dNBt)-int(phi/dNBt))),2))*NBz2/(zB3-zB2)/Rbo); 
	//		i<=(NVr3+NVr4)/2+int(sqrt(rV2*rV2-pow(R*dNBt*(j-(1+NVt1+NVt2+NVt3+NVt4+NVt5-int(rV2/R/dNBt)-int(phi/dNBt))),2))*NBz2/(zB3-zB2)/Rbo); i++){
	//		k = j;
	//		if(j<1) k = j + NNBt;
	//		Ppc[i][k] = Pc[N+2];	// pressure in suction valve port
	//	}
	//}
	//	
	///* Barrel Kidney(cylinder) pressure */
	//NBKt = int(thetaK/dNBt+1.5);		// Barrel Kidney 각도방향 grid 수
	//dNBKt = int(NNBt/N+0.5);			// Barrel Kidney간격 각도방향 grid 수
	//for(i=1; i<N; i++){
	//	for(k=NBz1+1; k<=NBz1+1+NBz2; k++){
	//		for(j=1+dNBKt*(i-1); j<=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j++){			
	//			Ppc[k][j] = Pc[i];
	//		}			
	//		for(j=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j<=dNBKt*i; j++){
	//			Ppc[k][j] = Pc[i+1];
	//		}
	//	}
	//	for(j=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(i-1); j++){
	//		for(k=1+NBz1+NBz2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBz2/(zB3-zB2)/Rbo); 
	//			k<=1+NBz1+NBz2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBz2/(zB3-zB2)/Rbo); k++){
	//				Ppc[k][j] = Pc[i];
	//		}
	//	}
	//	for(j=1+dNBKt*i-NBKt/2; j<=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j++){
	//		for(k=1+NBz1+NBz2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBz2/(zB3-zB2)/Rbo); 
	//			k<=1+NBz1+NBz2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBz2/(zB3-zB2)/Rbo); k++){
	//				Ppc[k][j] = Pc[i+1];
	//		}
	//	}
	//}
	//	// Barrel rounding
	//for(k=NBz1+1; k<=NBz1+1+NBz2; k++){
	//	for(j=1+dNBKt*(N-1); j<=1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt); j++){			
	//		Ppc[k][j] = Pc[N];
	//	}
	//	for(j=1+dNBKt*N-NBKt/2+int(rB/R/dNBt); j<=dNBKt*N; j++){
	//		Ppc[k][j] = Pc[1];
	//	}
	//}
	//for(j=1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(N-1); j++){
	//	for(k=1+NBz1+NBz2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt))),2))*NBz2/(zB3-zB2)/Rbo); 
	//		k<=1+NBz1+NBz2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt))),2))*NBz2/(zB3-zB2)/Rbo); k++){
	//			Ppc[k][j] = Pc[N];
	//	}
	//}
	//for(j=1+dNBKt*N-NBKt/2; j<=1+dNBKt*N-NBKt/2+int(rB/R/dNBt); j++){
	//	for(k=1+NBz1+NBz2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*N-NBKt/2))-rB,2))*NBz2/(zB3-zB2)/Rbo); 
	//		k<=1+NBz1+NBz2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*N-NBKt/2))-rB,2))*NBz2/(zB3-zB2)/Rbo); k++){
	//			Ppc[k][j] = Pc[1];
	//	}
	//}
	


}
void error(float **P, float **Pold, float *Psum, float *Esum)
{
	int i,j;

	*Esum=0.0;
	*Psum=0.0;
	for(j=1; j<=NNBt; j++){
		for(i=2; i<NNBz; i++){
			*Esum = *Esum + fabs(P[i][j]-Pold[i][j]);
			*Psum = *Psum + fabs(P[i][j]);
		}
	}
}