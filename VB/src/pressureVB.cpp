#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

extern int Nphi, NBr1, NBr2, NBr3, NBr4, NBr5, NBr, NNBr, NNBt, Nch, Imin, Jmin;
extern int NVt1, NVt2, NVt3, NVt4, NVt5, NVt6, NVr1, NVr2, NVr3, NVr4;
extern float dNBt, thetach;
extern float R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, revf, phif, phiout_c, phiout_vb, phiout_FM, phiout_LK, w;
extern float V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
extern float integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, rB1, rB2, rB3, rB4, rB5, rB6, Fs, hv0, theta0, psi0, dhv0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
extern float Lambda, Loc, Lpj, Lpg, Lbj, Rbo;
extern float FpressD, MpressxD, MpressyD, FliftD, MliftxD, MliftyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, FzD, MxD, MyD, MzD, leakageRinD, leakageRoutD, Hmin;
double gamH, zetaH;

/* Film thickness */
void filmthickness(double phi, double *ym, double **H, double **Hn, double **Hs, double **Hw, double **He, double **H_gro, float **Rvb, float **T, double gamma, double zeta)
{
	int i,j;
	//gamH = 0.005*Pi/180*Rbo/hv0;
	//zetaH = 40.0*Pi/180;

	// Solid contact
	//if(ym[1]-gamma <= hvbmin/hv0){
	//	printf("Minus Film Thickness < Hmin at Hmin = %f \n", (ym[1]-gamma)*hv0);	
	//	//gamma = ym[1] - hvbmin/hv0;
	//	for(j=1; j<=NNBt; j++){
	//		for(i=1; i<=NNBr; i++){
	//			H[i][j] = ym[1] - Rvb[i][j]*gamma*cos(T[i][j]-(Pi/2-(phi+zeta)))  + H_gro[i][j];
	//		
	//		}
	//	}		
	//}
	//else{
		//	FILM THICKNESS
		Hmin = H[1][1];
		for(j=1; j<=NNBt; j++){
			for(i=1; i<=NNBr; i++){
				H[i][j] = ym[1] - Rvb[i][j]*gamma*cos(T[i][j]-(Pi/2-(phi+zeta))) + H_gro[i][j];
				if(Hmin > H[i][j]){
					Hmin = H[i][j];
					Imin = i;
					Jmin = j;
				}
			}			 
		}
	//}

		
	for(j=2; j<NNBt; j++){
		for(i=2; i<NNBr; i++){
			Hn[i][j] = (H[i][j]+H[i+1][j])/2.0;
			Hs[i][j] = (H[i][j]+H[i-1][j])/2.0;
			Hw[i][j] = (H[i][j]+H[i][j-1])/2.0;
			He[i][j] = (H[i][j]+H[i][j+1])/2.0;
		}
	}
	// j=1
	for(i=2; i<NNBr; i++){
		Hn[i][1] = (H[i][1]+H[i+1][1])/2.0;
		Hs[i][1] = (H[i][1]+H[i-1][1])/2.0;
		Hw[i][1] = (H[i][1]+H[i][NNBt])/2.0;
		He[i][1] = (H[i][1]+H[i][2])/2.0;
	}
	// j=NNBt
	for(i=2; i<NNBr; i++){
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
}
/* FVM coefficient */
void descretize(double phi, double *ym, double gamma, double zeta, double dgamma, double dzeta, float **Rvb, float **dT, float **dR, double **Hn, double **Hs, double **Hw, double **He, 
				float **Rn, float **Rs, float **Rw, float **Re, float **Tn, float **Ts, float **Tw, float **Te, float **dN, float **dS, float **dW, float **dE, 
				float **aN, float **aS, float **aW, float **aE, float **aC, float **B)
{
	int i,j;
	float **Bwe, **Bsq1, **Bsq2;
	Bwe = matrix(1,NNBr,1,NNBt);
	Bsq1 = matrix(1,NNBr,1,NNBt);
	Bsq2 = matrix(1,NNBr,1,NNBt);

	for(j=1; j<=NNBt; j++){
		for(i=2; i<NNBr; i++){
			// aN, aS, aW, aE
			aN[i][j] = pow(Hn[i][j],3)*Rn[i][j]*dT[i][j]/dN[i][j];
			aS[i][j] = pow(Hs[i][j],3)*Rs[i][j]*dT[i][j]/dS[i][j];
			aW[i][j] = pow(Hw[i][j],3)/Rw[i][j]*dR[i][j]/dW[i][j];
			aE[i][j] = pow(He[i][j],3)/Re[i][j]*dR[i][j]/dE[i][j];

			aC[i][j] = aN[i][j] + aS[i][j] + aW[i][j] + aE[i][j]; //aC

			//if(aC[i][j] < 1.0e-18) aC[i][j] = TINY; 
			
			// Bwe, Bsq
			Bwe[i][j] = Rvb[i][j]*(He[i][j]-Hw[i][j])*dR[i][j];  // wedge
			Bsq1[i][j] = ym[4]*dT[i][j]*(Rn[i][j]*Rn[i][j]-Rs[i][j]*Rs[i][j])/2; // Normal squeeze
			Bsq2[i][j] = (pow(Rn[i][j],3)-pow(Rs[i][j],3))/3*((sin(Te[i][j]-(Pi/2-(phi+zeta)))-sin(Tw[i][j]-(Pi/2-(phi+zeta))))*dgamma
										+(cos(Te[i][j]-(Pi/2-(phi+zeta)))-cos(Tw[i][j]-(Pi/2-(phi+zeta))))*gamma*dzeta
										/*-(cos(Te[i][j]-(Pi/2-(phi+zetaH)))-cos(Tw[i][j]-(Pi/2-(phi+zetaH))))*gamH*0.0*/);  // Geometric squeeze

			B[i][j] = -Bwe[i][j] + 2*(Bsq1[i][j] - Bsq2[i][j]);  // b
		}
	}
	//float temp=0.0;
	//for(j=1; j<=NNBt; j+=1) 
	//	for(i=1; i<=NNBr; i++) 
	//		//if(P[i][j] > 1.0) 
	//		temp = temp + Bsq2[i][j]; 				
	//printf("%e \n", temp);

	free_matrix(Bwe,1,NNBr,1,NNBt);
	free_matrix(Bsq1,1,NNBr,1,NNBt);
	free_matrix(Bsq2,1,NNBr,1,NNBt);
}
void solver(double phi, double *Pc, float **P, float **aN, float **aS, float **aW, float **aE, float **aC, float **B)
{
	void equation(float **P, float **Pold, float **aN, float **aS, float **aW, float **aE, float **aC, float **B);
	void boundary(double phi, float **P, double *Pc);
	void error(float **P, float **Pold, float *Psum, float *Esum);

	int i, j, inumber;
	float **Pold, Psum, Esum;

	Pold = matrix(1,NNBr,1,NNBt);

	/* SOR */
	inumber = 0;
	float es = efs;
	//if(phi*180/Pi >45) es = 0.000002;
	//if(phi*180/Pi >55) es = 0.0000001;
	
	do{		 
		for(j=1; j<=NNBt; j++){				// Substitute OLD VALUE
			for(i=2; i<NNBr; i++){
				Pold[i][j] = P[i][j];
			}
		}		
		inumber = inumber + 1;				// Iteration number
		// Re-calculation
		equation(P, Pold, aN, aS, aW, aE, aC, B);
		boundary(phi, P, Pc);	// Maintain Boundary Condition					
		error(P, Pold, &Psum, &Esum);
		//if(inumber > 4e3) es = es*1.5;
	}while((Esum/Psum) >= es);
	//float temp=0.0;
	//for(j=1; j<=NNBt; j+=1) 
	//	for(i=1; i<=NNBr; i++) 
	//		//if(P[i][j] > 1.0) 
	//		temp = temp + P[i][j]; 				
	printf("================== %e  %d   %e\n", Esum/Psum, inumber, es);
	//es = efs;

	free_matrix(Pold,1,NNBr,1,NNBt);
}
void equation(float **P, float **Pold, float **aN, float **aS, float **aW, float **aE, float **aC, float **B)
{
	int i,j;
	float W = 1.8;

	//	j = 0.75*NNBt ~ NNBt-1
	for(j=int(NNBt*0.75); j<NNBt; j++){
		for(i=2; i<=1+NBr1+NBr2/2; i++){
			P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
						+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
		}
		for(i=NNBr-1; i>1+NBr1+NBr2/2; i--){
			P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
						+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
		}
	}
	//	j = NNBt	
	for(i=2; i<=1+NBr1+NBr2/2; i++){
		P[i][NNBt] = (1-W)*Pold[i][NNBt] + W*(aN[i][NNBt]*P[i+1][NNBt] + aS[i][NNBt]*P[i-1][NNBt] + aW[i][NNBt]*P[i][NNBt-1]   
						+ aE[i][NNBt]*P[i][1] - B[i][NNBt]) / aC[i][NNBt];
	}
	for(i=NNBr-1; i>1+NBr1+NBr2/2; i--){
		P[i][NNBt] = (1-W)*Pold[i][NNBt] + W*(aN[i][NNBt]*P[i+1][NNBt] + aS[i][NNBt]*P[i-1][NNBt] + aW[i][NNBt]*P[i][NNBt-1]   
						+ aE[i][NNBt]*P[i][1] - B[i][NNBt]) / aC[i][NNBt];
	}
	//for(i=2; i<NNBr; i++){
	//	P[i][NNBt] = (1-W)*Pold[i][NNBt] + W*(aN[i][NNBt]*P[i+1][NNBt] + aS[i][NNBt]*P[i-1][NNBt] + aW[i][NNBt]*P[i][NNBt-1]   
	//					+ aE[i][NNBt]*P[i][1] - B[i][NNBt]) / aC[i][NNBt];
	//}
	
	// j = 1
	for(i=2; i<=1+NBr1+NBr2/2; i++){
		P[i][1] = (1-W)*Pold[i][1] + W*(aN[i][1]*P[i+1][1] + aS[i][1]*P[i-1][1] + aW[i][1]*P[i][NNBt]   
					+ aE[i][1]*P[i][2] - B[i][1]) / aC[i][1];
	}
	for(i=NNBr-1; i>1+NBr1+NBr2/2; i--){
		P[i][1] = (1-W)*Pold[i][1] + W*(aN[i][1]*P[i+1][1] + aS[i][1]*P[i-1][1] + aW[i][1]*P[i][NNBt]   
					+ aE[i][1]*P[i][2] - B[i][1]) / aC[i][1];
	}
	//for(i=2; i<NNBr; i++){
	//	P[i][1] = (1-W)*Pold[i][1] + W*(aN[i][1]*P[i+1][1] + aS[i][1]*P[i-1][1] + aW[i][1]*P[i][NNBt]   
	//				+ aE[i][1]*P[i][2] - B[i][1]) / aC[i][1];
	//}

	//	j = 2 ~ 0.75*NNBt
	for(j=2; j<int(NNBt*0.75); j++){
		for(i=2; i<=1+NBr1+NBr2/2; i++){
			P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
						+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
		}
		for(i=NNBr-1; i>1+NBr1+NBr2/2; i--){
			P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
						+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
		}
	}	

	//	j = 2~NNBt-1
	//for(j=2; j<NNBt; j++){
	//	for(i=2; i<=1+NBr1+NBr2/2; i++){
	//		P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
	//					+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
	//	}
	//	for(i=NNBr-1; i>1+NBr1+NBr2/2; i--){
	//		P[i][j] = (1-W)*Pold[i][j] + W*(aN[i][j]*P[i+1][j] + aS[i][j]*P[i-1][j]	+ aW[i][j]*P[i][j-1]
	//					+ aE[i][j]*P[i][j+1] - B[i][j]) / aC[i][j];
	//	}
	//}	

	

	//	Pressure should be bigger than 0
	for(j=1; j<NNBt; j++){
		for(i=2; i<NNBr; i++){
			if(P[i][j]<=Pca/Lambda) P[i][j] = Pca/Lambda;       // Pca (MPa) = ABSOLUTE PRESSURE
		}
	}
}
/* Boundary Condition */
void boundary(double phi, float **Pvb, double *Pc)
{
	int i,j,k;
	int NBKt, dNBKt;

	/* Valve Port pressure */
	for(j=1+NVt1-int(phi/dNBt); j<=1+NVt1+int(rV1/R/dNBt)-int(phi/dNBt); j++){
		for(i=(NVr1+NVr2)/2-int(sqrt(rV1*rV1-pow(R*dNBt*(j-(1+NVt1-int(phi/dNBt)))-rV1,2))*NBr2/(rB3-rB2)/Rbo); 
			i<=(NVr1+NVr2)/2+int(sqrt(rV1*rV1-pow(R*dNBt*(j-(1+NVt1-int(phi/dNBt)))-rV1,2))*NBr2/(rB3-rB2)/Rbo); i++){
			k = j;
			if(j<1) k = j + NNBt;
			Pvb[i][k] = Pc[N+1];   // pressure in discharge valve port
		}
	}
	for(j=1+NVt1+int(rV1/R/dNBt)-int(phi/dNBt); j<=1+NVt1+NVt2-int(rV1/R/dNBt)-int(phi/dNBt); j++){
		for(i=NVr1; i<=NVr2; i++){
			k = j;
			if(j<1) k = j + NNBt;
			Pvb[i][k] = Pc[N+1];   // pressure in discharge valve port
		}
	}
	for(j=1+NVt1+NVt2-int(rV1/R/dNBt)-int(phi/dNBt); j<=1+NVt1+NVt2-int(phi/dNBt); j++){
		for(i=(NVr1+NVr2)/2-int(sqrt(rV1*rV1-pow(R*dNBt*(j-(1+NVt1+NVt2-int(rV1/R/dNBt)-int(phi/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); 
			i<=(NVr1+NVr2)/2+int(sqrt(rV1*rV1-pow(R*dNBt*(j-(1+NVt1+NVt2-int(rV1/R/dNBt)-int(phi/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); i++){
			k = j;
			if(j<1) k = j + NNBt;
			Pvb[i][k] = Pc[N+1];   // pressure in discharge valve port
		}
	}
	for(j=1+NVt1+NVt2+NVt3+NVt4-int(phi/dNBt); j<=1+NVt1+NVt2+NVt3+NVt4+int(rV2/R/dNBt)-int(phi/dNBt); j++){
		for(i=(NVr3+NVr4)/2-int(sqrt(rV2*rV2-pow(R*dNBt*(j-(1+NVt1+NVt2+NVt3+NVt4-int(phi/dNBt)))-rV2,2))*NBr2/(rB3-rB2)/Rbo); 
			i<=(NVr3+NVr4)/2+int(sqrt(rV2*rV2-pow(R*dNBt*(j-(1+NVt1+NVt2+NVt3+NVt4-int(phi/dNBt)))-rV2,2))*NBr2/(rB3-rB2)/Rbo); i++){
			k = j;
			if(j<1) k = j + NNBt;
			Pvb[i][k] = Pc[N+2];	// pressure in suction valve port
		}
	}
		// Valve rounding
	for(j=1+NVt1+NVt2+NVt3+NVt4+int(rV2/R/dNBt)-int(phi/dNBt); j<=1+NVt1+NVt2+NVt3+NVt4+NVt5-int(rV2/R/dNBt)-int(phi/dNBt); j++){
		for(i=NVr3; i<=NVr4; i++){
			k = j;
			if(j<1) k = j + NNBt;
			Pvb[i][k] = Pc[N+2];	// pressure in suction valve port
		}
	}
	for(j=1+NVt1+NVt2+NVt3+NVt4+NVt5-int(rV2/R/dNBt)-int(phi/dNBt); j<=1+NVt1+NVt2+NVt3+NVt4+NVt5-int(phi/dNBt); j++){
		for(i=(NVr3+NVr4)/2-int(sqrt(rV2*rV2-pow(R*dNBt*(j-(1+NVt1+NVt2+NVt3+NVt4+NVt5-int(rV2/R/dNBt)-int(phi/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); 
			i<=(NVr3+NVr4)/2+int(sqrt(rV2*rV2-pow(R*dNBt*(j-(1+NVt1+NVt2+NVt3+NVt4+NVt5-int(rV2/R/dNBt)-int(phi/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); i++){
			k = j;
			if(j<1) k = j + NNBt;
			Pvb[i][k] = Pc[N+2];	// pressure in suction valve port
		}
	}
		
	/* Barrel Kidney(cylinder) pressure */
	NBKt = int(thetaK/dNBt+1.5);		// Barrel Kidney 각도방향 grid 수
	dNBKt = int(NNBt/N+0.5);			// Barrel Kidney간격 각도방향 grid 수
	for(i=1; i<N; i++){
		for(k=NBr1+1; k<=NBr1+1+NBr2; k++){
			for(j=1+dNBKt*(i-1); j<=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j++){			
				Pvb[k][j] = Pc[i];
			}			
			for(j=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j<=dNBKt*i; j++){
				Pvb[k][j] = Pc[i+1];
			}
		}
		for(j=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(i-1); j++){
			for(k=1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); 
				k<=1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); k++){
					Pvb[k][j] = Pc[i];
			}
		}
		for(j=1+dNBKt*i-NBKt/2; j<=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j++){
			for(k=1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo); 
				k<=1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo); k++){
					Pvb[k][j] = Pc[i+1];
			}
		}
	}
		// Barrel rounding
	for(k=NBr1+1; k<=NBr1+1+NBr2; k++){
		for(j=1+dNBKt*(N-1); j<=1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt); j++){			
			Pvb[k][j] = Pc[N];
		}
		for(j=1+dNBKt*N-NBKt/2+int(rB/R/dNBt); j<=dNBKt*N; j++){
			Pvb[k][j] = Pc[1];
		}
	}
	for(j=1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(N-1); j++){
		for(k=1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); 
			k<=1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); k++){
				Pvb[k][j] = Pc[N];
		}
	}
	for(j=1+dNBKt*N-NBKt/2; j<=1+dNBKt*N-NBKt/2+int(rB/R/dNBt); j++){
		for(k=1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*N-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo); 
			k<=1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*N-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo); k++){
				Pvb[k][j] = Pc[1];
		}
	}
	

	//// Barrel Kidney boundary - circumferential
	//for(int i=1; i<=N; i++){
	//	for(j=1+dNBKt*(i-1); j<=1+NBKt/2+dNBKt*(i-1); j++){
	//		P[NBr1+1][j] = Pc[i];
	//		P[NBr1+1+NBr2][j] = Pc[i];
	//	}
	//}
	//for(i=1; i<N; i++){
	//	for(j=1+dNBKt*i-NBKt/2; j<=dNBKt*i; j++){
	//		P[NBr1+1][j] = Pc[i+1];
	//		P[NBr1+1+NBr2][j] = Pc[i+1];
	//	}		
	//}
	//for(j=1+dNBKt*N-NBKt/2; j<=dNBKt*N; j++){
	//	P[NBr1+1][j] = Pc[1];
	//	P[NBr1+1+NBr2][j] = Pc[1];
	//}
	////Barrel Kidney boundary - radial
	//for(int i=1; i<N; i++){
	//	for(j=NBr1+1; j<=NBr1+1+NBr2; j++){
	//		P[j][1+NBKt/2+dNBKt*(i-1)] = Pc[i];
	//		P[j][1+dNBKt*i-NBKt/2] = Pc[i+1];
	//	}
	//}
	//for(j=NBr1+1; j<=NBr1+1+NBr2; j++){
	//	P[j][1+NBKt/2+dNBKt*(N-1)] = Pc[N];
	//	P[j][1+dNBKt*N-NBKt/2] = Pc[1];
	//}		
}
void error(float **P, float **Pold, float *Psum, float *Esum)
{
	int i,j;

	*Esum=0.0;
	*Psum=0.0;
	for(j=1; j<=NNBt; j++){
		for(i=2; i<NNBr; i++){
			*Esum = *Esum + fabs(P[i][j]-Pold[i][j]);
			*Psum = *Psum + fabs(P[i][j]);
		}
	}
}