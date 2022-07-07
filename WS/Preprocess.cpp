#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

extern int Nphi, NBr1, NBr2, NBr3, NBr4, NBr5, NBr, NNBr, NNBt, Nch;
extern int NVt1, NVt2, NVt3, NVt4, NVt5, NVt6, NVr1, NVr2, NVr3, NVr4;
extern double dNBt, thetach;
extern double R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, revf, phif, phiout_c, hm, phiout_cy, phiout_dy, phiout_vb, phiout_FM, phiout_LK, w;
extern double V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
extern double integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, rB1, rB2, rB3, rB4, rB5, rB6, Fs, hv0, theta0, psi0, dhv0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
extern double Lambda, Loc, Lpj, Lpg, Lbj, Ls, Rbo;
extern float FpressD, MpressxD, MpressyD, FliftD, MliftxD, MliftyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, FzD, MxD, MyD, MzD, leakageRinD, leakageRoutD;

/* input data */
void Input()
{
	FILE *ep;

	if((ep = fopen("Input_pressure_asiatribo.dat", "r"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}

	fscanf(ep, "%*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*c %lf %*c %lf %*c %lf %*c %lf %*c %lf %*s %lf %*s %lf %*c %lf %*s %lf %*c %lf	\
				%*s %lf %*c %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*c %lf %*s %lf %*s %lf %*c %lf %*s %lf %*c %lf %*s %lf %*c %lf %*s %lf %*s %lf %*s %lf	\
				%*s %lf %*s %lf %*c %lf %*c %lf %*c %lf %*c %lf %*s %lf %*s %d %*c %d %*c %d %*c %d %*c %d %*s %d %*s %lf %*c %lf %*c %lf %*c %lf %*c %lf %*c %lf %*s %d %*s %lf %*s %lf %*s %lf %*c %lf %*c %lf %*c %lf %*c %lf %*c %lf %*s %lf %*s %lf %*s %lf %*c %lf %*c %lf %*s %lf %*c %lf %*c %lf %*c %lf %*c %lf", 
		&R, &Rp, &Rpi, &Lpi, &Lc, &Ld, &thetaK, &thetaV1, &thetaV2, &thetaV3, &thetaV4, &thetaV5, &thetaV6, &rB, &rV1, &rV2, &LB1, &LB2, 
		&Cvo, &Cvi, &Cs, &Pin, &Pout, &Ph, &rho, &eta, &K, &Rs1, &Rs2, &Rso, &beta, &beta2, &Vvo, &Vvi, &Avout, &Avin, &Lvo, &Kvo, &Cvin, 
		&revf, &phiout_c, &hm, &phiout_cy, &phiout_dy, &phiout_vb, &w, &NBr1, &NBr2, &NBr3, &NBr4, &NBr5, &NNBt, &rB1, &rB2, &rB3, &rB4, &rB5, &rB6, &Nch, &thetach, 
		&Fs, &hv0, &theta0, &psi0, &dhv0, &wx0, &wy0, &Hg, &Mp, &Mb, &Iba, &Ibt, &Loc, &Lpj, &Lpg, &Lbj, &Ls);
	printf("%e %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e	\
			%e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e	\
			%e  %e  %e  %e  %e  %e  %e  %d  %d  %d  %d  %d  %d  %e  %e  %e  %e  %e  %e  %d  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e %e \n",
			R, Rp, Rpi, Lpi, Lc, Ld, thetaK, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6, rB, rV1, rV2, LB1, LB2, 
			Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, 
			revf, phiout_c, hm, phiout_cy, phiout_dy, phiout_vb, w, NBr1, NBr2, NBr3, NBr4, NBr5, NNBt, rB1, rB2, rB3, rB4, rB5, rB6, Nch, thetach, 
			Fs, hv0, theta0, psi0, dhv0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt, Loc, Lpj, Lpg, Lbj, Ls); 

	fclose(ep);

}
/* Parameters & Unit conversion & Non-dimensionalization */
void parameter()
{
	Ap = Pi*Rp*Rp;
	Vpi = Pi*Rpi*Rpi*Lpi;
	hp = 10.e-6;
	hs = 10.e-6;
	hv = hv0;	

	// NBr: number of CV's , NNBr: Number of Node
	NBr = NBr1 + NBr2 + NBr3 + NBr4 + NBr5;
	NNBr = NBr + 1;
	dNBt = 2*Pi/NNBt;

	rB = 0.5*(rB3-rB2);
	thetaK = thetaK*Pi/180;
	thetaB = thetaK/2;
	thetaV1 = thetaV1*Pi/180;
	thetaV2 = thetaV2*Pi/180;
	thetaV3 = thetaV3*Pi/180;
	thetaV4 = thetaV4*Pi/180;
	thetaV5 = thetaV5*Pi/180;
	thetaV6 = thetaV6*Pi/180;
	w = w*2*Pi/60;
	beta = beta*Pi/180;
	beta2 = beta2*Pi/180;
	phif = revf*2*Pi;
	theta0 = theta0*Pi/180;
	psi0 = psi0*Pi/180;
	hm = hm*Pi/180;

	NVt1 = int(NNBt*thetaV1/2/Pi+1.5); // The number of grids of Valve port 1st angle
	NVt2 = int(NNBt*thetaV2/2/Pi+1.5);
	NVt3 = int(NNBt*thetaV3/2/Pi+1.5);
	NVt4 = int(NNBt*thetaV4/2/Pi+1.5);
	NVt5 = int(NNBt*thetaV5/2/Pi+1.5);
	NVt6 = int(NNBt*thetaV6/2/Pi+1.5);
	NVr1 = 1+NBr1 + int(NBr2*((R-rV1)-rB2)/(rB3-rB2)+0.5); // Valve port inner radius in discharge port
	NVr2 = 1+NBr1 + int(NBr2*((R+rV1)-rB2)/(rB3-rB2)+0.5); // Valve port outer radius in discharge port
	NVr3 = 1+NBr1 + int(NBr2*((R-rV2)-rB2)/(rB3-rB2)+0.5); 
	NVr4 = 1+NBr1 + int(NBr2*((R+rV2)-rB2)/(rB3-rB2)+0.5);

	//printf("%d %d %d \n", NVt1, NVt5, NVr4);

	
	// Nondimensionalize Radii in the Barrel 
	Rbo = rB6;
	Lambda = 6*eta*w*pow(Rbo/hv0,2);
	rB1 = rB1/Rbo;
	rB2 = rB2/Rbo;
	rB3 = rB3/Rbo;
	rB4 = rB4/Rbo;
	rB5 = rB5/Rbo;
	rB6 = 1;
	Fs = Fs/Lambda/Rbo/Rbo;
	Hg = Hg/hv0;
	Mb = Mb*w*w/Lambda/Rbo/(Rbo/hv0);
	Mp = Mp*w*w/Lambda/Rbo/(Rbo/hv0);
	Iba = Iba*w*w/Lambda/pow(Rbo,3);
	Ibt = Ibt*w*w/Lambda/pow(Rbo,3);
	dhv0 = dhv0/hv0/w;
	wx0 = wx0*Pi/180/w;
	wy0 = wy0*Pi/180/w;

}
/* Grid & Coordinate */
void grid(float **Rvb, float **Rn, float **Rs, float **Rw, float **Re, float **T, float **Tn, float **Ts, float **Tw, float **Te, float **dR, float **dT, float **dN, float **dS, float **dW, float **dE)
{
	int i,j;	
		
	// dr1
	for(j=1; j<=NNBt; j++){
		for(i=2; i<NBr1+1; i++){	
			dN[i][j] = (rB2-rB1)/NBr1;
			dS[i][j] = (rB2-rB1)/NBr1;
			dW[i][j] = 2*Pi/NNBt;
			dE[i][j] = 2*Pi/NNBt;
			dR[i][j] = (dN[i][j]+dS[i][j])/2.0;
			dT[i][j] = (dW[i][j]+dE[i][j])/2.0;
		}	
		dN[1][j] = (rB2-rB1)/NBr1;
		dS[1][j] = TINY;
		dW[1][j] = 2*Pi/NNBt;
		dE[1][j] = 2*Pi/NNBt;
		dR[1][j] = (dN[1][j]+dS[1][j])/2.0;
		dT[1][j] = (dW[1][j]+dE[1][j])/2.0;
		dN[NBr1+1][j] = (rB3-rB2)/NBr2;
		dS[NBr1+1][j] = (rB2-rB1)/NBr1;
		dW[NBr1+1][j] = 2*Pi/NNBt;
		dE[NBr1+1][j] = 2*Pi/NNBt;
		dR[NBr1+1][j] = (dN[NBr1+1][j]+dS[NBr1+1][j])/2.0;
		dT[NBr1+1][j] = (dW[NBr1+1][j]+dE[NBr1+1][j])/2.0;		
	}		
	
	// dr2
	for(j=1; j<=NNBt; j++){
		for(i=1; i<NBr2; i++){
			dN[NBr1+1+i][j] = (rB3-rB2)/NBr2;
			dS[NBr1+1+i][j] = (rB3-rB2)/NBr2;
			dW[NBr1+1+i][j] = 2*Pi/NNBt;
			dE[NBr1+1+i][j] = 2*Pi/NNBt;
			dR[NBr1+1+i][j] = (dN[NBr1+1+i][j]+dS[NBr1+1+i][j])/2.0;
			dT[NBr1+1+i][j] = (dW[NBr1+1+i][j]+dE[NBr1+1+i][j])/2.0;
		}
		dN[NBr1+1+NBr2][j] = (rB4-rB3)/NBr3;
		dS[NBr1+1+NBr2][j] = (rB3-rB2)/NBr2;
		dW[NBr1+1+NBr2][j] = 2*Pi/NNBt;
		dE[NBr1+1+NBr2][j] = 2*Pi/NNBt;
		dR[NBr1+1+NBr2][j] = (dN[NBr1+1+NBr2][j]+dS[NBr1+1+NBr2][j])/2.0;
		dT[NBr1+1+NBr2][j] = (dW[NBr1+1+NBr2][j]+dE[NBr1+1+NBr2][j])/2.0;
	}
	// dr3
	for(j=1; j<=NNBt; j++){
		for(i=1; i<NBr3; i++){
			dN[NBr1+1+NBr2+i][j] = (rB4-rB3)/NBr3;
			dS[NBr1+1+NBr2+i][j] = (rB4-rB3)/NBr3;
			dW[NBr1+1+NBr2+i][j] = 2*Pi/NNBt;
			dE[NBr1+1+NBr2+i][j] = 2*Pi/NNBt;
			dR[NBr1+1+NBr2+i][j] = (dN[NBr1+1+NBr2+i][j]+dS[NBr1+1+NBr2+i][j])/2.0;
			dT[NBr1+1+NBr2+i][j] = (dW[NBr1+1+NBr2+i][j]+dE[NBr1+1+NBr2+i][j])/2.0;
		}
		dN[NBr1+1+NBr2+NBr3][j] = (rB5-rB4)/NBr4;
		dS[NBr1+1+NBr2+NBr3][j] = (rB4-rB3)/NBr3;
		dW[NBr1+1+NBr2+NBr3][j] = 2*Pi/NNBt;
		dE[NBr1+1+NBr2+NBr3][j] = 2*Pi/NNBt;
		dR[NBr1+1+NBr2+NBr3][j] = (dN[NBr1+1+NBr2+NBr3][j]+dS[NBr1+1+NBr2+NBr3][j])/2.0;
		dT[NBr1+1+NBr2+NBr3][j] = (dW[NBr1+1+NBr2+NBr3][j]+dE[NBr1+1+NBr2+NBr3][j])/2.0;
	}
	// dr4
	for(j=1; j<=NNBt; j++){
		for(i=1; i<NBr4; i++){
			dN[NBr1+1+NBr2+NBr3+i][j] = (rB5-rB4)/NBr4;
			dS[NBr1+1+NBr2+NBr3+i][j] = (rB5-rB4)/NBr4;
			dW[NBr1+1+NBr2+NBr3+i][j] = 2*Pi/NNBt;
			dE[NBr1+1+NBr2+NBr3+i][j] = 2*Pi/NNBt;
			dR[NBr1+1+NBr2+NBr3+i][j] = (dN[NBr1+1+NBr2+NBr3+i][j]+dS[NBr1+1+NBr2+NBr3+i][j])/2.0;
			dT[NBr1+1+NBr2+NBr3+i][j] = (dW[NBr1+1+NBr2+NBr3+i][j]+dE[NBr1+1+NBr2+NBr3+i][j])/2.0;
		}
		dN[NBr1+1+NBr2+NBr3+NBr4][j] = (rB6-rB5)/NBr5;
		dS[NBr1+1+NBr2+NBr3+NBr4][j] = (rB5-rB4)/NBr4;
		dW[NBr1+1+NBr2+NBr3+NBr4][j] = 2*Pi/NNBt;
		dE[NBr1+1+NBr2+NBr3+NBr4][j] = 2*Pi/NNBt;
		dR[NBr1+1+NBr2+NBr3+NBr4][j] = (dN[NBr1+1+NBr2+NBr3+NBr4][j]+dS[NBr1+1+NBr2+NBr3+NBr4][j])/2.0;
		dT[NBr1+1+NBr2+NBr3+NBr4][j] = (dW[NBr1+1+NBr2+NBr3+NBr4][j]+dE[NBr1+1+NBr2+NBr3+NBr4][j])/2.0;
	}
	// dr5
	for(j=1; j<=NNBt; j++){
		for(i=1; i<NBr5; i++){
			dN[NBr1+1+NBr2+NBr3+NBr4+i][j] = (rB6-rB5)/NBr5;
			dS[NBr1+1+NBr2+NBr3+NBr4+i][j] = (rB6-rB5)/NBr5;
			dW[NBr1+1+NBr2+NBr3+NBr4+i][j] = 2*Pi/NNBt;
			dE[NBr1+1+NBr2+NBr3+NBr4+i][j] = 2*Pi/NNBt;
			dR[NBr1+1+NBr2+NBr3+NBr4+i][j] = (dN[NBr1+1+NBr2+NBr3+NBr4+i][j]+dS[NBr1+1+NBr2+NBr3+NBr4+i][j])/2.0;
			dT[NBr1+1+NBr2+NBr3+NBr4+i][j] = (dW[NBr1+1+NBr2+NBr3+NBr4+i][j]+dE[NBr1+1+NBr2+NBr3+NBr4+i][j])/2.0;
		}
		dN[NNBr][j] = TINY;
		dS[NNBr][j] = (rB6-rB5)/NBr5;
		dW[NNBr][j] = 2*Pi/NNBt;
		dE[NNBr][j] = 2*Pi/NNBt;
		dR[NNBr][j] = (dN[NNBr][j]+dS[NNBr][j])/2.0;
		dT[NNBr][j] = (dW[NNBr][j]+dE[NNBr][j])/2.0;
	}

	// Radius of each node
	for(j=1; j<=NNBt; j++){
		Rvb[1][j] = rB1;
		Rw[1][j] = Rvb[1][j];            
		Re[1][j] = Rvb[1][j];           
		Rn[1][j] = Rvb[1][j] + dN[1][j]/2.0; 
		Rs[1][j] = Rvb[1][j] - dS[1][j]/2.0; 
	}
	for(j=1; j<=NNBt; j++){
		for(i=2; i<=NNBr; i++){
			Rvb[i][j] = Rvb[i-1][j] + dS[i][j];	 // Radius of node[i][j]
			Rw[i][j] = Rvb[i][j];            
			Re[i][j] = Rvb[i][j];           
			Rn[i][j] = Rvb[i][j] + dN[i][j]/2.0; 
			Rs[i][j] = Rvb[i][j] - dS[i][j]/2.0; 
		}
	}
	// angle of each node
	for(i=1; i<=NNBr; i++){
		T[i][1] = 0.0;
		Tw[i][1] = T[i][j] - dW[i][j]/2.0;
		Te[i][1] = T[i][j] + dE[i][j]/2.0;
		Tn[i][1] = T[i][j]; 
		Ts[i][1] = T[i][j]; 
	}
	for(j=2; j<=NNBt; j++){
		for(i=1; i<=NNBr; i++){
			T[i][j] = T[i][j-1] + dW[i][j];	// Angle of node[i][j]
			Tw[i][j] = T[i][j] - dW[i][j]/2.0;
			Te[i][j] = T[i][j] + dE[i][j]/2.0;
			Tn[i][j] = T[i][j]; 
			Ts[i][j] = T[i][j]; 
		}
	}
	//printf("%e \n", Rvb[NNBr][1]);
}
/*	Area of each control volume */
void area(float **Acv, float *Acvt, float **Rn, float **Rs, float **dT)
{
	int i, j, k, NBKt, dNBKt;

	//	Area of each control volume
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBr; i++){
			Acv[i][j] = 0.5*(Rn[i][j]*Rn[i][j]-Rs[i][j]*Rs[i][j])*dT[i][j]; 
		}
		//Acv[1][j] = 0.5*(Rn[i][j]*Rn[i][j]-Rs[i][j]*Rs[i][j])*dT[i][j]*0.5;   // half-CV
		//Acv[NNBr][j] = 0.5*(Rn[i][j]*Rn[i][j]-Rs[i][j]*Rs[i][j])*dT[i][j]*0.5;  // half-CV
	}	

	// Area in Barrel Kidney
	NBKt = int(thetaK/dNBt+1.5);		// Barrel Kidney 각도방향 grid 수	
	dNBKt = int(NNBt/N+0.5);			// Barrel Kidney 사이의 각도방향 grid 간격
	//printf("%d %d \n", dNBKt, NBKt);
	for(i=1; i<N; i++){
		for(k=NBr1+2; k<=NBr1+NBr2; k++){
			for(j=1+dNBKt*(i-1); j<=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j++){			
				Acv[k][j] = 0.0;
			}			
			for(j=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j<=dNBKt*i; j++){
				Acv[k][j] = 0.0;
			}
		}
		for(j=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(i-1); j++){
			for(k=2+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); 
				k<=NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); k++){
					Acv[k][j] =0.0;
			}
		}
		for(j=1+dNBKt*i-NBKt/2; j<=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j++){
			for(k=2+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo); 
				k<=NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo); k++){
					Acv[k][j] = 0.0;
			}
		}
	}
		// Barrel rounding
	for(k=NBr1+2; k<=NBr1+NBr2; k++){
		for(j=1+dNBKt*(N-1); j<=1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt); j++){			
			Acv[k][j] = 0.0;
		}
		for(j=1+dNBKt*N-NBKt/2+int(rB/R/dNBt); j<=dNBKt*N; j++){
			Acv[k][j] = 0.0;
		}
	}
	for(j=1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(N-1); j++){
		for(k=2+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); 
			k<=NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); k++){
				Acv[k][j] = 0.0;
		}
	}
	for(j=1+dNBKt*N-NBKt/2; j<=1+dNBKt*N-NBKt/2+int(rB/R/dNBt); j++){
		for(k=2+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*N-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo); 
			k<=NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*N-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo); k++){
				Acv[k][j] = 0.0;
		}
	}

	//	Total Nondimentsional face area
	*Acvt = 0.0;		
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBr; i++){
			*Acvt = *Acvt + Acv[i][j];
		}
	}
	//printf("%e \n", *Acvt);
}

/* Initial Condition & Configuration */
void initial(float **Pvb, double *ym, double **Hvb, double **Hvb_gro, float **Rvb, float **Tvb, double *Pc, float *Baleakage)
{
	void groove(double **H_gro, float **Tvb);

	int i,j, k;

	V0 = Vpi + (Ld+2*R*tan(beta))*Ap;	// Initial cylinder volume

	// Initial cilinder pressure
	Pc[1] = 0.0e5;
	Pc[2] = Pout*1.;
	Pc[3] = Pout*1.;
	Pc[4] = Pout*1.;
	Pc[5] = Pout*1.;
	Pc[6] = Pin*0;
	Pc[7] = Pin*1;
	Pc[8] = Pin*1;
	Pc[9] = Pin*0;
	Pc[10] = Pout*1.0;
	Pc[11] = Pin*1.0;
	Pc[12] = 0.0000000015/*0.1*w*R*(tan(beta)*sin(2*Pi/9)+tan(beta2)/cos(beta)*cos(2*Pi/9*2))*Ap*/;
	//Pc[13] = Pout*1.0;

	// Initial Pressure Distribution
	for(j=1; j<NNBt/2; j++){
		for(i=2; i<NNBr; i++){
			Pvb[i][j] = Pout/2/Lambda;
		}
	}
	for(j=NNBt-49; j<=NNBt; j++){
		for(i=2; i<NNBr; i++){
			Pvb[i][j] = Pout/2/Lambda;
		}
	}
	for(j=NNBt/2; j<=NNBt-50; j++){
		for(i=2; i<NNBr; i++){
			Pvb[i][j] = Pin/Lambda;
		}
	}
	// Inner & Outer boundary condition
	for(j=1; j<=NNBt; j++){
		Pvb[NNBr][j] = Ph/Lambda;
		Pvb[1][j] = Ph/Lambda;
	}

	double NBKt = int(thetaK/dNBt+1.5);		// Barrel Kidney 각도방향 grid 수	
	double dNBKt = int(NNBt/N+0.5);			// Barrel Kidney 사이의 각도방향 grid 간격

	for(i=1; i<N; i++){
		for(k=NBr1+1; k<=NBr1+1+NBr2; k++){
			for(j=1+dNBKt*(i-1); j<=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j++){			
				Pvb[k][j] = Pc[i]/Lambda;
			}			
			for(j=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j<=dNBKt*i; j++){
				Pvb[k][j] = Pc[i+1]/Lambda;
			}
		}
		for(j=1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(i-1); j++){
			for(k=1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); 
				k<=1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(i-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); k++){
					Pvb[k][j] = Pc[i+1]/Lambda;
			}
		}
		for(j=1+dNBKt*i-NBKt/2; j<=1+dNBKt*i-NBKt/2+int(rB/R/dNBt); j++){
			for(k=1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo); 
				k<=1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*i-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo); k++){
					Pvb[k][j] = Pc[i]/Lambda;
			}
		}
	}
		// Barrel rounding
	for(k=NBr1+1; k<=NBr1+1+NBr2; k++){
		for(j=1+dNBKt*(N-1); j<=1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt); j++){			
			Pvb[k][j] = Pc[N]/Lambda;
		}
		for(j=1+dNBKt*N-NBKt/2+int(rB/R/dNBt); j<=dNBKt*N; j++){
			Pvb[k][j] = Pc[1]/Lambda;
		}
	}
	for(j=1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt); j<=1+NBKt/2+dNBKt*(N-1); j++){
		for(k=1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); 
			k<=1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+NBKt/2+dNBKt*(N-1)-int(rB/R/dNBt))),2))*NBr2/(rB3-rB2)/Rbo); k++){
				Pvb[k][j] = Pc[N]/Lambda;
		}
	}
	for(j=1+dNBKt*N-NBKt/2; j<=1+dNBKt*N-NBKt/2+int(rB/R/dNBt); j++){
		for(k=1+NBr1+NBr2/2-int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*N-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo); 
			k<=1+NBr1+NBr2/2+int(sqrt(rB*rB-pow(R*dNBt*(j-(1+dNBKt*N-NBKt/2))-rB,2))*NBr2/(rB3-rB2)/Rbo); k++){
				Pvb[k][j] = Pc[1]/Lambda;
		}
	}

	// Initial motion variables
	ym[1] = 1.0;
	ym[2] = theta0;
	ym[3] = psi0;
	ym[4] = dhv0/w;
	ym[5] = wx0/w;
	ym[6] = wy0/w;
	double gamma0 = acos(cos(ym[2])*cos(ym[3]))*Rbo/hv0;
	double zeta0 = -atan2(tan(ym[2]),sin(ym[3]));

	// Initial film thickness
	groove(Hvb_gro, Tvb);
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBr; i++){
			Hvb[i][j] = ym[1] - Rvb[i][j]*gamma0*cos(Tvb[i][j]-(Pi/2-zeta0)) + Hvb_gro[i][j];			
		}
	}	

	// Initialize leakages
	for(i=1; i<=N;i++){		
		Baleakage[i] = 0.0;
	}
}
/* Ring Groove depth */
void groove(double **H_gro, float **Tvb)
{
	int i,j,k,l;
	int NBch, dNBch;

	double **H_g;

	H_g = dmatrix(1,NNBr,1,NNBt);

	// Initialize groove depth
	for(j=1; j<=NNBt; j++){
		for(i=1; i<=NNBr; i++){
			H_gro[i][j] = 0.0;
		}
	}
	// groove ring
	for(j=1; j<=NNBt; j++){
		for(i=3+NBr1+NBr2+NBr3; i<NBr1+NBr2+NBr3+NBr4; i++){
			H_gro[i][j] = Hg;
		}
		H_gro[2+NBr1+NBr2+NBr3][j] = Hg/2.0;
		H_gro[NBr1+NBr2+NBr3+NBr4][j] = Hg/2.0;
	}
	// Nch Channels into the housing 
	NBch = int(thetach/(2*Pi/NNBt)+0.5);
	dNBch = int(NNBt/Nch+0.5);
	for(i=1; i<=Nch; i++){
		for(k=NBr1+NBr2+NBr3+NBr4; k<=NNBr; k++){
			for(j=NNBt/12+1+dNBch*(i-1); j<=NNBt/12-1+NBch/2+dNBch*(i-1); j++){			
				H_gro[k][j] = Hg;
			}
			H_gro[k][NNBt/12+NBch/2+dNBch*(i-1)] = Hg/2.0;
			for(j=NNBt/12+3+dNBch*(i-1)-NBch/2; j<=NNBt/12+dNBch*(i-1); j++){				
				H_gro[k][j] = Hg;
			}
			H_gro[k][NNBt/12+2+dNBch*(i-1)-NBch/2] = Hg/2.0;
		}
	}

	// pocket groove
	//int Ngr = 24;
	//float Hgr = 6.0e-6/hv0;
	//float thetagr = 5.0*Pi/180;
	//float drG = 1.6e-3/Rbo;
	//int NBgrdR = int(drG/(rB6-rB5)*NBr5+0.5);
	//int NBgrR1 = int((NBr5-NBgrdR)/2+0.5);
	//int NBgrR2 = NBgrR1 + NBgrdR;
	//int NBgrT = int(thetagr/(2*Pi/NNBt)+0.5);
	//int dNBgrT = int(NNBt/Ngr+0.5);
	//for(i=1; i<=Ngr; i++){
	//	for(k=1+NBr1+NBr2+NBr3+NBr4+NBgrR1+2; k<=1+NBr1+NBr2+NBr3+NBr4+NBgrR2-2; k++){
	//		for(j=1+dNBgrT*(i-1)+2; j<=1+NBgrT+dNBgrT*(i-1)-2; j++){
	//			H_gro[k][j] = Hgr;
	//		}
	//	}
	//	for(k=1+NBr1+NBr2+NBr3+NBr4+NBgrR1+1; k<=1+NBr1+NBr2+NBr3+NBr4+NBgrR2-1; k++){
	//		H_gro[k][1+dNBgrT*(i-1)+1] = Hgr/2;
	//		H_gro[k][1+dNBgrT*i-1] = Hgr/2;
	//	}
	//	for(j=1+dNBgrT*(i-1)+1; j<=1+NBgrT+dNBgrT*(i-1)-1; j++){
	//		H_gro[1+NBr1+NBr2+NBr3+NBr4+NBgrR1+1][j] = Hgr/2;
	//		H_gro[1+NBr1+NBr2+NBr3+NBr4+NBgrR2-1][j] = Hgr/2;
	//	}					
	//}


	// Wavy groove - Circumerential
	//float Wh1 = 0.6e-6/hv0, Wh2 = 0.8e-6/hv0, Wh3 = 0.8e-6/hv0, Wh5 = 0.4e-6/hv0;
	//int Nwa1 = 30, Nwa2 = 38, Nwa3 = 46, Nwa5 = 54;
	//float thetawa1 = 6.0*Pi/180, thetawa2 = 5.0*Pi/180, thetawa3 = 4.0*Pi/180, thetawa5 = 3.0*Pi/180;
	//int NBwaT1 = int(thetawa1/(2*Pi/NNBt)+0.5);
	//int NBwaT2 = int(thetawa2/(2*Pi/NNBt)+0.5);
	//int NBwaT3 = int(thetawa3/(2*Pi/NNBt)+0.5);
	//int NBwaT5 = int(thetawa5/(2*Pi/NNBt)+0.5);
	//int dNBwaT1 = int(NNBt/Nwa1+0.5);
	//int dNBwaT2 = int(NNBt/Nwa2+0.5);
	//int dNBwaT3 = int(NNBt/Nwa3+0.5);
	//int dNBwaT5 = int(NNBt/Nwa5+0.5);

	//for(k=1; k<=Nwa1; k++){
	//	for(i=1; i<=1+NBr1; i++){
	//		for(j=1+dNBwaT1*(k-1); j<=1+NBwaT1+dNBwaT1*(k-1); j++){
	//			H_g[i][j] = Wh1*sin((j-dNBwaT1*(k-1)-1)*dNBt/thetawa1*2*Pi);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//}
	//for(k=1; k<=Nwa2; k++){
	//	for(i=2+NBr1; i<=1+NBr1+NBr2; i++){
	//		for(j=1+dNBwaT2*(k-1); j<=1+NBwaT2+dNBwaT2*(k-1); j++){
	//			H_g[i][j] = Wh2*sin((j-dNBwaT2*(k-1)-1)*dNBt/thetawa2*2*Pi);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//}
	//for(k=1; k<=Nwa3; k++){
	//	for(i=2+NBr1+NBr2; i<=1+NBr1+NBr2+NBr3; i++){
	//		for(j=1+dNBwaT3*(k-1); j<=1+NBwaT3+dNBwaT3*(k-1); j++){
	//			H_g[i][j] = Wh3*sin((j-dNBwaT3*(k-1)-1)*dNBt/thetawa3*2*Pi);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//}
	//for(k=1; k<=Nwa5; k++){
	//	for(i=2+NBr1+NBr2+NBr3+NBr4; i<=NNBr; i++){
	//		for(j=1+dNBwaT5*(k-1); j<=1+NBwaT5+dNBwaT5*(k-1); j++){
	//			H_g[i][j] = Wh5*sin((j-dNBwaT5*(k-1)-1)*dNBt/thetawa5*2*Pi);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//}

	// Wavy groove - Radial
	//float WhR1 = 1.6e-6/hv0, WhR2 = 0.4e-6/hv0, WhR3 = 2.0e-6/hv0, WhR5 = 0.6e-6/hv0;
	//int NwaR1 = 1, NwaR2 = 2, NwaR3 = 1, NwaR5 = 2;
	//float drwa1=2*(rB2-rB1), drwa2=0.5*(rB3-rB2), drwa3=2*(rB4-rB3), drwa5=0.5*(rB5-rB4);
	//int NBwaR1 = int(drwa1/((rB2-rB1)/NBr1)+0.5);
	//int NBwaR2 = int(drwa2/((rB3-rB2)/NBr2)+0.5);
	//int NBwaR3 = int(drwa3/((rB4-rB3)/NBr3)+0.5);
	//int NBwaR5 = int(drwa5/((rB6-rB5)/NBr5)+0.5);
	//int dNBwaR1 = int(NBr1/float(NwaR1)+0.5);
	//int dNBwaR2 = int(NBr2/float(NwaR2)+0.5);
	//int dNBwaR3 = int(NBr3/float(NwaR3)+0.5);
	//int dNBwaR5 = int(NBr5/float(NwaR5)+0.5);

	//for(j=1; j<=NNBt; j++){
	//	for(k=1; k<=NwaR1; k++){
	//		for(i=1+dNBwaR1*(k-1); i<=1+NBwaR1+dNBwaR1*(k-1); i++){
	//			H_g[i][j] = WhR1*sin(float(i-1)/NBwaR1*2.0*Pi*NwaR1);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//	for(k=1; k<=NwaR2; k++){
	//		for(i=1+NBr1+dNBwaR2*(k-1); i<=1+NBr1+NBwaR2+dNBwaR2*(k-1); i++){
	//			H_g[i][j] = WhR2*sin(float(i-(1+NBr1))/NBwaR2*2.0*Pi*NwaR2);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//	for(k=1; k<=NwaR3; k++){
	//		for(i=1+NBr1+NBr2+dNBwaR3*(k-1); i<=1+NBr1+NBr2+NBwaR3+dNBwaR3*(k-1); i++){
	//			H_g[i][j] = WhR3*sin(float(i-(1+NBr1+NBr2))/NBwaR3*2.0*Pi*NwaR3);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//	for(k=1; k<=NwaR5; k++){
	//		for(i=1+NBr1+NBr2+NBr3+NBr4+dNBwaR5*(k-1); i<=1+NBr1+NBr2+NBr3+NBr4+NBwaR5+dNBwaR5*(k-1); i++){
	//			H_g[i][j] = WhR5*sin(float(i-(1+NBr1+NBr2+NBr3+NBr4))/NBwaR5*2.0*Pi*NwaR5);
	//			H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//		}
	//	}
	//}

	// Wavy groove - 3D type 4
	float Wh1 = -0.4e-6/hv0, Wh2 = -0.4e-6/hv0, Wh3 = -0.4e-6/hv0, Wh5 = -4.0e-6/hv0;
	int NwaC1 = 108, NwaC2 = 108, NwaC3 = 108, NwaC5 = 9;
	float thetawa1 = 3.33333*Pi/180, thetawa2 = 3.33333*Pi/180, thetawa3 = 3.33333*Pi/180, thetawa5 = 40*Pi/180;
	int NBwaT1 = int(thetawa1/(2*Pi/NNBt)+0.5);
	int NBwaT2 = int(thetawa2/(2*Pi/NNBt)+0.5);
	int NBwaT3 = int(thetawa3/(2*Pi/NNBt)+0.5);
	int NBwaT5 = int(thetawa5/(2*Pi/NNBt)+0.5);
	int dNBwaT1 = int(NNBt/NwaC1+0.5);
	int dNBwaT2 = int(NNBt/NwaC2+0.5);
	int dNBwaT3 = int(NNBt/NwaC3+0.5);
	int dNBwaT5 = int(NNBt/NwaC5+0.5);
	//float WhR1 = -0.8e-6/hv0, WhR2 = -0.8e-6/hv0, WhR3 = -0.8e-6/hv0, WhR5 = -0.8e-6/hv0;
	int NwaR1 = 2, NwaR2 = 2, NwaR3 = 2, NwaR5 = 2;
	float drwa1=(rB2-rB1)/NwaR1, drwa2=(rB3-rB2)/NwaR2, drwa3=(rB4-rB3)/NwaR3, drwa5=(rB6-rB5)/NwaR5;
	int NBwaR1 = int(drwa1/((rB2-rB1)/NBr1)+0.5);
	int NBwaR2 = int(drwa2/((rB3-rB2)/NBr2)+0.5);
	int NBwaR3 = int(drwa3/((rB4-rB3)/NBr3)+0.5);
	int NBwaR5 = int(drwa5/((rB6-rB5)/NBr5)+0.5);
	int dNBwaR1 = int(NBr1/float(NwaR1)+0.5);
	int dNBwaR2 = int(NBr2/float(NwaR2)+0.5);
	int dNBwaR3 = int(NBr3/float(NwaR3)+0.5);
	int dNBwaR5 = int(NBr5/float(NwaR5)+0.5);

	//for(l=1; l<=NwaR1; l++){
	//	for(k=1; k<=NwaC1; k++){
	//		for(i=1+dNBwaR1*(l-1); i<=1+NBwaR1+dNBwaR1*(l-1); i++){
	//			for(j=1+dNBwaT1*(k-1); j<=1+NBwaT1+dNBwaT1*(k-1); j++){
	//				H_g[i][j] = Wh1/4*( (1+cos(float(i-1-dNBwaR1*(l-1))/NBwaR1*2.0*Pi)) * (1-cos((j-1-dNBwaT1*(k-1))*dNBt/thetawa1*2.0*Pi)) - 1 );
	//				H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//			}
	//		}
	//	}
	//}
	//for(l=1; l<=NwaR2; l++){
	//	for(k=1; k<=NwaC2; k++){
	//		for(i=1+NBr1+dNBwaR2*(l-1); i<=1+NBr1+NBwaR2+dNBwaR2*(l-1); i++){
	//			for(j=1+dNBwaT2*(k-1); j<=1+NBwaT2+dNBwaT2*(k-1); j++){
	//				H_g[i][j] = Wh2/4*( (1+cos(float(i-1-NBr1-dNBwaR2*(l-1))/NBwaR2*2.0*Pi)) * (1-cos((j-1-dNBwaT2*(k-1))*dNBt/thetawa2*2.0*Pi)) - 1 );
	//				H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//			}
	//		}
	//	}
	//}
	//for(l=1; l<=NwaR3; l++){
	//	for(k=1; k<=NwaC3; k++){
	//		for(i=1+NBr1+NBr2+dNBwaR3*(l-1); i<=1+NBr1+NBr2+NBwaR3+dNBwaR3*(l-1); i++){
	//			for(j=1+dNBwaT3*(k-1); j<=1+NBwaT3+dNBwaT3*(k-1); j++){
	//				H_g[i][j] = Wh3/4*( (1+cos(float(i-1-NBr1-NBr2-dNBwaR3*(l-1))/NBwaR3*2.0*Pi)) * (1-cos((j-1-dNBwaT3*(k-1))*dNBt/thetawa3*2.0*Pi)) - 1 );
	//				H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
	//			}
	//		}
	//	}
	//}

		for(k=1; k<=NwaC5; k++){
			for(i=1+NBr1+NBr2+NBr3+NBr4; i<=1+NBr1+NBr2+NBr3+NBr4+NBr5; i++){
				for(j=1+dNBwaT5*(k-1); j<1+NBwaT5+dNBwaT5*(k-1); j++){
					H_g[i][j] = Wh5/4*( 2 * (1-cos((j-1-dNBwaT5*(k-1))*dNBt/thetawa5*2.0*Pi)) - 1 );
					H_gro[i][j] = H_gro[i][j] + H_g[i][j]; 
				}
			}
		}

	free_dmatrix(H_g,1,NNBr,1,NNBt);

}


