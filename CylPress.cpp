#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

extern double Avo(double );

extern int Nd;
extern int Nphi, NBz1, NBz2, NBz3, NBz4, NBz5, NBz6, NBz, NNBz, NNBt, Nch;
extern double dNBt, thetach;
extern double R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, revf, phif, phiout_c, phiout_vb, phiout_FM, phiout_LK, w;
extern double V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
extern double integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1, zB1, zB2, zB3, zB4, zB5, zB6, Fs, e0, theta0, psi0, de0, wx0, wy0, Hg, Mp, Mb, Iba, Ibt;
extern double Lambda, Loc, Lpj, Lpg, Lbj, Rbo;
//extern float FxD, MyxD, MyyD, FyD, MxxD, MxyD, FlatxD, FlatyD, MlatxD, MlatyD, FfricD, MfricxD, MfricyD, frD, MxD, MyD, MzD, leakageRinD, leakageRoutD;

/* R-K 4th for cylinder pressure */
void RK4c(double *P, double *dP, int n, double *phiR, double htry, double *Pscal, double *hdid, double *hnext, void (*Derivs)(double, double *, double *, float *), float *PCleakage)
{
	int i, jtry;	
	double *K1, *K2, *K3, *K4, *K5, *K6, phisav, *Pt, *Psav;	
	double *err, errmax, h;

	K1 = dvector(1,n);
	K2 = dvector(1,n);
	K3 = dvector(1,n);
	K4 = dvector(1,n);
	K5 = dvector(1,n);
	K6 = dvector(1,n);
	err = dvector(1,n);
	Pt = dvector(1,n);
	Psav = dvector(1,n);
	phisav = (*phiR);
	//save initial values	
	for(i=1; i<=n; i++){
		Psav[i] = P[i];
		K1[i] = dP[i];
	}
	//set step size to the initial trial value
	h = htry;
	for(jtry=1; jtry<=MAXTRY; jtry++){
		for(i=1; i<=n; i++)	Pt[i] = Psav[i] + A21*K1[i]*h;	
		(*Derivs)(phisav+A2X*h, Pt, K2, PCleakage);
		for(i=1; i<=n; i++)	Pt[i] = Psav[i] + (A31*K1[i]+A32*K2[i])*h;
		(*Derivs)(phisav+A3X*h, Pt, K3, PCleakage);
		for(i=1; i<=n; i++)	Pt[i] = Psav[i] + (A41*K1[i]+A42*K2[i]+A43*K3[i])*h;	
		(*Derivs)(phisav+A4X*h, Pt, K4, PCleakage);
		for(i=1; i<=n; i++)	Pt[i] = Psav[i] + (A51*K1[i]+A52*K2[i]+A53*K3[i]+A54*K4[i])*h;
		(*Derivs)(phisav+A5X*h, Pt, K5, PCleakage);
		for(i=1; i<=n; i++)	Pt[i] = Psav[i] + (A61*K1[i]+A62*K2[i]+A63*K3[i]+A64*K4[i]+A65*K5[i])*h;
		(*Derivs)(phisav+A6X*h, Pt, K6, PCleakage);
		for(i=1; i<=n; i++){
			P[i] = Psav[i] + (B1*K1[i]+B3*K3[i]+B4*K4[i]+B6*K6[i])*h;
			err[i] = (E1*K1[i] + E3*K3[i] + E4*K4[i] + E5*K5[i] + E6*K6[i])*h;			
		}		
		*phiR = phisav + h;
		if(*phiR == phisav/2) nrerror("stepsize not significant in RK4c");
		errmax = 0.0;			//Evaluate accuracy		
		for(i=1; i<=n; i++) errmax = FMAX(errmax, fabs(err[i]/Pscal[i]));
		//for(i=1; i<=n; i++) printf("%e \n", h);
		errmax /= eps;			//scale Zelative to the ZequiZed tolerance 		
		if(errmax <= 1.0){
			*hdid = h;
			*hnext = (errmax > ERRCON ? SAFETY*h*pow(errmax,PGROW) : GROW*h);
			if(*hnext>hmax*Pi/180) *hnext = hmax*Pi/180;
			//printf("%lf %e \n",*phiR,*hnext);
			free_dvector(Pt,1,n);
			free_dvector(Psav,1,n);
			free_dvector(err,1,n);
			free_dvector(K1,1,n);
			free_dvector(K2,1,n);
			free_dvector(K3,1,n);
			free_dvector(K4,1,n);
			free_dvector(K5,1,n);
			free_dvector(K6,1,n);
			return;
		} 
		else{				//Truncation error too large, Zeduce stepsize
			*hnext = SAFETY*h*pow(errmax,PSHZnK);
			h = (h >= 0.0 ? FMAX(*hnext,SHZnK*h) : FMIN(*hnext,SHZnK*h));
		}
	}
	nrerror("exceeded MAXTRY in RKstep");
}
/* dP/dphi */
void Derivs(double phid, double *Pt, double *dP, float *PCleakage)
{
	int i, sgn1, sgn2;
	double Cv, *Pr, Pv, *phiD, *Lpc, Ao1, Ao2, Qvo, Qvi, Lvp=10e-3, Ain, Aout, bv1, bv2, bvo=10e-3, Lv1=10e-3, Rpipei=10e-3, Rpipe=5e-3, Lv2=20e-3, Lv3=150e-3, Lv4=10e-3, tanv, Bv;
	float *Qcb;

	phiD = dvector(1,N);
	Pr = dvector(1,N);
	Lpc = dvector(1,N);
	Qcb = vector(1,N);
	
	while(phid >= 2*Pi) phid = phid - 2*Pi;

	for(i=1; i<=N; i++){
		phiD[i] = phid + 2*Pi/N*(i-1);
		if(phiD[i] > 2*Pi) phiD[i] = phiD[i] - 2*Pi;	
		//hv = hv - 40e-6*sin(phiD[i]);
		Lpc[i] = Lc - Ld - R*(tan(beta)*(1+cos(phiD[i]))-tan(beta2)/cos(beta)*sin(phiD[i]));
		if(Pt[i]>Ph) Pr[i] = Ph + 2*(Pt[i]-Ph)/(1+sqrt(1+pow(Pi*pow(hs,3)/6/eta/log(Rs2/Rs1),2)*2*rho*(Pt[i]-Ph)/pow(Cs,2)/pow(Pi*Rso*Rso,2)));
		else Pr[i] = Ph + 2*(Pt[i]-Ph)/(1+sqrt(1+pow(Pi*pow(hs,3)/6/eta/log(Rs2/Rs1),2)*2*rho*(Ph-Pt[i])/pow(Cs,2)/pow(Pi*Rso*Rso,2)));
		//theta12(phiD[i]);
		Cv = (phiD[i]<=Pi?Cvo:Cvi);
		sgn1 = (Pt[i]>=Pv?1.1:-1);
		sgn2 = (Pt[i]>=Pr[i]?1:-1);
		if(phiD[i]<=Pi) Pv = Pt[N+1];
		else Pv = Pt[N+2];
		if(Nd >= 1){
			if(i == 1) Qcb[i] = PCleakage[1];
			else Qcb[i] = 2*Pi*Rp*pow(hp,3)/12/eta/Lpc[i]*(Pt[i]-Ph)-2*Pi*Rp/2*hp*w*R*(tan(beta)*sin(phiD[i])+tan(beta2)/cos(beta)*cos(phiD[i]));			
		}
		else Qcb[i] = 2*Pi*Rp*pow(hp,3)/12/eta/Lpc[i]*(Pt[i]-Ph)-2*Pi*Rp/2*hp*w*R*(tan(beta)*sin(phiD[i])+tan(beta2)/cos(beta)*cos(phiD[i]));		
		
		dP[i] = -K/w/(V0-R*(tan(beta)*(1-cos(phiD[i]))+tan(beta2)/cos(beta)*sin(phiD[i]))*Ap)
				*( -w*R*(tan(beta)*sin(phiD[i])+tan(beta2)/cos(beta)*cos(phiD[i]))*Ap + Cv*Avo(phiD[i])*sqrt(2/rho)*sqrt(fabs(Pt[i]-Pv))*sgn1
				+ Qcb[i]/*2*Pi*Rp*pow(hp,3)/12/eta/Lpc[i]*(Pt[i]-Ph)-2*Pi*Rp/2*hp*w*R*(tan(beta)*sin(phiD[i])+tan(beta2)/cos(beta)*cos(phiD[i]))*/
				+ /*Qcb[i]*/pow(hv,3)/12/eta*(Pt[i]-Ph)*(1/log(1+LB1/(R-rB))+1/log(1+LB2/(R+rB)))*2*(thetaB)  
				+ Cs*Pi*Rso*Rso*sqrt(2/rho)*sqrt(fabs(Pt[i]-Pr[i]))*sgn2 );
		//printf("%f  ",dP[i]);
	}
	Qvo = 0., Qvi = 0.;
	Ao1 = 0., Ao2 = 0.;
	for(i=1; i<=N; i++){
		if(phiD[i]<=Pi){
			if(Pt[i]>Pt[N+1]) Qvo = Qvo + Cvo*Avo(phiD[i])*sqrt(2/rho)*sqrt(Pt[i]-Pt[N+1]);
			else Qvo = Qvo - Cvo*Avo(phiD[i])*sqrt(2/rho)*sqrt(Pt[N+1]-Pt[i]);
			Ao1 = Ao1 + Avo(phiD[i]);
		}
		else{
			if(Pt[i]>Pt[N+2]) Qvi = Qvi + Cvi*Avo(phiD[i])*sqrt(2/rho)*sqrt(Pt[i]-Pt[N+2]);
			else Qvi = Qvi - Cvi*Avo(phiD[i])*sqrt(2/rho)*sqrt(Pt[N+2]-Pt[i]);
			Ao2 = Ao2 + Avo(phiD[i]);
		}		
	}
	Qvo = Qvo - pow(hv,3)/12/eta*(Pt[N+1]-Ph)*(1/log(1+LB1/(R-rV1))+1/log(1+LB2/(R+rV1)))*(2*R*rV1*thetaV2-Ao1)/2/R/rV1;
	Qvi = Qvi - pow(hv,3)/12/eta*(Pt[N+2]-Ph)*(1/log(1+LB1/(R-rV2))+1/log(1+LB2/(R+rV2)))*(2*R*rV2*thetaV5-Ao2)/2/R/rV2;

	//dP[N+3] = 2*Ap/50e-3/w/rho*((Pt[N+1]-Pt[N+4])-Pt[N+3]*Cp/pow(2*Ap,2));  //Qc
	bv1 = R*thetaV2;
	tanv = (bv1-bvo)/2/Lvp;	
	bv2 = bv1 - 2*Lv1*tanv;	
	//Bv = log(bv1/(bv1-2*Lv1*tanv))/2/tanv;
	Aout = Pi*Rpipe*Rpipe;
	Ain = Pi*Rpipei*Rpipei;
	dP[N+1] = 0.0/*K/w/((bv1-Lv1*tanv)*Lv1*2*rV1)*(Qvo-Pt[N+3])*/;  //P1	
	dP[N+2] = 0.0/*K/w/Vvi*(Pt[N+4]-Qvi)*/;  //P2
	//dP[N+3] = Aout/Lv3/w/rho*((Pt[N+1]-Pout*0.92) -64/(rho*2*Rpipe/eta)/2*rho*Pt[N+3]/Aout*Lv3/2/Rpipe -50*rho/2*Pt[N+3]*fabs(Pt[N+3])/pow(Aout*0.5,2) );
	//dP[N+4] = Ain/Lv2/w/rho*((Pin-Pt[N+2]) -64/(rho*2*Rpipei/eta)/2*rho*Pt[N+4]/Ain*Lv2/2/Rpipei -3*rho/2*Pt[N+4]*fabs(Pt[N+4])/pow(Ain,2) );
	//if(Pt[N+5]>Ph) dP[N+5] = K/w/Aout/Lv1*(Pt[N+3]-Cvo*Aout*0.0365*sqrt(2/rho)*sqrt(Pt[N+5]-Ph));
	//else dP[N+5] = K/w/Aout/Lv1*(Pt[N+3]+Cvo*Aout*0.0365*sqrt(2/rho)*sqrt(Ph-Pt[N+5]));
	//if(Pt[N+4]>Pout) dP[N+4] = K/w/Aout/Lv4*(Pt[N+3]-Cvo*Avout*sqrt(2/rho)*sqrt(Pt[N+4]-Pout));  //Pt2 orifce model
	//else dP[N+4] = K/w/Aout/Lv4*(Pt[N+3]+Cvo*Avout*sqrt(2/rho)*sqrt(Pout-Pt[N+4]));	
	//if(Pt[N+2]>Pin) dP[N+2] = K/w/Vvi*(Qvi-Cvin*Avin*sqrt(2/rho)*sqrt(Pt[N+2]-Pin));  
	//else dP[N+2] = K/w/Vvi*(Qvi+Cvin*Avin*sqrt(2/rho)*sqrt(Pin-Pt[N+2]));

	//dP[N+3] = -2*N*(0.1*Pout)*sin(N*phid);

	//printf("%lf %.6e\n", phiD[1]*180/Pi, dP[N+5]);

	free_dvector(phiD,1,N);
	free_dvector(Pr,1,N);
	free_dvector(Lpc,1,N);
	free_vector(Qcb,1,N);
}
/* Leakage */
void Leakage(double phi, double *Pt, double *Qp, double *Qv, double *Qs, double *Qpt, double *Qvt, double *Qst)
{
	double *Lpc, Pv, *Pr;
	double Qvo, Qvi, Ao1, Ao2;	
	double *phiL;

	Pr = dvector(1,N);
	Lpc = dvector(1,N);
	phiL = dvector(1,N);
	//leakage in a cylinder
	Lpc[1] = Lc - Ld - R*(tan(beta)*(1+cos(phi))-tan(beta2)/cos(beta)*sin(phi));
	if(Pt[1]>Ph) Pr[1] = Ph + 2*(Pt[1]-Ph)/(1+sqrt(1+pow(Pi*pow(hs,3)/6/eta/log(Rs2/Rs1),2)*2*rho*(Pt[1]-Ph)/pow(Cs,2)/pow(Pi*Rso*Rso,2)));
	else Pr[1] = Ph + 2*(Pt[1]-Ph)/(1+sqrt(1+pow(Pi*pow(hs,3)/6/eta/log(Rs2/Rs1),2)*2*rho*(Ph-Pt[1])/pow(Cs,2)/pow(Pi*Rso*Rso,2)));
	if(phi<=Pi) Pv = Pt[N+1];
	else Pv = Pt[N+2];
	*Qp = 2*Pi*Rp*pow(hp,3)/12/eta/Lpc[1]*(Pt[1]-Ph)-2*Pi*Rp/2*hp*w*R*(tan(beta)*sin(phi)+tan(beta2)/cos(beta)*cos(phi));
	*Qv = pow(hv,3)/12/eta*(Pt[1]-Ph)*(1/log(1+LB1/(R-rB))+1/log(1+LB2/(R+rB)))*2*thetaB;
	*Qs = Cs*Pi*Rso*Rso*sqrt(2/rho)*sqrt(fabs(Pt[1]-Pr[1]));	
	//total leakage 
	*Qpt = 0., *Qst = 0.;
	Qvo = 0., Qvi = 0.;
	Ao1 = 0., Ao2 = 0.;
	for(int i=1; i<=N; i++){
		phiL[i] = phi + 2*Pi/N*(i-1);
		if(phiL[i] > 2*Pi) phiL[i] = phiL[i] - 2*Pi;
		if(Pt[i]>Ph) Pr[i] = Ph + 2*(Pt[i]-Ph)/(1+sqrt(1+pow(Pi*pow(hs,3)/6/eta/log(Rs2/Rs1),2)*2*rho*(Pt[i]-Ph)/pow(Cs,2)/pow(Pi*Rso*Rso,2)));
		else Pr[i] = Ph + 2*(Pt[i]-Ph)/(1+sqrt(1+pow(Pi*pow(hs,3)/6/eta/log(Rs2/Rs1),2)*2*rho*(Ph-Pt[i])/pow(Cs,2)/pow(Pi*Rso*Rso,2)));
		Lpc[i] = Lc - Ld - R*(tan(beta)*(1+cos(phiL[i]))-tan(beta2)/cos(beta)*sin(phiL[i]));
		*Qpt = *Qpt + 2*Pi*Rp*pow(hp,3)/12/eta/Lpc[i]*(Pt[i]-Ph)-2*Pi*Rp/2*hp*w*R*(tan(beta)*sin(phiL[i])+tan(beta2)/cos(beta)*cos(phiL[i]));
		*Qst = *Qst + Cs*Pi*Rso*Rso*sqrt(2/rho)*sqrt(fabs(Pt[i]-Pr[i]));
		if(phiL[i]<=Pi){
			Qvo = Qvo + pow(hv,3)/12/eta*(Pt[i]-Ph)*(1/log(1+LB1/(R-rB))+1/log(1+LB2/(R+rB)))*2*thetaB;				
			Ao1 = Ao1 + Avo(phiL[i]);
		}
		else{
			Qvi = Qvi + pow(hv,3)/12/eta*(Pt[i]-Ph)*(1/log(1+LB1/(R-rB))+1/log(1+LB2/(R+rB)))*2*thetaB;				
			Ao2 = Ao2 + Avo(phiL[i]);
		}
	}
	Qvo = Qvo - pow(hv,3)/12/eta*(Pt[N+1]-Ph)*(1/log(1+LB1/(R-rV1))+1/log(1+LB2/(R+rV1)))*(2*R*rV1*thetaV2-Ao1)/2/R/rV1;
	Qvi = Qvi - pow(hv,3)/12/eta*(Pt[N+2]-Ph)*(1/log(1+LB1/(R-rV2))+1/log(1+LB2/(R+rV2)))*(2*R*rV2*thetaV5-Ao2)/2/R/rV2;
	*Qvt = Qvo + Qvi;

	free_dvector(Pr,1,N);
	free_dvector(Lpc,1,N);
	free_dvector(phiL,1,N);
}
/* Point of action of the Zesultant force & Moment */
void Point_resul(double phi, double *Pt, float *eBVx, float *eBVy)
{
	float Pe=0.0;	
	double *phiP;
	
	phiP = dvector(1,N);
	*eBVx = 0.0, *eBVy = 0.0;
	for(int i=1; i<=N; i++){
		phiP[i] = phi + 2*Pi/N*(i-1);
		if(phiP[i] > 2*Pi) phiP[i] = phiP[i] - 2*Pi;
		*eBVx = *eBVx + Pt[i]*R*sin(phiP[i]);
		*eBVy = *eBVy + Pt[i]*R*cos(phiP[i]);
		Pe = Pe + Pt[i];
	}
	*eBVx = *eBVx/Pe;
	*eBVy = *eBVy/Pe;

	free_dvector(phiP,1,N);
}