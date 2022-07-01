/* Consider the Momentum Effect in the Valve Port */
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"


double Avo(double );
double Avi(double );
//void theta12(double );
//double integB(double );

int Nphi;
double R, Rp, Rpi, Lpi, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta0, beta2, phif, phiout, w;
double V0, Vpi, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6; 
double integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1;
double g=9.8, Rol=12.7e-3, Rpipe=12.7e-3, Aol, Lol=1000e-3, Lpipe=1000e-3, Apipe;
int PCYCLE=0, Nol=10000, Np=10000;
double c, c2, Ep=200e9, ep=20.5e-3, ep2=20.5e-3, pphi, olphi, phiout_p=30;

int main()
{
	void integBpre();
	void Input();
	void Init(double *, double *, double *, double *, double *);
	void Derivs(double , double *, double *, double *, double *);
	void RK4(double *, double *, double *, double *, int , double *, double , double *, double *, double *, void (*Derivs)(double, double *, double *, double *, double *));
	void Leakage(double , double *, double *, double *, double *, double *, double *, double *);
	void Point_resul(double , double *, double *, double *);
	void STorq(double *, double *, double );
	void Pipe(double , double *, double *, double *, double *, double * );
	void Outlet(double , double *, double *, double *, double *, double * );

	FILE *fp, *gp;

	int i;
	double phi=0.0, phiend;
	double h, hout, hdid, hnext;
	double *P, *dP, *Pscal, Qvo, Qvi, *ph;
	double Qpi, Qlp, Qv, Qs, Qpt, Qvt, Qst, eBVx, eBVy;
	double Ao1, Ao2;
	double *ST;
	double *Pp, *Qp, hpipe, hout_p, houtlet;
	double *Pol, *Qol;

	P = dvector(1,NRK);
	dP = dvector(1,NRK);
	Pscal = dvector(1,NRK);
	ph = dvector(1,N);
	ST = dvector(1,3);
	Pp = dvector(1,Np);
	Qp = dvector(1,Np);
	Pol = dvector(1,Nol);
	Qol = dvector(1,Nol);

	if((fp = fopen("Output_pressure2.dat", "w"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}
	if((gp = fopen("Output_Pipe.dat", "w"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}
	/* Read parameters */
	Input();
	/* Initial condition & Unit conversion*/
	Init(P, Pp, Qp, Pol, Qol);	
	h = 0.000001*Pi/180;

	/* Solve equations & Print the results*/
	//fprintf(fp, "*Input File: Input_pressure2.dat, *Output File: Output_pressure2.dat \n");
	fprintf(fp, "  Phi		P1		P2		P3	     P4		   P5		  P6		P7		P8		P9		Po		Pi \n");
	fprintf(fp, "%4.1f  ", phi*180/Pi);
	for(i=1; i<=NRK; i++) fprintf(fp, "%e  ", P[i]);
	fprintf(fp,"\n");
	hout = 0.0;
	hpipe = 0.0;
	houtlet = 0.0;
	hout_p = 0.0;
	while(phi<phif){		
		printf("%lf \n", phi*180/Pi);
		if(phi>2*Pi){		
			if(hout>=phiout){
				fprintf(fp,"%4.1f  ", (phi-2*Pi)*180/Pi);
				for(i=1; i<=NRK; i++) fprintf(fp,"%e  ", P[i]);
				//Qvo = P[N+3];
				//printf("======================== %lf \n", phi*180/Pi);
				Qvo = 0.0, Qvi = 0.0;		// Outlet flow
				Ao1 = 0., Ao2 = 0.;
				for(i=1; i<=N; i++){
					ph[i] = phi + 2*Pi/N*(i-1);
					while(ph[i] > 2*Pi) ph[i] = ph[i] - 2*Pi;
					if(ph[i]<=Pi){
						if(P[i]>P[N+1]) Qvo = Qvo + Cvo*Avo(ph[i])*sqrt(2/rho)*sqrt(P[i]-P[N+1]);
						else Qvo = Qvo - Cvo*Avo(ph[i])*sqrt(2/rho)*sqrt(P[N+1]-P[i]);
						Ao1 = Ao1 + Avo(ph[i]);
					}
					else{
						if(P[i]>P[N+2]) Qvi = Qvi + Cvi*Avi(ph[i])*sqrt(2/rho)*sqrt(P[i]-P[N+2]);
						else Qvi = Qvi - Cvi*Avi(ph[i])*sqrt(2/rho)*sqrt(P[N+2]-P[i]);
						Ao2 = Ao2 + Avi(ph[i]);
					}		
				}
				Qvo = Qvo - pow(hv,3)/12/eta*(P[N+1]-Ph)*(1/log(1+LB1/(R-rV1))+1/log(1+LB2/(R+rV1)))*(2*R*rV1*thetaV2-Ao1)/2/R/rV1;
				Qvi = Qvi - pow(hv,3)/12/eta*(P[N+2]-Ph)*(1/log(1+LB1/(R-rV2))+1/log(1+LB2/(R+rV2)))*(2*R*rV2*thetaV5-Ao2)/2/R/rV2;


				//fprintf(fp,"%e %e ", Pol[1], Pp[1]);
				fprintf(fp,"%e %e ", -Qvi, Qvo);
				fprintf(fp,"%e %e ", Avo(phi), Avi(phi));
				//Leakage(phi, P, &Qlp, &Qv, &Qs, &Qpt, &Qvt, &Qst);	// Leakage in a Cylinder & Total leakage
				//Point_resul(phi, P, &eBVx, &eBVy);					// Resultant force action point
				//fprintf(fp,"%e %e %e %e %e %e %e ", Qlp, Qv, Qs, Qpt, Qvt, Qst, Qpt+Qvt+Qst);
				//fprintf(fp,"%f %f ", eBVx, eBVy);				
				//fprintf(fp,"%e %e %e ", ST[1], ST[2], ST[3]);
				fprintf(fp,"\n");
				hout = 0.0;
			}	
			//if(hout_p>=phiout_p){
			//	fprintf(gp, "NCYCLE=%d\n VARIABLES=Z,PRESSURE(Pa)\n ZONE T=P_Pipe, I=         %d , F=POINT\n", PCYCLE, Np);
			//	for(i=1; i<=Nol; i++) fprintf(gp,"%d  %e  \n", i, Pol[i]/*, Qol[i]*/);
			//	PCYCLE++;
			//	hout_p = 0.0;
			//}			
		}
		hout += h*180/Pi;
		hpipe += h*180/Pi;
		houtlet += h*180/Pi; 
		hout_p += h*180/Pi;
		Derivs(phi, P, dP, Pol, Qol);
		for(i=1; i<=NRK; i++) Pscal[i] = fabs(P[i]) + fabs(h*dP[i]) + TINY;
		RK4(Pol, Qol, P, dP, NRK, &phi, h, Pscal, &hdid, &hnext, Derivs);
		for(i=1; i<=N+2; i++) if(P[i] < 0) P[i] = 0;


		if(phi>2*Pi){
			Qvo = 0.0, Qvi = 0.0;		// Outlet flow
			Ao1 = 0., Ao2 = 0.;
			for(i=1; i<=N; i++){
				ph[i] = phi + 2*Pi/N*(i-1);
				while(ph[i] > 2*Pi) ph[i] = ph[i] - 2*Pi;
				if(ph[i]<=Pi){
					if(P[i]>Pol[1]) Qvo = Qvo + Cvo*Avo(ph[i])*sqrt(2/rho)*sqrt(P[i]-Pol[1]);
					else Qvo = Qvo - Cvo*Avo(ph[i])*sqrt(2/rho)*sqrt(Pol[1]-P[i]);
					Ao1 = Ao1 + Avo(ph[i]);
				}
				else{
					if(P[i]>P[N+2]) Qvi = Qvi + Cvi*Avo(ph[i])*sqrt(2/rho)*sqrt(P[i]-P[N+2]);
					else Qvi = Qvi - Cvi*Avo(ph[i])*sqrt(2/rho)*sqrt(P[N+2]-P[i]);
					Ao2 = Ao2 + Avo(ph[i]);
				}		
			}
			Qvo = Qvo - pow(hv,3)/12/eta*(Pol[1]-Ph)*(1/log(1+LB1/(R-rV1))+1/log(1+LB2/(R+rV1)))*(2*R*rV1*thetaV2-Ao1)/2/R/rV1;
			Qvi = Qvi - pow(hv,3)/12/eta*(P[N+2]-Ph)*(1/log(1+LB1/(R-rV2))+1/log(1+LB2/(R+rV2)))*(2*R*rV2*thetaV5-Ao2)/2/R/rV2;

			//if(houtlet>=olphi*180/Pi){
			//	Outlet(phi, Pol, Qol, P, Pp, Qp);
			//	houtlet = 0.0;
			//}
		}


		//if(hpipe>=pphi*180/Pi){			
		//	Pipe(phi, Pp, Qp, P, Pol, Qol);
		//	hpipe = 0.0;
		//}
		
		//STorq(P,ST,phi);
		h = hnext;
	}
		
	fclose(fp);
	fclose(gp);

	free_dvector(P,1,NRK);
	free_dvector(dP,1,NRK);
	free_dvector(Pscal,1,NRK);
	free_dvector(ph,1,N);
	free_dvector(ST,1,3);
	free_dvector(Pp,1,Np);
	free_dvector(Qp,1,Np);
	free_dvector(Pol,1,Nol);
	free_dvector(Qol,1,Nol);

	return 0;
}
/* Read parameters */
void Input()
{
	FILE *ep;

	if((ep = fopen("Input_pressure_475.dat", "r"))==NULL) {
		printf("Cannot open file.\n");
		exit(1);
	}

	fscanf(ep, "%*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*c %lf %*c %lf %*c %lf %*c %lf %*c %lf %*s %lf %*s %lf %*c %lf %*s %lf %*c %lf	%*s %lf %*c %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*c %lf %*s %lf %*s %lf %*c %lf %*s %lf %*c %lf %*s %lf %*c %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf", 
		&R, &Rp, &Rpi, &Lpi, &Lc, &Ld, &thetaK, &thetaV1, &thetaV2, &thetaV3, &thetaV4, &thetaV5, &thetaV6, &rB, &rV1, &rV2, &LB1, &LB2, 
		&Cvo, &Cvi, &Cs, &Pin, &Pout, &Ph, &rho, &eta, &K, &Rs1, &Rs2, &Rso, &beta0, &beta2, &Vvo, &Vvi, &Avout, &Avin, &Lvo, &Kvo, &Cvin, &phif, &phiout, &w);
	printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n", 
			R, Rp, Rpi, Lpi, Lc, Ld, thetaK, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6, rB, rV1, rV2, LB1, LB2, 
			Cvo, Cvi, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta0, beta2, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvin, phif, phiout, w);

	fclose(ep);
}
/*Initial condition */
void Init(double *P, double *Pp, double *Qp, double *Pol, double *Qol)
{
	int i;

	Ap = Pi*Rp*Rp;
	Vpi = Pi*Rpi*Rpi*Lpi;
	hp = 0.e-6;
	hs = 0.e-6;
	hv = 0.e-6;
	thetaB = thetaK/2*Pi/180;
	thetaV1 = thetaV1*Pi/180;
	thetaV2 = thetaV2*Pi/180;
	thetaV3 = thetaV3*Pi/180;
	thetaV4 = thetaV4*Pi/180;
	thetaV5 = thetaV5*Pi/180;
	thetaV6 = thetaV6*Pi/180;
	beta0 = beta0*Pi/180;
	beta2 = beta2*Pi/180;
	phif = phif*2*Pi;
	V0 = Vpi + (Ld+2*R*tan(beta0))*Ap;	
	//phiout = phiout*Pi/180;	
	Apipe = Pi*Rpipe*Rpipe;	

	P[1] = Pout*0.;
	P[2] = Pout*1;
	P[3] = Pout*1;
	P[4] = Pout*1;
	P[5] = Pout*1;
	P[6] = Pin*0.01;
	P[7] = Pin*1;
	P[8] = Pin*1;
	P[9] = Pin*0.01;
	P[10] = Pout*1.0;
	P[11] = Pin*1.0;
	P[12] = Ap*2*R*tan(beta0)*N*w/60;
	P[13] = Pout*0.0;
	P[14] = Ap*2*R*tan(beta0)*N*w/60/*0.0000000015*/;
	P[15] = Pout*0.98;
	
	for(i=1; i<=Np; i++){
		Pp[i] = Pout+0.01*Pout*sin((i-1)*2*Pi);
		Qp[i] = Ap*2*R*tan(beta0)*N*w/60*(1+0.01*sin((i-1)*2*Pi));
	}
	for(i=1; i<=Nol; i++){
		Pol[i] = Pout+0.001*Pout*sin((i-1)*2*Pi);
		Qol[i] = Ap*2*R*tan(beta0)*N*w/60*(1+0.001*sin((i-1)*2*Pi));
	}
	
	w = w*2*Pi/60;	
	c = sqrt(K/rho)/sqrt(1+K*2*Rol/Ep/ep);
	c2 = sqrt(K/rho)/sqrt(1+K*2*Rpipe/Ep/ep2);
	olphi = w*Lol/c/Nol;
	pphi = w*Lpipe/c2/Np;
}

/* Bulk modulus */
double Kt(double Pt)
{
	double Kt;

	if(Pt<=5e6) Kt = (14.5*5e6+(K-14.5*19.5e6))/5e6*Pt;
	else if(Pt>5e6 && Pt<=19.5e6) Kt = 14.5*Pt+(K-14.5*19.5e6);
	else if(Pt>19.5e6) Kt = K;

	return K;
}

/* R-K 4th */
void RK4(double *Pol, double *Qol, double *P, double *dP, int n, double *phiR, double htry, double *Pscal, double *hdid, double *hnext, void (*Derivs)(double, double *, double *, double *, double *))
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
	Pt =  dvector(1,n);
	Psav =  dvector(1,n);
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
		(*Derivs)(phisav+A2X*h, Pt, K2, Pol, Qol);
		for(i=1; i<=n; i++)	Pt[i] = Psav[i] + (A31*K1[i]+A32*K2[i])*h;
		(*Derivs)(phisav+A3X*h, Pt, K3, Pol, Qol);
		for(i=1; i<=n; i++)	Pt[i] = Psav[i] + (A41*K1[i]+A42*K2[i]+A43*K3[i])*h;	
		(*Derivs)(phisav+A4X*h, Pt, K4, Pol, Qol);
		for(i=1; i<=n; i++)	Pt[i] = Psav[i] + (A51*K1[i]+A52*K2[i]+A53*K3[i]+A54*K4[i])*h;
		(*Derivs)(phisav+A5X*h, Pt, K5, Pol, Qol);
		for(i=1; i<=n; i++)	Pt[i] = Psav[i] + (A61*K1[i]+A62*K2[i]+A63*K3[i]+A64*K4[i]+A65*K5[i])*h;
		(*Derivs)(phisav+A6X*h, Pt, K6, Pol, Qol);
		for(i=1; i<=n; i++){
			P[i] = Psav[i] + (B1*K1[i]+B3*K3[i]+B4*K4[i]+B6*K6[i])*h;
			err[i] = (E1*K1[i] + E3*K3[i] + E4*K4[i] + E5*K5[i] + E6*K6[i])*h;			
		}
		//printf("%e \n",h);
		*phiR = phisav + h;
		if(*phiR == phisav/2) nrerror("stepsize not significant in stiff");
		errmax = 0.0;			//Evaluate accuracy
		for(i=1; i<=n; i++) errmax = FMAX(errmax, fabs(err[i]/Pscal[i]));
		errmax /= eps;			//scale relative to the required tolerance 
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
		else{				//Truncation error too large, reduce stepsize
			*hnext = SAFETY*h*pow(errmax,PSHRNK);
			h = (h >= 0.0 ? FMAX(*hnext,SHRNK*h) : FMIN(*hnext,SHRNK*h));
		}
	}
	nrerror("exceeded MAXTRY in RKstep");
}
/* dP/dphi */
void Derivs(double phid, double *Pt, double *dP, double *Pol, double *Qol)
{
	double Kt(double );

	int i, sgn1, sgn2, sgn3;
	double Cv, *Pr, Pv, *phiD, *Lpc, Ao1, Ao2, Qvo, Qvi,Apo, Lvp=50e-3, Ain, Aout, bv1, bv2, bvo=100e-3, Lv1=300e-3, Rpipei=10.0e-3, Lv2=200e-3, Lv3=1000e-3, Lv4=1000e-3, tanv, Bv, Cvpo;
	double beta, dbeta, wd, wn=60*2*Pi, zeta=1/1.414/2;
	double Fin, Fp;

	phiD = dvector(1,N);
	Pr = dvector(1,N);
	Lpc = dvector(1,N);	

	//if(phid>2*Pi){
	//	wd = wn*sqrt(1-zeta*zeta);
	//	beta = beta0*(1+0.2*exp(-zeta*wn*phid/w)*(cos(wd*phid/w)-zeta*wn/wd*sin(wd*phid/w)));
	//	dbeta = -beta0*0.2*(zeta*wn*exp(-zeta*wn*phid/w)*(cos(wd*phid/w)-zeta*wn/wd*sin(wd*phid/w))+exp(-zeta*wn*phid/w)*(wd*sin(wd*phid/w)+zeta*wn*cos(wd*phid/w)));
	//}
	//else{
	//	beta = beta0*1.2;
	//	dbeta = 0.0;
	//}
	beta = beta0*1.0;
	dbeta = 0.0;


	while(phid > 2*Pi) phid = phid - 2*Pi;	

	for(i=1; i<=N; i++){
		phiD[i] = phid + 2*Pi/N*(i-1);
		while(phiD[i] > 2*Pi) phiD[i] = phiD[i] - 2*Pi;	
		//hv = hv - 40e-6*sin(phiD[i]);
		Lpc[i] = Lc - Ld - R*(tan(beta)*(1+cos(phiD[i]))-tan(beta2)/cos(beta)*sin(phiD[i]));
		if(Pt[i]>Ph) Pr[i] = Ph + 2*(Pt[i]-Ph)/(1+sqrt(1+pow(Pi*pow(hs,3)/6/eta/log(Rs2/Rs1),2)*2*rho*(Pt[i]-Ph)/pow(Cs,2)/pow(Pi*Rso*Rso,2)));
		else Pr[i] = Ph + 2*(Pt[i]-Ph)/(1+sqrt(1+pow(Pi*pow(hs,3)/6/eta/log(Rs2/Rs1),2)*2*rho*(Ph-Pt[i])/pow(Cs,2)/pow(Pi*Rso*Rso,2)));
		//theta12(phiD[i]);
		Cv = (phiD[i]<=Pi?Cvo:Cvi);
		sgn1 = (Pt[i]>=Pt[N+1]?1:-1);
		sgn2 = (Pt[i]>=Pr[i]?1:-1);
		sgn3 = (Pt[i]>=Pt[N+2]?1:-1);
		//dP[i] = -K/w/(V0-R*(tan(beta)*(1-cos(phiD[i]))+tan(beta2)/cos(beta)*sin(phiD[i]))*Ap)
		//		*( -w*R*(tan(beta)*sin(phiD[i])+tan(beta2)/cos(beta)*cos(phiD[i]))*Ap + Cv*Avo(phiD[i])*sqrt(2/rho)*sqrt(fabs(Pt[i]-Pv))*sgn1
		//		+ 2*Pi*Rp*pow(hp,3)/12/eta/Lpc[i]*(Pt[i]-Ph)-2*Pi*Rp/2*hp*w*R*(tan(beta)*sin(phiD[i])+tan(beta2)/cos(beta)*cos(phiD[i]))
		//		+ pow(hv,3)/12/eta*(Pt[i]-Ph)*(1/log(1+LB1/(R-rB))+1/log(1+LB2/(R+rB)))*2*thetaB  
		//		+ Cs*Pi*Rso*Rso*sqrt(2/rho)*sqrt(fabs(Pt[i]-Pr[i]))*sgn2 );
		dP[i] = -Kt(Pt[i])/w/(V0-R*(tan(beta)*(1-cos(phiD[i]))+tan(beta2)/cos(beta)*sin(phiD[i]))*Ap)
				*( -R*(dbeta/pow(cos(beta),2)*(1-cos(phiD[i]))+dbeta*tan(beta2)*sin(beta)/cos(beta)/cos(beta)*sin(phiD[i])+w*(tan(beta)*sin(phiD[i])+tan(beta2)/cos(beta)*cos(phiD[i])))*Ap 
				+ Cv*Avo(phiD[i])*sqrt(2/rho)*sqrt(fabs(Pt[i]-Pt[N+1]))*sgn1 + Cv*Avi(phiD[i])*sqrt(2/rho)*sqrt(fabs(Pt[i]-Pt[N+2]))*sgn3
				/*+ 2*Pi*Rp*pow(hp,3)/12/eta/Lpc[i]*(Pt[i]-Ph)-2*Pi*Rp/2*hp*R*(dbeta/pow(cos(beta),2)*(1-cos(phiD[i]))+w*(tan(beta)*sin(phiD[i])+tan(beta2)/cos(beta)*cos(phiD[i])))
				+ pow(hv,3)/12/eta*(Pt[i]-Ph)*(1/log(1+LB1/(R-rB))+1/log(1+LB2/(R+rB)))*2*thetaB 
				+ Cs*Pi*Rso*Rso*sqrt(2/rho)*sqrt(fabs(Pt[i]-Pr[i]))*sgn2 */);
		//printf("%f  ",dP[i]);
	}

	Qvo = 0., Qvi = 0.;
	Ao1 = 0., Ao2 = 0.;
	Fin = 0., Fp = 0.;	

	for(i=1; i<=N; i++){
		if(phiD[i]<=Pi){
			if(Pt[i]>Pt[N+1]) Qvo = Qvo + Cvo*Avo(phiD[i])*sqrt(2/rho)*sqrt(Pt[i]-Pt[N+1]);
			else Qvo = Qvo - Cvo*Avo(phiD[i])*sqrt(2/rho)*sqrt(Pt[N+1]-Pt[i]);
			Ao1 = Ao1 + Avo(phiD[i]);			 
			Fp = Fp + Pt[i]*Avo(phiD[i]);
			Fin = Fin + Cvo*Cvo*Avo(phiD[i])*2*(Pt[i]-Pt[N+1]);
		}
		else{
			if(Pt[i]>Pt[N+2]) Qvi = Qvi + Cvi*Avi(phiD[i])*sqrt(2/rho)*sqrt(Pt[i]-Pt[N+2]);
			else Qvi = Qvi - Cvi*Avi(phiD[i])*sqrt(2/rho)*sqrt(Pt[N+2]-Pt[i]);
			Ao2 = Ao2 + Avi(phiD[i]);
		}
	}
	Qvo = Qvo - pow(hv,3)/12/eta*(Pt[N+1]-Ph)*(1/log(1+LB1/(R-rV1))+1/log(1+LB2/(R+rV1)))*(2*R*rV1*thetaV2-Ao1)/2/R/rV1;
	Qvi = Qvi - pow(hv,3)/12/eta*(Pt[N+2]-Ph)*(1/log(1+LB1/(R-rV2))+1/log(1+LB2/(R+rV2)))*(2*R*rV2*thetaV5-Ao2)/2/R/rV2;

	
	bv1 = R*thetaV2;
	tanv = (bv1-bvo)/2/Lvp;	
	Apo = bvo*2*rV1;
	Aout = Pi*Rol*Rol;
	Ain = Pi*Rpipei*Rpipei;
	Cvpo = 0.7;

	/* Single Piston Model */
	//dP[N+1] = 0.0;
	//dP[N+2] = 0.0;
	//printf("%f %f \n", 1/rho/w/((bv1-Lvp*tanv)*Lvp*2*rV1)*(Fin-rho/Aout*Pt[N+3]*fabs(Pt[N+3])+Fp-Pp[1]*Aout)*Aout, (Qvo-Pt[N+3])/w/((bv1-Lvp*tanv)*Lvp*2*rV1)*Pt[N+3]);

	///* Simplest Shape with Piping */
	//sgn1 = (Pt[N+3]>=0?0.1:2);
	dP[N+1] = Kt(Pt[N+1])/w/((bv1-Lvp*tanv)*Lvp*2*rV1)*(Qvo - /*Qol[1]*//*Pt[N+3]*/Cvo*1.75e-6*sqrt(2/rho)*sqrt(fabs(Pt[N+1]-Ph)));
	//dP[N+3] = Ain/Lv3/w/rho*((Pt[N+1]-0.965*Pout/*Pol[1]*/) /*- 5/pow(Ain,2)*rho/2*Pt[N+3]*fabs(Pt[N+3])*/ - 3/pow(Aout*0.02,2)*rho/2*Pt[N+3]*fabs(Pt[N+3]) );
	//dP[N+3] = 1/rho/w/((bv1-Lvp*tanv)*Lvp*2*rV1)*(Fin-rho/Aout*Pt[N+3]*fabs(Pt[N+3])+Fp-Pol[1]*Aout+Pt[N+1]*(Aout-Ao1))*Aout - (Qvo-Pt[N+3])/w/((bv1-Lvp*tanv)*Lvp*2*rV1)*Pt[N+3];
	if(Pt[N+2]>Pin) dP[N+2] = Kt(Pt[N+2])/w/Vvi*(Qvi-Cvin*10.5e-6*sqrt(2/rho)*sqrt(Pt[N+2]-Pin));  
	else dP[N+2] = Kt(Pt[N+2])/w/Vvi*(Qvi+Cvin*10.5e-6*sqrt(2/rho)*sqrt(Pin-Pt[N+2]));

	/* Internal impedance modeling with piping */
	

	/* Compressibility in front of orifice */
	//float k = 1e-12;
	//dP[N+1] = Kt(Pt[N+1])/w/((bv1-Lvp*tanv)*Lvp*2*rV1)*(Qvo-Pt[N+3]);
	//dP[N+3] = Aout/Lv3/w/rho*((Pt[N+1]-Pt[N+4]) /*- 1/pow(Aout,2)*rho/2*Pt[N+3]*fabs(Pt[N+3])*/);
	//if(Pt[N+4]>Ph) dP[N+4] = Kt(Pt[N+4])/w/Aout/Lv4*(Pt[N+3]-Cvo*10e-6*sqrt(2/rho)*sqrt(Pt[N+4]-Ph)); 
	//else dP[N+4] = Kt(Pt[N+4])/w/Aout/Lv4*(Pt[N+3]+Cvo*10e-6*sqrt(2/rho)*sqrt(Ph-Pt[N+4]));
	////dP[N+2] = 0.0;
	//if(Pt[N+2]>Pin) dP[N+2] = Kt(Pt[N+2])/w/Vvi*(Qvi-Cvin*Avin*sqrt(2/rho)*sqrt(Pt[N+2]-Pin));  
	//else dP[N+2] = Kt(Pt[N+2])/w/Vvi*(Qvi+Cvin*Avin*sqrt(2/rho)*sqrt(Pin-Pt[N+2]));

	/* Relief Valve Modeling */
	//float k = 1e-12;
	//dP[N+1] = K/w/((bv1-Lvp*tanv)*Lvp*2*rV1)*(Qvo-Pt[N+3]);
	//dP[N+3] = Aout/Lv3/w/rho*((Pt[N+1]-Pt[N+4]) /*- 80/pow(Aout,2)*rho/2*Pt[N+3]*fabs(Pt[N+3])*/);
	//if(Pt[N+4]>Pout*0.6) dP[N+4] = K/w/Aout/Lv4*(Pt[N+3]-Cvo*k*(Pt[N+4]-Pout)*sqrt(2/rho)*sqrt(Pt[N+4]-Ph)-0.9*Aout*sqrt(2/rho)*sqrt(Pt[N+4]-Pout*0.6)); 
	//else dP[N+4] = K/w/Aout/Lv4*(Pt[N+3]/*+0.9*Aout*sqrt(2/rho)*sqrt(Pt[N+4])*/);	
	//if(Pt[N+2]>Pin) dP[N+2] = K/w/Vvi*(Qvi-Cvin*Avin*sqrt(2/rho)*sqrt(Pt[N+2]-Pin));  
	//else dP[N+2] = K/w/Vvi*(Qvi+Cvin*Avin*sqrt(2/rho)*sqrt(Pin-Pt[N+2]));


	/* Port Momentum Modeling */
	//dP[N+1] = K/w/((bv1-Lv1*tanv)*Lv1*2*rV1)*(Qvo-Pt[N+6]);  //P1
	//if(Pt[N+5]>Pt[N+4]) dP[N+5] = K/w/((bv1-(2*Lv1+2*Lv2+Lv4)*tanv)*Lv4*2*rV1)*(-Cvpo*Apo*sqrt(2/rho)*sqrt(Pt[N+5]-Pt[N+4])+Pt[N+6]);
	//else dP[N+5] = K/w/((bv1-(2*Lv1+2*Lv2+Lv4)*tanv)*Lv4*2*rV1)*(Cvpo*Apo*sqrt(2/rho)*sqrt(Pt[N+4]-Pt[N+5])+Pt[N+6]);
	//if(Pt[N+5]>Pt[N+4]) dP[N+4] = K/w/Aout/Lv5*(Cvpo*Apo*sqrt(2/rho)*sqrt(Pt[N+5]-Pt[N+4])-Pt[N+3]);
	//else dP[N+4] = K/w/Aout/Lv5*(-Cvpo*Apo*sqrt(2/rho)*sqrt(Pt[N+4]-Pt[N+5])-Pt[N+3]);
	//dP[N+6] = ((bv1-(2*Lv1+Lv2)*tanv)*2*rV1)/Lv2/w/rho*((Pt[N+1]-Pt[N+5]) -100*rho/2*Pt[N+6]*fabs(Pt[N+6])/pow(((bv1-(2*Lv1+Lv2)*tanv)*2*rV1),2) ) ;
	//dP[N+3] = Aout/Lv3/w/rho*((Pt[N+4]-Pout) - 150*rho/2*Pt[N+3]*fabs(Pt[N+3])/pow(Aout,2) );
	//dP[N+4] = Ain/Lv2/w/rho*((Pin-Pt[N+2]) -64/(rho*2*Rpipei/eta)/2*rho*Pt[N+4]/Ain*Lv2/2/Rpipei -3*rho/2*Pt[N+4]*fabs(Pt[N+4])/pow(Ain,2) );
	//if(Pt[N+5]>Ph) dP[N+5] = K/w/Aout/Lv1*(Pt[N+3]-Cvo*Aout*0.0365*sqrt(2/rho)*sqrt(Pt[N+5]-Ph));
	//else dP[N+5] = K/w/Aout/Lv1*(Pt[N+3]+Cvo*Aout*0.0365*sqrt(2/rho)*sqrt(Ph-Pt[N+5]));
	//if(Pt[N+4]>Pout) dP[N+4] = K/w/Aout/Lv4*(Pt[N+3]-Cvo*Avout*sqrt(2/rho)*sqrt(Pt[N+4]-Pout));  //Pt2 orifce model
	//else dP[N+4] = K/w/Aout/Lv4*(Pt[N+3]+Cvo*Avout*sqrt(2/rho)*sqrt(Pout-Pt[N+4]));	
	//dP[N+2] = K/w/Vvi*(Pt[N+4]-Qvi);  //P2

	//dP[N+3] = -2*N*(0.1*Pout)*sin(N*phid);

	//printf("%lf %.6e\n", phiD[1]*180/Pi, Qvo);

	free_dvector(phiD,1,N);
	free_dvector(Pr,1,N);
	free_dvector(Lpc,1,N);
}
/* Outlet Pathway Pulse */
void Outlet(double phi, double *Pol, double *Qol, double *P, double *Pp, double *Qp)
{
	int i;
	double alpha, Cp, Cm, BL, zeta=0.8, Res, BL2, zeta2=0.8, Res2;
	//Rol = Rpipe;
	Aol = Pi*Rol*Rol;

	printf("**** \n");	
	alpha = 0.6*11.8e-6*sqrt(2/rho);
	BL2 = c/g/Apipe;
	Res2 = zeta2*Lpipe/Np/2/g/2/Rpipe/Apipe;
	BL = c/g/Aol;
	Res = zeta*Lol/Nol/2/g/2/Rol/Aol;
	for(i=2; i<Nol; i++){
		Cp = Pol[i-1]/rho/g+BL*Qol[i-1]-Res*Qol[i-1]*fabs(Qol[i-1]);
		Cm = Pol[i+1]/rho/g-BL*Qol[i+1]+Res*Qol[i+1]*fabs(Qol[i+1]);
		Qol[i] = (Cp-Cm)/2/BL;
		Pol[i] = rho*g*(Cp-BL*Qol[i]);
	}
	Qol[1] = P[N+3];
	Pol[1] = rho*g*(Pol[2]/rho/g-BL*Qol[2]+Res*Qol[2]*fabs(Qol[2])+BL*Qol[1]);	
	//Pol[1] = P[N+1];
	//Qol[1] = Pol[1]/rho/g/BL-(Pol[2]/rho/g-BL*Qol[2]+Res*Qol[2]*fabs(Qol[2]))/BL;
	//Qol[Nol] = Qp[1];
	//Pol[Nol] = rho*g*(Pol[Nol-1]/rho/g+BL*Qol[Nol-1]-Res*Qol[Nol-1]*fabs(Qol[Nol-1])-BL*Qol[Nol]);
	//Pol[Nol] = Pout;
	//Qol[Nol] = -Pol[Nol]/rho/g/BL+(Pol[Nol-1]/rho/g+BL*Qol[Nol-1]-Res*Qol[Nol-1]*fabs(Qol[Nol-1]))/BL;
	Qol[Nol] = (-alpha*alpha*BL*rho*g+sqrt(pow(alpha*alpha*BL*rho*g,2)-4*alpha*alpha*(Pin-(Pol[Nol-1]/rho/g+BL*Qol[Nol-1]-Res*Qol[Nol-1]*fabs(Qol[Nol-1]))*rho*g)))/2;
	Pol[Nol] = rho*g*(Pol[Nol-1]/rho/g+BL*Qol[Nol-1]-Res*Qol[Nol-1]*fabs(Qol[Nol-1])-BL*Qol[Nol]);

	//Cp = Pol[Nol-1]/rho/g+BL*Qol[Nol-1]-Res*Qol[Nol-1]*fabs(Qol[Nol-1]);
	//Cm = Pp[2]/rho/g-BL2*Qp[2]+Res2*Qp[2]*fabs(Qp[2]);
	//Qol[Nol] = (Cp-Cm)/(BL+BL2);
	//Pol[Nol] = rho*g*(Cp-BL*Qol[Nol]);
	
}
/* Pipe Pulse */
void Pipe(double phi, double *Pp, double *Qp, double *P, double *Pol, double *Qol)
{
	int i;	
	double alpha, Cp, Cm, BL, zeta=0.8, Res;

	printf("---- \n");	
	alpha = 0.6*11.8e-6*sqrt(2/rho); /*15.8, 11.4,11.8, 8.8, 1.65*/
	BL = c2/g/Apipe;
	Res = zeta*Lpipe/Np/2/g/2/Rpipe/Apipe;
	for(i=2; i<Np; i++){
		Cp = Pp[i-1]/rho/g+BL*Qp[i-1]-Res*Qp[i-1]*fabs(Qp[i-1]);
		Cm = Pp[i+1]/rho/g-BL*Qp[i+1]+Res*Qp[i+1]*fabs(Qp[i+1]);
		Qp[i] = (Cp-Cm)/2/BL;
		Pp[i] = rho*g*(Cp-BL*Qp[i]);
	}
	//Pp[1] = Pol[Nol];
	//Qp[1] = Qol[Nol];
	Pp[1] = P[N+1];
	Qp[1] = Pp[1]/rho/g/BL-(Pp[2]/rho/g-BL*Qp[2]+Res*Qp[2]*fabs(Qp[2]))/BL;
	//Qp[1] = Qol[Nol];
	//Pp[1] = rho*g*(Pp[2]/rho/g-BL*Qp[2]+Res*Qp[2]*fabs(Qp[2])+BL*Qp[1]);
	//Pp[1] = Pol[Nol];
	//Qp[1] = Pp[1]/rho/g/BL-(Pp[2]/rho/g-BL*Qp[2]+Res*Qp[2]*fabs(Qp[2]))/BL;
	Qp[Np] = (-alpha*alpha*BL*rho*g+sqrt(pow(alpha*alpha*BL*rho*g,2)-4*alpha*alpha*(Pin-(Pp[Np-1]/rho/g+BL*Qp[Np-1]-Res*Qp[Np-1]*fabs(Qp[Np-1]))*rho*g)))/2;
	Pp[Np] = rho*g*(Pp[Np-1]/rho/g+BL*Qp[Np-1]-Res*Qp[Np-1]*fabs(Qp[Np-1])-BL*Qp[Np]);
	//Pp[Np] = Pout;
	//Qp[Np] = -Pout/rho/g/BL+(Pp[Np-1]/rho/g+BL*Qp[Np-1]-Res*Qp[Np-1]*fabs(Qp[Np-1]))/BL;

}

/* Leakage */
void Leakage(double phi, double *Pt, double *Qlp, double *Qv, double *Qs, double *Qpt, double *Qvt, double *Qst)
{
	double *Lpc, Pv, *Pr;
	double Qvo, Qvi, Ao1, Ao2;
	double *phiL, beta;

	beta = beta0;
	Pr = dvector(1,N);
	Lpc = dvector(1,N);
	phiL = dvector(1,N);
	//leakage in a cylinder
	Lpc[1] = Lc - Ld - R*(tan(beta)*(1+cos(phi))-tan(beta2)/cos(beta)*sin(phi));
	if(Pt[1]>Ph) Pr[1] = Ph + 2*(Pt[1]-Ph)/(1+sqrt(1+pow(Pi*pow(hs,3)/6/eta/log(Rs2/Rs1),2)*2*rho*(Pt[1]-Ph)/pow(Cs,2)/pow(Pi*Rso*Rso,2)));
	else Pr[1] = Ph + 2*(Pt[1]-Ph)/(1+sqrt(1+pow(Pi*pow(hs,3)/6/eta/log(Rs2/Rs1),2)*2*rho*(Ph-Pt[1])/pow(Cs,2)/pow(Pi*Rso*Rso,2)));
	if(phi<=Pi) Pv = Pt[N+1];
	else Pv = Pt[N+2];
	*Qlp = 2*Pi*Rp*pow(hp,3)/12/eta/Lpc[1]*(Pt[1]-Ph)-2*Pi*Rp/2*hp*w*R*(tan(beta)*sin(phi)+tan(beta2)/cos(beta)*cos(phi));
	*Qv = pow(hv,3)/12/eta*(Pt[1]-Ph)*(1/log(1+LB1/(R-rB))+1/log(1+LB2/(R+rB)))*2*thetaB;
	*Qs = Cs*Pi*Rso*Rso*sqrt(2/rho)*sqrt(fabs(Pt[1]-Pr[1]));	
	//total leakage 
	*Qpt = 0., *Qst = 0.;
	Qvo = 0., Qvi = 0.;
	Ao1 = 0., Ao2 = 0.;
	for(int i=1; i<=N; i++){
		phiL[i] = phi + 2*Pi/N*(i-1);
		while(phiL[i] > 2*Pi) phiL[i] = phiL[i] - 2*Pi;
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
	Qvo = Qvo - pow(hv,3)/12/eta*(Pt[N+1]-Ph)*(1/log(1+LB1/(R-rV1))+1/log(1+LB2/(R+rV1)))*(2*R*rV1*thetaV2)/2/R/rV1;
	Qvi = Qvi - pow(hv,3)/12/eta*(Pt[N+2]-Ph)*(1/log(1+LB1/(R-rV2))+1/log(1+LB2/(R+rV2)))*(2*R*rV2*thetaV5)/2/R/rV2;
	*Qvt = Qvo + Qvi;

	free_dvector(Pr,1,N);
	free_dvector(Lpc,1,N);
	free_dvector(phiL,1,N);
}
/* Point of action of the resultant force & Moment */
void Point_resul(double phi, double *Pt, double *eBVx, double *eBVy)
{
	double Pe=0.0;	
	double *phiP;
	
	phiP = dvector(1,N);
	*eBVx = 0.0, *eBVy = 0.0;
	for(int i=1; i<=N; i++){
		phiP[i] = phi + 2*Pi/N*(i-1);
		while(phiP[i] > 2*Pi) phiP[i] = phiP[i] - 2*Pi;
		*eBVx = *eBVx + Pt[i]*R*sin(phiP[i]);
		*eBVy = *eBVy + Pt[i]*R*cos(phiP[i]);
		Pe = Pe + Pt[i];
	}
	*eBVx = *eBVx/Pe;
	*eBVy = *eBVy/Pe;

	free_dvector(phiP,1,N);
}
//Torque of the Swash-plate
//void STorq(double *Pt, double *STt, double phis)
//{
//	int i;
//	double *phiS;
//
//	phiS = dvector(1,N);
//	
//	for(i=1; i<=3; i++) STt[i] = 0.0;
//	//while(phis >= 2*Pi) phis = phis - 2*Pi;
//
//	for(i=1; i<=N; i++){
//		phiS[i] = phis + 2*Pi/N*(i-1);
//		if(phiS[i] > 2*Pi) phiS[i] = phiS[i] - 2*Pi;
//		STt[1] = STt[1]  + Pt[i]*Ap/cos(beta0)*/*(cos(beta)+tan(beta)*sin(beta))**/R*cos(phiS[i])/cos(beta0);
//		STt[2] = STt[2]  + Pt[i]*Ap/cos(beta0)*R*sin(phiS[i]);
//		STt[3] = 0.0;
//		//printf("%e \n", STt[1]);
//	}
//
//	free_dvector(phiS,1,N);
//}