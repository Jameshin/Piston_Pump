#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "nrutil.h"
#include "var.h"

extern void theta12(double );
extern double Avo(double );
extern double Avi(double );
extern double integB(double );

extern double R, Rp, Lc, Ld, rB, rV1, rV2, LB1, LB2, Cv, Cs, Pin, Pout, Ph, rho, eta, K, Rs1, Rs2, Rso, beta, beta2, phif, phiout, w;
extern double V0, Ap, Ao, Vvo, Vvi, Avout, Avin, Lvo, Kvo, Cvi, thetaK, thetaB, thetaV1, thetaV2, thetaV3, thetaV4, thetaV5, thetaV6;
extern double integB1, integB2, integB3, hp, hs, hv, theta1, theta2, theta_1;

/* Valve opening area */
//double Avo(double phiA)
//{
//	double Ao = 0., Ai = 0., rV;
//
//	if(phiA < 0) phiA = phiA + 2*Pi;
//	while(phiA > 2*Pi) phiA = phiA - 2*Pi;
//	if(phiA <= Pi) rV = rV1;
//	else rV = rV2;
//
//	theta_1 = asin(rV/rB);
//
//	//printf("%lf %.10e \n", phiA*180/Pi, Ao);
//
//	if(phiA >= 0 && phiA <= thetaV1-thetaB) Ao = 0.;
//	//else if(phiA >= thetaV1-thetaB-10*Pi/180 && phiA <= thetaV1-thetaB+1*Pi/180) Ao = 2e-6;
//	else if(phiA >= thetaV1-thetaB && phiA <= thetaV1+thetaB+(rV1-rB*(1-cos(theta_1)))/R)	
//		Ao = (4*R*rV1*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2))/(2*thetaB+(rV1-rB*(1-cos(theta_1)))/R)*(phiA-thetaV1+thetaB);
//	else if(phiA >= thetaV1+thetaB+(rV1-rB*(1-cos(theta_1)))/R && phiA <= thetaV1+thetaV2-thetaB-(rV1-rB*(1-cos(theta_1)))/R)
//		Ao = 4*R*rV1*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2);
//	else if(phiA >= thetaV1+thetaV2-thetaB-(rV1-rB*(1-cos(theta_1)))/R && phiA <= Pi-(thetaV3-thetaB))
//		Ao = (4*R*rV1*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2))/(thetaV1+thetaV2-thetaB-(rV1-rB*(1-cos(theta_1)))/R-Pi+(thetaV3-thetaB))*(phiA-Pi+(thetaV3-thetaB));
//	else if(phiA >= Pi-(thetaV3-thetaB) && phiA <= Pi+(thetaV4-thetaB)) Ao = 0.;
//
//	else if(phiA >= Pi+(thetaV4-thetaB) && phiA <= Pi+thetaV4+thetaB+(rV2-rB*(1-cos(theta_1)))/R)	
//		Ao = (4*R*rV2*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2))/(2*thetaB+(rV2-rB*(1-cos(theta_1)))/R)*(phiA-Pi-(thetaV4-thetaB));
//	else if(phiA >= Pi+thetaV4+thetaB+(rV2-rB*(1-cos(theta_1)))/R && phiA <= Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R)
//		Ao = 4*R*rV2*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2);
//	else if(phiA >= Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R && phiA <= 2*Pi-(thetaV6-thetaB))
//		Ao = (4*R*rV2*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2))/(Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R-2*Pi+(thetaV6-thetaB))*(phiA-2*Pi+(thetaV6-thetaB));
//	else if(phiA >= 2*Pi-(thetaV6-thetaB) && phiA <= 2*Pi) Ao = 0.;
//
//	return Ao;
//}
double Avo(double phiA)
{
	double Ao = 0., Ai = 0., rV;

	if(phiA < 0) phiA = phiA + 2*Pi;
	while(phiA > 2*Pi) phiA = phiA - 2*Pi;
	rV = rV1;

	theta_1 = asin(rV/rB);
		
	// Ao: outlet area
	if(phiA >= 0 && phiA <= thetaV1-thetaB/*-5*Pi/180*/) Ao = 0.;
	//else if(phiA >= thetaV1-thetaB-5*Pi/180 && phiA <= thetaV1-thetaB) Ao = 1e-6;
	else if(phiA >= thetaV1-thetaB && phiA <= thetaV1+thetaB+(rV1-rB*(1-cos(theta_1)))/R)	
		Ao = (4*R*rV1*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2))/(2*thetaB+(rV1-rB*(1-cos(theta_1)))/R)*(phiA-thetaV1+thetaB);
	else if(phiA >= thetaV1+thetaB+(rV1-rB*(1-cos(theta_1)))/R && phiA <= thetaV1+thetaV2-thetaB-(rV1-rB*(1-cos(theta_1)))/R)
		Ao = 4*R*rV1*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2);
	else if(phiA >= thetaV1+thetaV2-thetaB-(rV1-rB*(1-cos(theta_1)))/R && phiA <= Pi-(thetaV3-thetaB))
		Ao = (4*R*rV1*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2))/(thetaV1+thetaV2-thetaB-(rV1-rB*(1-cos(theta_1)))/R-Pi+(thetaV3-thetaB))*(phiA-Pi+(thetaV3-thetaB));
	else if(phiA >= Pi-(thetaV3-thetaB) && phiA <= 2*Pi+(thetaV1-thetaB)) Ao = 0.;
	else if(phiA >= 2*Pi+(thetaV1-thetaB) && phiA <= 2*Pi) 
		Ao = (4*R*rV1*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2))/(2*thetaB+(rV1-rB*(1-cos(theta_1)))/R)*(phiA-2*Pi-thetaV1+thetaB);

	//printf("%lf %.10e \n", phiA*180/Pi, Ao);

	return Ao;
}

double Avi(double phiA)
{
	double Ao = 0., Ai = 0., rV;

	if(phiA < 0) phiA = phiA + 2*Pi;
	while(phiA > 2*Pi) phiA = phiA - 2*Pi;
	rV = rV2;

	theta_1 = asin(rV/rB);

	// Ai: inlet area
	if(phiA >= 0 && phiA <= thetaB-thetaV6) Ai = (4*R*rV2*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2))/(Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R-2*Pi+(thetaV6-thetaB))*(phiA-(thetaB-thetaV6));
	else if(phiA >= thetaB-thetaV6 && phiA <= Pi+(thetaV4-thetaB)/*-5*Pi/180*/) Ai = 0.;
	//else if(phiA >= Pi+(thetaV4-thetaB)-5*Pi/180 && phiA <= Pi+(thetaV4-thetaB)) Ao = 1e-6;
	else if(phiA >= Pi+(thetaV4-thetaB) && phiA <= Pi+thetaV4+thetaB+(rV2-rB*(1-cos(theta_1)))/R)	
		Ai = (4*R*rV2*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2))/(2*thetaB+(rV2-rB*(1-cos(theta_1)))/R)*(phiA-Pi-(thetaV4-thetaB));
	else if(phiA >= Pi+thetaV4+thetaB+(rV2-rB*(1-cos(theta_1)))/R && phiA <= Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R)
		Ai = 4*R*rV2*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2);
	else if(phiA >= Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R && phiA <= 2*Pi-(thetaV6-thetaB))
		Ai = (4*R*rV2*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2))/(Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R-2*Pi+(thetaV6-thetaB))*(phiA-2*Pi+(thetaV6-thetaB));
	else if(phiA >= 2*Pi-(thetaV6-thetaB) && phiA <= 2*Pi) Ai = 0.;

	//printf("%lf %.10e \n", phiA*180/Pi, Ai);

	return Ai;
}


//double Avo(double phiA)
//{
//	double Ao=0.;
//	
//	if(phiA < 0) phiA = phiA + 2*Pi;
//	else if(phiA > 2*Pi) phiA = phiA - 2*Pi;
//
//	theta12(phiA);
//
//	if(phiA >= 0 && phiA <= thetaV1-thetaB) Ao = 0;
//	else if(phiA >= thetaV1-thetaB && phiA <= thetaV1-thetaB+(rV1+rB*(1-cos(theta_1)))/R)	
//		Ao = rB*rB*(theta1-sin(2*theta1)/2)+rV1*rV1*(theta2-sin(2*theta2)/2);
//	else if(phiA >= thetaV1-thetaB+(rV1+rB*(1-cos(theta_1)))/R && phiA <= thetaV1+thetaB) 
//		Ao = rB*rB*(theta_1-sin(2*theta_1)/2)+Pi/2*rV1*rV1+2*R*rV1*(phiA-thetaV1+thetaB-(rV1+rB*(1-cos(theta_1)))/R);
//	else if(phiA >= thetaV1+thetaB && phiA <= thetaV1+thetaB+(rV1-rB*(1-cos(theta_1)))/R) 
//		Ao = rB*rB*(theta_1-sin(2*theta_1)/2)+Pi/2*rV1*rV1+2*R*rV1*(phiA-thetaV1+thetaB-(rV1+rB*(1-cos(theta_1)))/R)
//			-(rV1*rV1*(theta2-sin(2*theta2)/2)-rB*rB*(theta1-sin(2*theta1)/2));
//	else if(phiA >= thetaV1+thetaB+(rV1-rB*(1-cos(theta_1)))/R && phiA <= thetaV1+thetaV2-thetaB-(rV1-rB*(1-cos(theta_1)))/R)
//		Ao = 4*R*rV1*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2);
//	else if(phiA >= thetaV1+thetaV2-thetaB-(rV1-rB*(1-cos(theta_1)))/R && phiA <= thetaV1+thetaV2-thetaB)
//		Ao = 4*R*rV1*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2)
//		    -2*R*rV1*(phiA-(thetaV1+thetaV2-thetaB-(rV1-rB*(1-cos(theta_1)))/R))+Pi/2*rV1*rV1-rB*rB*(theta_1-sin(2*theta_1)/2)
//			-(rV1*rV1*(theta2-sin(2*theta2)/2)-rB*rB*(theta1-sin(2*theta1)/2));
//	else if(phiA >= thetaV1+thetaV2-thetaB && phiA <= Pi-(thetaV3-thetaB)-(rV1+rB*(1-cos(theta_1)))/R)
//		Ao = Pi/2*rV1*rV1+4*R*rV1*(thetaB-rB*(1-cos(theta_1))/R)+rB*rB*(theta_1-sin(2*theta_1)/2)-2*R*rV1*((rV1-rB*(1-cos(theta_1)))/R+phiA-(thetaV1+thetaV2-thetaB));
//	else if(phiA >= Pi-(thetaV3-thetaB)-(rV1+rB*(1-cos(theta_1)))/R && phiA <= Pi-(thetaV3-thetaB))	
//		Ao = rB*rB*(theta1-sin(2*theta1)/2)+rV1*rV1*(theta2-sin(2*theta2)/2);
//	else if(phiA >= Pi-(thetaV3-thetaB) && phiA <= Pi+(thetaV4-thetaB)) Ao = 0;
//
//	else if(phiA >= Pi+(thetaV4-thetaB) && phiA <= Pi+(thetaV4-thetaB)+(rV2+rB*(1-cos(theta_1)))/R)	
//		Ao = rB*rB*(theta1-sin(2*theta1)/2)+rV2*rV2*(theta2-sin(2*theta2)/2);
//	else if(phiA >= Pi+thetaV4-thetaB+(rV2+rB*(1-cos(theta_1)))/R && phiA <= Pi+thetaV4+thetaB) 
//		Ao = rB*rB*(theta_1-sin(2*theta_1)/2)+Pi/2*rV2*rV2+2*R*rV2*(phiA-Pi-thetaV4+thetaB-(rV2+rB*(1-cos(theta_1)))/R);
//	else if(phiA >= Pi+thetaV4+thetaB && phiA <= Pi+thetaV4+thetaB+(rV2-rB*(1-cos(theta_1)))/R) 
//		Ao = rB*rB*(theta_1-sin(2*theta_1)/2)+Pi/2*rV2*rV2+2*R*rV2*(phiA-Pi-thetaV4+thetaB-(rV2+rB*(1-cos(theta_1)))/R)
//			-(rV2*rV2*(theta2-sin(2*theta2)/2)-rB*rB*(theta1-sin(2*theta1)/2));
//	else if(phiA >= Pi+thetaV4+thetaB+(rV2-rB*(1-cos(theta_1)))/R && phiA <= Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R)
//		Ao = 4*R*rV2*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2);
//	else if(phiA >= Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R && phiA <= Pi+thetaV4+thetaV5-thetaB)
//		Ao = 4*R*rV2*(thetaB-rB*(1-cos(theta_1))/R)+2*rB*rB*(theta_1-sin(2*theta_1)/2)
//			-2*R*rV2*(phiA-(Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R))+Pi/2*rV2*rV2-rB*rB*(theta_1-sin(2*theta_1)/2)
//			-(rV2*rV2*(theta2-sin(2*theta2)/2)-rB*rB*(theta1-sin(2*theta1)/2));
//	else if(phiA >= Pi+thetaV4+thetaV5-thetaB && phiA <= 2*Pi-(thetaV6-thetaB)-(rV2+rB*(1-cos(theta_1)))/R)
//		Ao = Pi/2*rV2*rV2+4*R*rV2*(thetaB-rB*(1-cos(theta_1))/R)+rB*rB*(theta_1-sin(2*theta_1)/2)-2*R*rV2*((rV2-rB*(1-cos(theta_1)))/R+phiA-(Pi+thetaV4+thetaV5-thetaB));
//	else if(phiA >= 2*Pi-(thetaV6-thetaB)-(rV2+rB*(1-cos(theta_1)))/R && phiA <= 2*Pi-(thetaV6-thetaB))	
//		Ao = rB*rB*(theta1-sin(2*theta1)/2)+rV2*rV2*(theta2-sin(2*theta2)/2);
//	else if(phiA >= 2*Pi-(thetaV6-thetaB) && phiA <= 2*Pi) Ao = 0;
//
//	//printf("%e %.10e \n", phi*180/Pi, Ao);
//	
//	return Ao;
//}

/* calculation of overlap angles */
void theta12(double phit)
{
	double rV, u, v, dudtheta1, dudtheta2, dvdtheta1, dvdtheta2, J, ea1, ea2;

	if(phit <= Pi) rV = rV1;
	else rV = rV2;

	theta_1 = asin(rV/rB);

	if(phit >= thetaV1-thetaB && phit <= thetaV1-thetaB+(rV+rB*(1-cos(theta_1)))/R){
		theta1 = 45*Pi/180;
		theta2 = 45*Pi/180;
		do{
			u = 2*R*sin((phit-(thetaV1-thetaB))/2)-rB*(1-cos(theta1))-rV*(1-cos(theta2));
			v = rB*sin(theta1)-rV*sin(theta2);
			dudtheta1 = -rB*sin(theta1);
			dudtheta2 = -rV*sin(theta2);
			dvdtheta1 = rB*cos(theta1);
			dvdtheta2 = -rV*cos(theta2);
			J = dudtheta1*dvdtheta2 - dudtheta2*dvdtheta1;
			theta1 = theta1 - (u*dvdtheta2-v*dudtheta2)/J;
			theta2 = theta2 - (v*dudtheta1-u*dvdtheta1)/J;
			ea1 = (u*dvdtheta2-v*dudtheta2)/J/theta1*100;
			ea2 = (v*dudtheta1-u*dvdtheta1)/J/theta2*100;
		} while(ea1 > eps && ea2 > eps);
		//fprintf(gp,"%e %.10e %.10e \n", phi*180/Pi, theta1, theta2);
	}
	else if(phit >= thetaV1+thetaB && phit <= thetaV1+thetaB+(rV-rB*(1-cos(theta_1)))/R){
		theta1 = 0.999;
		theta2 = 0.4;
		do{
			u = 2*R*sin((phit-(thetaV1+thetaB))/2)+rB*(1-theta1)-rV*(1-theta2);
			v = rB*rB*(1-theta1*theta1)-rV*rV*(1-theta2*theta2);
			dudtheta1 = -rB;
			dudtheta2 = rV;
			dvdtheta1 = -2*rB*rB*theta1;
			dvdtheta2 = 2*rV*rV*theta2;
			J = dudtheta1*dvdtheta2 - dudtheta2*dvdtheta1;
			theta1 = theta1 - (u*dvdtheta2-v*dudtheta2)/J;
			theta2 = theta2 - (v*dudtheta1-u*dvdtheta1)/J;
			ea1 = (u*dvdtheta2-v*dudtheta2)/J/theta1*100;
			ea2 = (v*dudtheta1-u*dvdtheta1)/J/theta2*100;
		} while(ea1 > eps && ea2 > eps);
		theta1 = acos(theta1);
		theta2 = acos(theta2);
		//fprintf(gp,"%e %.10e %.10e \n", phi*180/Pi, theta1, theta2);
	}
	else if(phit >= thetaV1+thetaV2-thetaB-(rV-rB*(1-cos(theta_1)))/R && phit <= thetaV1+thetaV2-thetaB){
		theta1 = 0.8;
		theta2 = 0.3;
		do{
			u = rV-rB*(1-cos(theta_1))-2*R*sin((phit-(thetaV1+thetaV2-thetaB-(rV-rB*(1-cos(theta_1)))/R))/2)+rB*(1-theta1)-rV*(1-theta2);
			v = rB*rB*(1-theta1*theta1)-rV*rV*(1-theta2*theta2);
			dudtheta1 = -rB;
			dudtheta2 = rV;
			dvdtheta1 = -2*rB*rB*theta1;
			dvdtheta2 = 2*rV*rV*theta2;
			J = dudtheta1*dvdtheta2 - dudtheta2*dvdtheta1;
			theta1 = theta1 - (u*dvdtheta2-v*dudtheta2)/J;
			theta2 = theta2 - (v*dudtheta1-u*dvdtheta1)/J;
			ea1 = (u*dvdtheta2-v*dudtheta2)/J/theta1*100;
			ea2 = (v*dudtheta1-u*dvdtheta1)/J/theta2*100;
		} while(ea1 > eps && ea2 > eps);
		theta1 = acos(theta1);
		theta2 = acos(theta2);
		//fprintf(gp,"%e %.10e %.10e \n", phi*180/Pi, theta1, theta2);
	}
	else if(phit >= Pi-(thetaV3-thetaB)-(rV+rB*(1-cos(theta_1)))/R && phit <= Pi-(thetaV3-thetaB)){
		theta1 = 20*Pi/180;
		theta2 = 80*Pi/180;
		do{
			u = rV+rB*(1-cos(theta_1))-2*R*sin((phit-(Pi-(thetaV3-thetaB)-(rV+rB*(1-cos(theta_1)))/R))/2)-rB*(1-cos(theta1))-rV*(1-cos(theta2));
			v = rB*sin(theta1)-rV*sin(theta2);
			dudtheta1 = -rB*sin(theta1);
			dudtheta2 = -rV*sin(theta2);
			dvdtheta1 = rB*cos(theta1);
			dvdtheta2 = -rV*cos(theta2);
			J = dudtheta1*dvdtheta2 - dudtheta2*dvdtheta1;
			theta1 = theta1 - (u*dvdtheta2-v*dudtheta2)/J;
			theta2 = theta2 - (v*dudtheta1-u*dvdtheta1)/J;
			//if(theta1>Pi/2 || theta2>Pi/2){
			//	theta1 = 
			ea1 = (u*dvdtheta2-v*dudtheta2)/J/theta1*100;
			ea2 = (v*dudtheta1-u*dvdtheta1)/J/theta2*100;
		}while(ea1 > eps && ea2 > eps);
		//fprintf(gp,"%e %.10e %.10e \n", phi*180/Pi, theta1, theta2);
	}
	else if(phit >= Pi+(thetaV4-thetaB) && phit <= Pi+(thetaV4-thetaB)+(rV2+rB*(1-cos(theta_1)))/R){
		theta1 = 45*Pi/180;
		theta2 = 45*Pi/180;
		do{
			u = 2*R*sin((phit-(Pi+thetaV4-thetaB))/2)-rB*(1-cos(theta1))-rV*(1-cos(theta2));
			v = rB*sin(theta1)-rV*sin(theta2);
			dudtheta1 = -rB*sin(theta1);
			dudtheta2 = -rV*sin(theta2);
			dvdtheta1 = rB*cos(theta1);
			dvdtheta2 = -rV*cos(theta2);
			J = dudtheta1*dvdtheta2 - dudtheta2*dvdtheta1;
			theta1 = theta1 - (u*dvdtheta2-v*dudtheta2)/J;
			theta2 = theta2 - (v*dudtheta1-u*dvdtheta1)/J;
			ea1 = (u*dvdtheta2-v*dudtheta2)/J/theta1*100;
			ea2 = (v*dudtheta1-u*dvdtheta1)/J/theta2*100;
		}while(ea1 > eps && ea2 > eps);
		//fprintf(gp,"%e %.10e %.10e \n", phi*180/Pi, theta1, theta2);
	}
	else if(phit >= Pi+thetaV4+thetaB && phit <= Pi+thetaV4+thetaB+(rV2-rB*(1-cos(theta_1)))/R){
		theta1 = 0.999;
		theta2 = 0.4;
		do{
			u = 2*R*sin((phit-(Pi+thetaV4+thetaB))/2)+rB*(1-theta1)-rV*(1-theta2);
			v = rB*rB*(1-theta1*theta1)-rV*rV*(1-theta2*theta2);
			dudtheta1 = -rB;
			dudtheta2 = rV;
			dvdtheta1 = -2*rB*rB*theta1;
			dvdtheta2 = 2*rV*rV*theta2;
			J = dudtheta1*dvdtheta2 - dudtheta2*dvdtheta1;
			theta1 = theta1 - (u*dvdtheta2-v*dudtheta2)/J;
			theta2 = theta2 - (v*dudtheta1-u*dvdtheta1)/J;
			ea1 = (u*dvdtheta2-v*dudtheta2)/J/theta1*100;
			ea2 = (v*dudtheta1-u*dvdtheta1)/J/theta2*100;
		} while(ea1 > eps && ea2 > eps);
		theta1 = acos(theta1);
		theta2 = acos(theta2);
		//fprintf(gp,"%e %.10e %.10e \n", phi*180/Pi, theta1, theta2);
	}
	else if(phit >= Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R && phit <= Pi+thetaV4+thetaV5-thetaB){
		theta1 = 0.8;
		theta2 = 0.3;
		do{
			u = rV-rB*(1-cos(theta_1))-2*R*sin((phit-(Pi+thetaV4+thetaV5-thetaB-(rV-rB*(1-cos(theta_1)))/R))/2)+rB*(1-theta1)-rV*(1-theta2);
			v = rB*rB*(1-theta1*theta1)-rV*rV*(1-theta2*theta2);
			dudtheta1 = -rB;
			dudtheta2 = rV;
			dvdtheta1 = -2*rB*rB*theta1;
			dvdtheta2 = 2*rV*rV*theta2;
			J = dudtheta1*dvdtheta2 - dudtheta2*dvdtheta1;
			theta1 = theta1 - (u*dvdtheta2-v*dudtheta2)/J;
			theta2 = theta2 - (v*dudtheta1-u*dvdtheta1)/J;
			ea1 = (u*dvdtheta2-v*dudtheta2)/J/theta1*100;
			ea2 = (v*dudtheta1-u*dvdtheta1)/J/theta2*100;
		} while(ea1 > eps && ea2 > eps);
		theta1 = acos(theta1);
		theta2 = acos(theta2);
		//fprintf(gp,"%e %.10e %.10e \n", phi*180/Pi, theta1, theta2);
	}
	else if(phit >= 2*Pi-(thetaV6-thetaB)-(rV2+rB*(1-cos(theta_1)))/R && phit <= 2*Pi-(thetaV6-thetaB)){
		theta1 = 45*Pi/180;
		theta2 = 65*Pi/180;
		do{
			u = rV+rB*(1-cos(theta_1))-2*R*sin((phit-(2*Pi-(thetaV6-thetaB)-(rV+rB*(1-cos(theta_1)))/R))/2)-rB*(1-cos(theta1))-rV*(1-cos(theta2));
			v = rB*sin(theta1)-rV*sin(theta2);
			dudtheta1 = -rB*sin(theta1);
			dudtheta2 = -rV*sin(theta2);
			dvdtheta1 = rB*cos(theta1);
			dvdtheta2 = -rV*cos(theta2);
			J = dudtheta1*dvdtheta2 - dudtheta2*dvdtheta1;
			theta1 = theta1 - (u*dvdtheta2-v*dudtheta2)/J;
			theta2 = theta2 - (v*dudtheta1-u*dvdtheta1)/J;
			ea1 = (u*dvdtheta2-v*dudtheta2)/J/theta1*100;
			ea2 = (v*dudtheta1-u*dvdtheta1)/J/theta2*100;
		}while(ea1 > eps && ea2 > eps);
		//fprintf(gp,"%e %.10e %.10e \n", phi*180/Pi, theta1, theta2);
	}
	else {
		theta1 = theta_1;
		theta2 = Pi/2;
	}
	//printf("%lf %.10e %.10e \n", phit*180/Pi, theta1, theta2);
}
/* Flow near roundings */
void integBpre()
{
	int i, j;
	double sum1, sum2, L;
	double theta_11, theta_12;

	theta_11 = asin(rV1/rB);
	theta_12 = asin(rV2/rB);

	//Precompress range
	sum1 = pow(hv,3)/12/eta/log(1+LB2/(R+rB));
	for(i=1; i<Ni; i++){
		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(i*Pi/Ni)/2)-rB;
		sum1 = sum1 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	}
	sum1 = sum1 + pow(hv,3)/12/eta/log(1+LB1/(R+rB));
	integB1 = Pi/Ni*2*sum1/2;
	//Outlet range
	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta_11)/2)-rB;
	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	for(j=int(Ni*(Pi/2-theta_11)/Pi+1.5); j<int(Ni*(Pi/2+theta_11)/Pi+0.5); j++){
		L = sqrt(((LB1-rB)*(LB1-rB)+(LB2-rB)*(LB2-rB))/2-((LB1-rB)*(LB1-rB)-(LB2-rB)*(LB2-rB))*cos(j*Pi/Ni)/2)-rB;
		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	}
	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta_11)/2)-rB;
	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	integB2 = integB1 - Pi/Ni*2*sum2/2;
	//Inlet range
	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta_12)/2)-rB;
	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	for(j=int(Ni*(Pi/2-theta_12)/Pi+1.5); j<int(Ni*(Pi/2+theta_12)/Pi+0.5); j++){
		L = sqrt(((LB1-rB)*(LB1-rB)+(LB2-rB)*(LB2-rB))/2-((LB1-rB)*(LB1-rB)-(LB2-rB)*(LB2-rB))*cos(j*Pi/Ni)/2)-rB;
		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	}
	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta_12)/2)-rB;
	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	integB3 = integB1 - Pi/Ni*2*sum2/2;
}
double integB(double phii)
{
	int i, j;
	double sum1, sum2, sum3, L;
	double integB=0.;

	sum1 = pow(hv,3)/12/eta/log(1+LB2/(R+rB));
	for(i=1; i<Ni; i++){
		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(i*Pi/Ni)/2)-rB;
		sum1 = sum1 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	}
	sum1 = sum1 + pow(hv,3)/12/eta/log(1+LB1/(R+rB));
	integB1 = Pi/Ni*2*sum1/2;

	//theta12(phii);

	//if(phii >= 0 && phii <= thetaV1-thetaB) integB = integB1;
	//else if(phii >= thetaV1-thetaB && phii <= thetaV1-thetaB+(rV1+rB*(1-cos(theta_1)))/R){
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta1)/2)-rB;
	//	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta1)/Pi+1.5); j<int(Ni*(Pi/2+theta1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta1)/2)-rB;
	//	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	integB = integB1 - Pi/Ni*sum2/2;
	//}	
	//else if(phii >= thetaV1-thetaB+(rV1+rB*(1-cos(theta_1)))/R && phii <= thetaV1+thetaB){
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta_1)/2)-rB;
	//	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta_1)/Pi+1.5); j<int(Ni*(Pi/2+theta_1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta_1)/2)-rB;
	//	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	integB = integB1 - Pi/Ni*sum2/2;
	//}			
	//else if(phii >= thetaV1+thetaB && phii <= thetaV1+thetaB+(rV1-rB*(1-cos(theta_1)))/R){
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta_1)/2)-rB;
	//	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta_1)/Pi+1.5); j<int(Ni*(Pi/2+theta_1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta_1)/2)-rB;
	//	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta1)/2)-rB;
	//	sum3 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta1)/Pi+1.5); j<int(Ni*(Pi/2+theta1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum3 = sum3 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta1)/2)-rB;
	//	sum3 = sum3 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	integB = integB1 - Pi/Ni*(sum2+sum3)/2;
	//}
	//else if(phii >= thetaV1+thetaB+(rV1-rB*(1-cos(theta_1)))/R && phii <= thetaV1+thetaV2-thetaB-(rV1-rB*(1-cos(theta_1)))/R)
	//	integB = integB2;	
	//else if(phii >= thetaV1+thetaV2-thetaB-(rV1-rB*(1-cos(theta_1)))/R && phii <= thetaV1+thetaV2-thetaB){
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta_1)/2)-rB;
	//	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta_1)/Pi+1.5); j<int(Ni*(Pi/2+theta_1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta_1)/2)-rB;
	//	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta1)/2)-rB;
	//	sum3 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta1)/Pi+1.5); j<int(Ni*(Pi/2+theta1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum3 = sum3 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta1)/2)-rB;
	//	sum3 = sum3 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	integB = integB1 - Pi/Ni*(sum2+sum3)/2;
	//}
	//else if(phii >= thetaV1+thetaV2-thetaB && phii <= Pi-(thetaV3-thetaB)-(rV1+rB*(1-cos(theta_1)))/R){
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta_1)/2)-rB;
	//	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta_1)/Pi+1.5); j<int(Ni*(Pi/2+theta_1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta_1)/2)-rB;
	//	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	integB = integB1 - Pi/Ni*sum2/2;
	//}	

	//else if(phii >= Pi-(thetaV3-thetaB)-(rV1+rB*(1-cos(theta_1)))/R && phii <= Pi-(thetaV3-thetaB)){	
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta1)/2)-rB;
	//	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta1)/Pi+1.5); j<int(Ni*(Pi/2+theta1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta1)/2)-rB;
	//	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	integB = integB1 - Pi/Ni*sum2/2;
	//}	
	//else if(phii >= Pi-(thetaV3-thetaB) && phii <= Pi+(thetaV4-thetaB))	integB = integB1;
	//else if(phii >= Pi+(thetaV4-thetaB) && phii <= Pi+(thetaV4-thetaB)+(rV2+rB*(1-cos(theta_1)))/R){ 
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta1)/2)-rB;
	//	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta1)/Pi+1.5); j<int(Ni*(Pi/2+theta1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta1)/2)-rB;
	//	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	integB = integB1 - Pi/Ni*sum2/2;
	//}	
	//else if(phii >= Pi+thetaV4-thetaB+(rV2+rB*(1-cos(theta_1)))/R && phii <= Pi+thetaV4+thetaB){
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta_1)/2)-rB;
	//	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta_1)/Pi+1.5); j<int(Ni*(Pi/2+theta_1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta_1)/2)-rB;
	//	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	integB = integB1 - Pi/Ni*sum2/2;
	//} 
	//else if(phii >= Pi+thetaV4+thetaB && phii <= Pi+thetaV4+thetaB+(rV2-rB*(1-cos(theta_1)))/R){
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta_1)/2)-rB;
	//	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta_1)/Pi+1.5); j<int(Ni*(Pi/2+theta_1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta_1)/2)-rB;
	//	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta1)/2)-rB;
	//	sum3 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta1)/Pi+1.5); j<int(Ni*(Pi/2+theta1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum3 = sum3 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta1)/2)-rB;
	//	sum3 = sum3 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	integB = integB1 - Pi/Ni*(sum2+sum3)/2;
	//} 
	//else if(phii >= Pi+thetaV4+thetaB+(rV2-rB*(1-cos(theta_1)))/R && phii <= Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R)
	//	integB = integB3;
	//else if(phii >= Pi+thetaV4+thetaV5-thetaB-(rV2-rB*(1-cos(theta_1)))/R && phii <= Pi+thetaV4+thetaV5-thetaB){
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta_1)/2)-rB;
	//	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta_1)/Pi+1.5); j<int(Ni*(Pi/2+theta_1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta_1)/2)-rB;
	//	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta1)/2)-rB;
	//	sum3 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta1)/Pi+1.5); j<int(Ni*(Pi/2+theta1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum3 = sum3 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta1)/2)-rB;
	//	sum3 = sum3 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	integB = integB1 - Pi/Ni*(sum2+sum3)/2;
	//}
	//else if(phii >= Pi+thetaV4+thetaV5-thetaB && phii <= 2*Pi-(thetaV6-thetaB)-(rV2+rB*(1-cos(theta_1)))/R){
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta_1)/2)-rB;
	//	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta_1)/Pi+1.5); j<int(Ni*(Pi/2+theta_1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta_1)/2)-rB;
	//	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	integB = integB1 - Pi/Ni*sum2/2;
	//}
	//else if(phii >= 2*Pi-(thetaV6-thetaB)-(rV2+rB*(1-cos(theta_1)))/R && phii <= 2*Pi-(thetaV6-thetaB)){	
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2-theta1)/2)-rB;
	//	sum2 = pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	for(j=int(Ni*(Pi/2-theta1)/Pi+1.5); j<int(Ni*(Pi/2+theta1)/Pi+0.5); j++){
	//		L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(j*Pi/Ni)/2)-rB;
	//		sum2 = sum2 + 2*pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	}
	//	L = sqrt(((LB1+rB)*(LB1+rB)+(LB2+rB)*(LB2+rB))/2-((LB1+rB)*(LB1+rB)-(LB2+rB)*(LB2+rB))*cos(Pi/2+theta1)/2)-rB;
	//	sum2 = sum2 + pow(hv,3)/12/eta/log(1+L/(R+rB));
	//	integB = integB1 - Pi/Ni*sum2/2;
	//}
	//else if(phii >= 2*Pi-(thetaV6-thetaB) && phii <= 2*Pi) integB = integB1;

	return integB;
}