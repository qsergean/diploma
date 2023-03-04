
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "riem.h"

#define e(d,u,p)  p/(g-1.)+0.5*d*u*u

void riemi_H2(double* rqp_io, int *alarm)
{

	double om = 0.;

	double c[2],rc[2],q[3][2],du,ff,df;
	double ss[2][2],rr[2],cc[2],rqp[3][2];
	//definition for ff_df
	double rc_1;
	double uo[3],f[3][2];

	int i;

	double e_pmin,u_vc,p0,p1,pp,uu,ai,tmp1,tmp2;
        double cz,rz,rcz,pz,uz,sl,sr,ql,qr;        
        int lom=0;
        *alarm=0;

        rqp[0][0]=rqp_io[0]; rqp[1][0]=rqp_io[1]; rqp[2][0]=rqp_io[2];
        rqp[0][1]=rqp_io[3]; rqp[1][1]=rqp_io[4]; rqp[2][1]=rqp_io[5];

	// i==0 //
	{
		c[0]=sqrt(g*rqp_io[2]/rqp_io[0]);
		rc[0]=rqp_io[0]*c[0];
	}
	// i==1 //
	{
		c[1]=sqrt(g*rqp_io[5]/rqp_io[3]);
		rc[1]=rqp_io[3]*c[1];
	}
	e_pmin=min(rqp_io[2],rqp_io[5])*ee;

	/* check existance */
	u_vc=-c[0]*gc-c[1]*gc;
	du=rqp_io[1]-rqp_io[4];
	if(du <= u_vc)
	{
	printf("\n vacuem - no solution \n");
	printf(" Left  %e  %e  %e  \n",rqp_io[0],rqp_io[1],rqp_io[2]);
	printf(" Right %e  %e  %e  \n",rqp_io[3],rqp_io[4],rqp_io[5]);
        *alarm=77;
        return;
	}

       cz=0.5*(c[0]+c[1]);
       rz=0.5*(rqp[0][0]+rqp[0][1]);
       rcz=cz*rz;
       pz=0.5*(rqp[2][0]+rqp[2][1])+0.5*(rqp[1][0]-rqp[1][1])*rcz;
       uz=0.5*(rqp[1][0]+rqp[1][1])+0.5*(rqp[2][0]-rqp[2][1])/rcz;

       if(rqp[2][0] >= pz) ql=1.;
       else ql=sqrt(1.+ga*(pz/rqp[2][0]-1.));

       if(rqp[2][1] >= pz) qr=1.;
       else qr=sqrt(1.+ga*(pz/rqp[2][1]-1.));

       sl=rqp[1][0]-c[0]*ql;  sr=rqp[1][1]+c[1]*qr;
//       sl=rqp[1][0]-c[0];  sr=rqp[1][1]+c[1];
        

        rqp_io[9]=sl;
	rqp_io[11]=sr;

        for(i=0; i < 2; i++)
         { 
          q[0][i]=rqp[0][i];
          q[1][i]=rqp[0][i]*rqp[1][i];
          q[2][i]=e(rqp[0][i],rqp[1][i],rqp[2][i]);
          f[0][i]=rqp[0][i]*rqp[1][i];
          f[1][i]=f[0][i]*rqp[1][i]+rqp[2][i];
          f[2][i]=(e(rqp[0][i],rqp[1][i],rqp[2][i])+rqp[2][i])*rqp[1][i];
         } 

        for(i=0; i < 3; i++)
         uo[i]=(sr*q[i][1]-sl*q[i][0]+f[i][0]-f[i][1])/(sr-sl);
//         uo[i]=(sr*q[i][1]-sl*q[i][0]+f[i][1]-f[i][0])/(sr-sl);

        rqp_io[10]=uo[1]/uo[0];

	/* 1 */
	if( om <= sl )
	{
	rqp_io[6]=rqp_io[0]; rqp_io[7]=rqp_io[1]; rqp_io[8]=rqp_io[2]; goto rtn_o; 
	}

	/* 2 */
	if( ( om > sl ) && ( om <= sr) )
	{
         rqp_io[6]=uo[0]; 
         rqp_io[7]=uo[1]/uo[0];
         rqp_io[8]=(g-1.)*(uo[2]-0.5*uo[1]*rqp_io[7]);
	goto rtn_o;
	}
         
	/* 3 */
	if( om > sr )	{ rqp_io[6]=rqp_io[3]; rqp_io[7]=rqp_io[4]; rqp_io[8]=rqp_io[5]; goto rtn_o; }

	rtn_o:

	if(rqp_io[8] <= e_pmin)
	{
		printf("\n Final pressure <= e_pmin - no solution \n");
		printf(" Left  %e  %e  %e  \n",rqp_io[0],rqp_io[1],rqp_io[2]);
		printf(" Right %e  %e  %e  \n",rqp_io[3],rqp_io[4],rqp_io[5]);
                printf("Waves %e  %e  \n",sl,sr);
                printf("Midle state %e  %e  %e  \n",uo[0],uo[1],uo[2]);
                printf("Solution %e  %e  %e \n",rqp_io[6],rqp_io[7],rqp_io[8]);
                *alarm=79;
                return;
	}
}







