/*
CONTENTS
void inp_mesh(void);
void inp_gdf(void);
void deriv_lr(void);
void deriv_du(void);
void step(void);
double flow_lr(void);
double flow_du(void);
void new_gdf(void);
*/

//void turb_lr(void);
//void turb_du(void);
//void mu_mol_func_lr(void);
//void mu_mol_func_du(void);


#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "func.h"

#define IN(i,j)   ( (i)*(JT+1)+(j) )
#define IC(i,j)   ( (i)*JT+(j) )
#define Ilr(i,j)  ( (i)*JT+(j) )
#define Idu(i,j)  ( (i)*(JT+1)+(j) )


#define  GAM  1.4
#define PI 3.141592653 


double min_mod(double a, double b);


double *gdf[6],*x,*y,*fllr[6],*fldu[6],*xc,*yc,
       *nlr[3],*tlr[2],*ndu[3],*tdu[2],*vol;
double *dxflr[6],*dyflr[6],*dxfdu[6],*dyfdu[6];
double t,dt,t_max,coord,turb,Mach,Re_ref,u_ref,mu_ref;
double lxmin,lymin,mu_max,H_Z;
int n_stp,max_stp,inp_par,IT,JT;

double *wlr[2],*wdu[2],*gdf_[6];
double tau;
int i_qgd,i_appr,i_Rey;
double RES;


double *mu_turb_lr, *mu_turb_du;
int Flag_turb;

double eforstep=1.e-13;


//!!! For UY equations 
double *fllr5,*fldu5,*dxflr5[5],*dyflr5[5],*dxfdu5[5],*dyfdu5[5],*rhs5_1;




//*****************************************************

int main(void)
{
 FILE *fpr;
 int lout,nout;

 fpr=fopen("par_i.dat","r");
 fscanf(fpr,"%d %d",&IT,&JT);
 fscanf(fpr,"%lg",&t_max);
 fscanf(fpr,"%d %d %d",&max_stp,&nout,&inp_par);
 fscanf(fpr,"%lg %lg",&Mach,&Re_ref);
 fscanf(fpr,"%lg %d",&coord,&Flag_turb);
 fscanf(fpr,"%lg",&H_Z);
 fscanf(fpr,"%d",&i_qgd);
 fscanf(fpr,"%d",&i_appr);
 fscanf(fpr,"%d",&i_Rey);


 fclose(fpr);

 t=0.;
 dt=1.e-17;

 inp_mesh();
 inp_gdf();
 n_stp=1;
 lout=1;

 while( t < t_max )
   {
   double lcell,ki,bi;

   printf(" N_step=%d t=%e dt=%e RES/dt=%e \n",n_stp,t,dt,RES/dt);
   step();

   if(lout > nout)
     {
      out();
      out_tec();
      lout=1;
     }


   t+=dt;
   n_stp++;
   lout++;

   if(n_stp > max_stp) break;
   }

 out();
 out_tec();
}

//*****************************************************
void step(void)
{
 int i,j,k;
 double dt_lr,dt_du,lcell,ki,bi;
 FILE *fpr;


// BOUNDARY CONDITION:  ghost cells
// right out

   for(j=0; j < JT; j++)
    {
     gdf[0][IC(IT-2,j)] =gdf[0][IC(IT-1,j)] =gdf[0][IC(IT-3,j)];
     gdf[1][IC(IT-2,j)] =gdf[1][IC(IT-1,j)] =gdf[1][IC(IT-3,j)];
     gdf[2][IC(IT-2,j)] =gdf[2][IC(IT-1,j)] =gdf[2][IC(IT-3,j)];
     gdf[3][IC(IT-2,j)] =gdf[3][IC(IT-1,j)] =gdf[3][IC(IT-3,j)];
     gdf[4][IC(IT-2,j)] =gdf[4][IC(IT-1,j)] =gdf[4][IC(IT-3,j)];
    }

// left in

// down  wall
   for(i=2; i < IT-2; i++)
    {
     gdf[0][IC(i,0)]=gdf[0][IC(i,1)]= gdf[0][IC(i,2)];
     gdf[1][IC(i,0)]=gdf[1][IC(i,1)]=-gdf[1][IC(i,2)];
     gdf[2][IC(i,0)]=gdf[2][IC(i,1)]=-gdf[2][IC(i,2)];
     gdf[3][IC(i,0)]=gdf[3][IC(i,1)]= gdf[3][IC(i,2)];
     gdf[4][IC(i,0)]=gdf[4][IC(i,1)]= gdf[4][IC(i,2)];
    }



// up  in
/*
    for(i=2; i < IT-2; i++)
    {
     gdf[0][IC(i,JT-1)]=gdf[0][IC(i,JT-2)]= gdf[0][IC(i,JT-3)];
     gdf[1][IC(i,JT-1)]=gdf[1][IC(i,JT-2)]= gdf[1][IC(i,JT-3)];
     gdf[2][IC(i,JT-1)]=gdf[2][IC(i,JT-2)]=-gdf[2][IC(i,JT-3)];
     gdf[3][IC(i,JT-1)]=gdf[3][IC(i,JT-2)]= gdf[3][IC(i,JT-3)];
     gdf[4][IC(i,JT-1)]=gdf[4][IC(i,JT-2)]= gdf[4][IC(i,JT-3)];
     gdf[5][IC(i,JT-1)]=gdf[5][IC(i,JT-2)]= gdf[5][IC(i,JT-3)];
    }
*/
 deriv_lr();
 deriv_du();
// deriv_lr_5();
// deriv_du_5();

/*
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   rhs5_1[IC(i,j)]=deriv_t(i,j);
//    rhs5_1[IC(i,j)]=0.;
*/

// if(Flag_turb == 1)  turb_lr();

 dt_lr=flow_lr();

// if(Flag_turb == 1)   turb_du();

 dt_du=flow_du();


 new_gdf();
// smoof();


 lcell=min(lxmin,lymin);
 ki=min(dt_lr,dt_du);
 bi=0.25*lcell*lcell/mu_max;
 dt=min(bi,ki);
 dt=H_Z*dt;
}

//*****************************************************

void inp_mesh(void)
{
 int i,j;
 double dx,dy;
 FILE *fpr;

 fpr=fopen("mesh.dat","r");
// inp coord 
 fscanf(fpr,"%d %d",&IT,&JT);
 memo_on();

//for(i=2; i <= IT-2; i++)
//  for(j=2; j <= JT-2; j++)

  for(j=2; j <= JT-2; j++)
   for(i=2; i <= IT-2; i++)

     fscanf(fpr,"%lg %lg",&x[IN(i,j)],&y[IN(i,j)]);

 lxmin=9.9e+9; lymin=9.9e+9;

// norms left_right
 for(i=2; i <= IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    dx=x[IN(i,j+1)]-x[IN(i,j)];
    dy=y[IN(i,j+1)]-y[IN(i,j)];
    nlr[0][Ilr(i,j)]=dy;
    nlr[1][Ilr(i,j)]=-dx;
    nlr[2][Ilr(i,j)]=sqrt(dx*dx+dy*dy);
    lxmin=min(lxmin,nlr[2][Ilr(i,j)]);
    tlr[0][Ilr(i,j)]=dx;
    tlr[1][Ilr(i,j)]=dy;
   }
// norms down_up
 for(i=2; i < IT-2; i++)
  for(j=2; j <= JT-2; j++)
   {
    dx=x[IN(i+1,j)]-x[IN(i,j)];
    dy=y[IN(i+1,j)]-y[IN(i,j)];
    ndu[0][Idu(i,j)]=-dy;
    ndu[1][Idu(i,j)]=dx;
    ndu[2][Idu(i,j)]=sqrt(dx*dx+dy*dy);
    lymin=min(lymin,ndu[2][Idu(i,j)]);
    tdu[0][Idu(i,j)]=dx;
    tdu[1][Idu(i,j)]=dy;
   }
// volumes
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    dx=x[IN(i+1,j+1)]-x[IN(i,j)];
    dy=y[IN(i,j+1)]-y[IN(i+1,j)];
    vol[IC(i,j)]=0.5*dx*dy;
    dx=x[IN(i,j+1)]-x[IN(i+1,j)];
    dy=y[IN(i+1,j+1)]-y[IN(i,j)];
    vol[IC(i,j)]=vol[IC(i,j)]-0.5*dx*dy;
    xc[IC(i,j)]=0.25*(x[IN(i,j)]+x[IN(i+1,j)]+x[IN(i+1,j+1)]+x[IN(i,j+1)]);
    yc[IC(i,j)]=0.25*(y[IN(i,j)]+y[IN(i+1,j)]+y[IN(i+1,j+1)]+y[IN(i,j+1)]);
   }
for( j=0; j < JT; j++)
{
 vol[IC(0,j)]=vol[IC(1,j)]=vol[IC(2,j)];
 vol[IC(IT-1,j)]=vol[IC(IT-2,j)]=vol[IC(IT-3,j)];
}

for( i=0; i < IT; i++)
{
 vol[IC(i,0)]=vol[IC(i,1)]=vol[IC(i,2)];
 vol[IC(i,JT-1)]=vol[IC(i,JT-2)]=vol[IC(i,JT-3)];
}

for(j=0; j < JT; j++)
  {
  xc[IC(IT-2,j)]=2.*xc[IC(IT-3,j)]-xc[IC(IT-4,j)];
  yc[IC(IT-2,j)]=2.*yc[IC(IT-3,j)]-yc[IC(IT-4,j)];
  xc[IC(IT-1,j)]=2.*xc[IC(IT-2,j)]-xc[IC(IT-3,j)];
  yc[IC(IT-1,j)]=2.*yc[IC(IT-2,j)]-yc[IC(IT-3,j)];

  xc[IC(1,j)]=2.*xc[IC(2,j)]-xc[IC(3,j)];
  yc[IC(1,j)]=2.*yc[IC(2,j)]-yc[IC(3,j)];
  xc[IC(0,j)]=2.*xc[IC(1,j)]-xc[IC(2,j)];
  yc[IC(0,j)]=2.*yc[IC(1,j)]-yc[IC(2,j)];
  }

for(i=0; i < IT; i++)
  {
  xc[IC(i,JT-2)]=2.*xc[IC(i,JT-3)]-xc[IC(i,JT-4)];
  yc[IC(i,JT-2)]=2.*yc[IC(i,JT-3)]-yc[IC(i,JT-4)];
  xc[IC(i,JT-1)]=2.*xc[IC(i,JT-2)]-xc[IC(i,JT-3)];
  yc[IC(i,JT-1)]=2.*yc[IC(i,JT-2)]-yc[IC(i,JT-3)];
                              
  xc[IC(i,1)]=2.*xc[IC(i,2)]-xc[IC(i,3)];
  yc[IC(i,1)]=2.*yc[IC(i,2)]-yc[IC(i,3)];
  xc[IC(i,0)]=2.*xc[IC(i,1)]-xc[IC(i,2)];
  yc[IC(i,0)]=2.*yc[IC(i,1)]-yc[IC(i,2)];
  }

 fclose(fpr);
}

//*****************************************************
void inp_gdf(void)
{
 int i,j,k;
 double d,u,v,p,e;
 double L_ref;


 FILE *fpr;

 d=1.; p=1.;  v=0.;
 u_ref=u=Mach*sqrt(GAM*p/d);
// !!!!!!!!!!!!!!!!
 L_ref=1.;
 mu_ref=u*L_ref*d/Re_ref;
 for(i=0; i < IT; i++)
  for(j=0; j < JT; j++)
   {
    k=IC(i,j);
    gdf[0][k]=d;
    gdf[1][k]=u;
    gdf[2][k]=v;
    gdf[3][k]=p;
    gdf[4][k]=p/d;
    gdf[5][k]=0.;
   }

 if(inp_par == 1)
  {
   fpr=fopen("reg_i.dat","rb");
   fread(&IT,sizeof(int),1,fpr);
   fread(&JT,sizeof(int),1,fpr);
   fread(&t,sizeof(double),1,fpr);
   fread(&dt,sizeof(double),1,fpr);
   for(i=2; i < IT-2; i++)
    for(j=2; j < JT-2; j++)
     for(k=0; k < 6; k++)
      fread(&gdf[k][IC(i,j)],sizeof(double),1,fpr);
   fclose(fpr);  
  }

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// BOUNDARY CONDITION:  ghost cells
// right out
   for(j=0; j < JT; j++)
    {
     gdf[0][IC(IT-2,j)] =gdf[0][IC(IT-1,j)] =gdf[0][IC(IT-3,j)];
     gdf[1][IC(IT-2,j)] =gdf[1][IC(IT-1,j)] =gdf[1][IC(IT-3,j)];
     gdf[2][IC(IT-2,j)] =gdf[2][IC(IT-1,j)] =gdf[2][IC(IT-3,j)];
     gdf[3][IC(IT-2,j)] =gdf[3][IC(IT-1,j)] =gdf[3][IC(IT-3,j)];
     gdf[4][IC(IT-2,j)] =gdf[4][IC(IT-1,j)] =gdf[4][IC(IT-3,j)];
    }


// left in

// down  wall
   for(i=2; i < IT-2; i++)
    {
     gdf[0][IC(i,0)]=gdf[0][IC(i,1)]= gdf[0][IC(i,2)];
     gdf[1][IC(i,0)]=gdf[1][IC(i,1)]=-gdf[1][IC(i,2)];
     gdf[2][IC(i,0)]=gdf[2][IC(i,1)]=-gdf[2][IC(i,2)];
     gdf[3][IC(i,0)]=gdf[3][IC(i,1)]= gdf[3][IC(i,2)];
     gdf[4][IC(i,0)]=gdf[4][IC(i,1)]= gdf[4][IC(i,2)];
     gdf[5][IC(i,0)]=gdf[5][IC(i,1)]=-gdf[5][IC(i,2)];
    }

// up  in
/*
   for(i=2; i < IT-2; i++)
    {
     gdf[0][IC(i,JT-1)]=gdf[0][IC(i,JT-2)]= gdf[0][IC(i,JT-3)];
     gdf[1][IC(i,JT-1)]=gdf[1][IC(i,JT-2)]= gdf[1][IC(i,JT-3)];
     gdf[2][IC(i,JT-1)]=gdf[2][IC(i,JT-2)]=-gdf[2][IC(i,JT-3)];
     gdf[3][IC(i,JT-1)]=gdf[3][IC(i,JT-2)]= gdf[3][IC(i,JT-3)];
     gdf[4][IC(i,JT-1)]=gdf[4][IC(i,JT-2)]= gdf[4][IC(i,JT-3)];
     gdf[5][IC(i,JT-1)]=gdf[5][IC(i,JT-2)]= gdf[5][IC(i,JT-3)];
    }
*/

// !!!!!!!!!!!!!!!!Tempo 
   for(i=0; i < IT; i++)
    for(j=0; j < JT; j++)                            
     for(k=0; k < 6; k++)
      gdf_[k][IC(i,j)]=gdf[k][IC(i,j)]; 

}

//*****************************************************
void deriv_lr(void)
{
 int i,j,k;
 double dx1,dy1,dx2,dy2,df1,df2;
 double DDx,DDy,DD;
 double xc1,yc1,xc2,yc2,f3,f4;

 double mum,c,gdb[6];


   for(i=2; i <= IT-2; i++)
    for(j=2; j < JT-2; j++)
     {
     xc1=xc[IC(i-1,j)];  
     yc1=yc[IC(i-1,j)];  
     xc2=xc[IC(i,j)];    
     yc2=yc[IC(i,j)];    
     dx1=xc2-xc1;
     dy1=yc2-yc1;
     dx2=x[IN(i,j+1)]-x[IN(i,j)];
     dy2=y[IN(i,j+1)]-y[IN(i,j)];
     DD=dx1*dy2-dx2*dy1;
     for(k=0; k < 6; k++)
      {
      df1=gdf[k][IC(i,j)]-gdf[k][IC(i-1,j)];
      f3=0.25*(gdf[k][IC(i-1,j-1)]+gdf[k][IC(i,j-1)]+
               gdf[k][IC(i,j)]+gdf[k][IC(i-1,j)]); 
      f4=0.25*(gdf[k][IC(i-1,j)]+gdf[k][IC(i,j)]+
               gdf[k][IC(i,j+1)]+gdf[k][IC(i-1,j+1)]); 

      df2=f4-f3;
      DDx=df1*dy2-df2*dy1;
      DDy=df2*dx1-df1*dx2;

      dxflr[k][Ilr(i,j)]=DDx/DD;
      dyflr[k][Ilr(i,j)]=DDy/DD;  
      gdb[k]=0.5*(gdf[k][IC(i,j)]+gdf[k][IC(i-1,j)]);
      }
// QGD 
mum=mu_ref*pow((gdb[3]/gdb[0]),0.76);
c=sqrt(GAM*gdb[3]/gdb[0]);
if(i_qgd) tau=mum/gdb[3];
//if(i_qgd) tau=mum/gdb[3]+sqrt(lxmin*lymin)/c;
else tau=0.;
/*     
      wlr[0][Ilr(i,j)]=0.5*tau*( gdb[1]*gdb[1]*dxflr[0][Ilr(i,j)]+
                            2.*gdb[0]*gdb[1]*dxflr[1][Ilr(i,j)]+
                            dxflr[3][Ilr(i,j)] + gdb[1]*gdb[2]*dyflr[0][Ilr(i,j)]+
                            gdb[0]*gdb[2]*dyflr[1][Ilr(i,j)]+
                            gdb[0]*gdb[1]*dyflr[2][Ilr(i,j)])/gdb[0];

      wlr[1][Ilr(i,j)]=0.5*tau*( gdb[1]*gdb[2]*dxflr[0][Ilr(i,j)]+
                            gdb[0]*gdb[2]*dxflr[1][Ilr(i,j)]+
                            gdb[0]*gdb[1]*dxflr[2][Ilr(i,j)]+
                            gdb[2]*gdb[2]*dyflr[0][Ilr(i,j)]+
                            2.*gdb[0]*gdb[2]*dyflr[2][Ilr(i,j)]+
                            dyflr[3][Ilr(i,j)] )/gdb[0];
*/
wlr[0][Ilr(i,j)]=0.;
wlr[1][Ilr(i,j)]=0.;

     }
}
//*****************************************************
void deriv_du(void)
{
 int i,j,k;
 double dx1,dy1,dx2,dy2,df1,df2;
 double DDx,DDy,DD;
 double xc3,yc3,xc4,yc4,f1,f2;

 double mum,c,gdb[6];

   for(i=2; i < IT-2; i++)
    for(j=2; j <= JT-2; j++)
     {
     xc3=xc[IC(i,j-1)];
     yc3=yc[IC(i,j-1)];
     xc4=xc[IC(i,j)]  ;
     yc4=yc[IC(i,j)]  ;
     dx2=xc4-xc3;
     dy2=yc4-yc3;
     dx1=x[IN(i+1,j)]-x[IN(i,j)];
     dy1=y[IN(i+1,j)]-y[IN(i,j)];
     DD=dx1*dy2-dx2*dy1;
     for(k=0; k < 6; k++)
      {
      df2=gdf[k][IC(i,j)]-gdf[k][IC(i,j-1)];
      f1=0.25*(gdf[k][IC(i-1,j-1)]+gdf[k][IC(i,j-1)]+
               gdf[k][IC(i,j)]+gdf[k][IC(i-1,j)]); 
      f2=0.25*(gdf[k][IC(i,j-1)]+gdf[k][IC(i+1,j-1)]+
               gdf[k][IC(i+1,j)]+gdf[k][IC(i,j)]);    
      df1=f2-f1;
      DDx=df1*dy2-df2*dy1;
      DDy=df2*dx1-df1*dx2;
      dxfdu[k][Idu(i,j)]=DDx/DD;
      dyfdu[k][Idu(i,j)]=DDy/DD;  
      gdb[k]=0.5*(gdf[k][IC(i,j)]+gdf[k][IC(i,j-1)]);
      }
// QGD 
mum=mu_ref*pow((gdb[3]/gdb[0]),0.76);
c=sqrt(GAM*gdb[3]/gdb[0]);
if(i_qgd) tau=mum/gdb[3];
//if(i_qgd) tau=mum/gdb[3]+sqrt(lxmin*lymin)/c;
else tau=0.;
/*      
      wdu[0][Idu(i,j)]=0.5*tau*( gdb[1]*gdb[1]*dxfdu[0][Idu(i,j)]+
                            2.*gdb[0]*gdb[1]*dxfdu[1][Idu(i,j)]+
                            dxfdu[3][Idu(i,j)] + gdb[1]*gdb[2]*dyfdu[0][Idu(i,j)]+
                            gdb[0]*gdb[2]*dyfdu[1][Idu(i,j)]+
                            gdb[0]*gdb[1]*dyfdu[2][Idu(i,j)])/gdb[0];

      wdu[1][Idu(i,j)]=0.5*tau*( gdb[1]*gdb[2]*dxfdu[0][Idu(i,j)]+
                            gdb[0]*gdb[2]*dxfdu[1][Idu(i,j)]+
                            gdb[0]*gdb[1]*dxfdu[2][Idu(i,j)]+
                            gdb[2]*gdb[2]*dyfdu[0][Idu(i,j)]+
                            2.*gdb[0]*gdb[2]*dyfdu[2][Idu(i,j)]+
                            dyfdu[3][Idu(i,j)] )/gdb[0];
*/
wdu[0][Idu(i,j)]=0.;
wdu[1][Idu(i,j)]=0.;


     }
}

//*****************************************************

double flow_lr(void)
{
 int i,j,k,l,lom,alarm;
 double dx,dy,db,ub,vb,pb,e,NV1,TV1,NV2,TV2;
 double t11,t12,t21,t22,q1,q2;
 double mum,mut,muef,kapm,kapt,kapef;
// static double dfx[6],dfy[6],ff[6];
 double dfx[6],dfy[6],ff[6];
 double xc1,yc1,xc2,yc2;
 double par[32],rqp[3][2],r_o,u_o,v_o,p_o,om,s1,s2,s3;
 double dksi1,dksi2;
 double delf,delf1,delf2,s;
 double w1,w2;
 double dt_lr,dt1,dt2,dt_max;
 double Sb,Sl,Sr;

int n;
double wzx,wzy,rz,dv; 

double pr1,pr2;

 lom=0;  om=0.; alarm=0;
 mu_max=0.;
 dt_lr=1.e17;
               
 for(i=2; i <= IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
// left  values for riem
	for( k=0; k < 6; k++)
	{

	dx = xc[IC(i-1,j)] - xc[IC(i - 2,j)];
	dy = yc[IC(i-1,j)] - yc[IC(i - 2,j)];
	dksi1 = sqrt(dx*dx+dy*dy);
	dx = xc[IC(i,j)] - xc[IC(i - 1,j)];
	dy = yc[IC(i,j)] - yc[IC(i - 1,j)];
	dksi2 = sqrt(dx*dx+dy*dy);
	delf1 = (gdf[k][IC(i-1,j)] - gdf[k][IC(i-2,j)]) / dksi1;
	delf2 = (gdf[k][IC(i,j)] - gdf[k][IC(i-1,j)]) / dksi2;
	if(i_appr)
        {
	delf = min_mod(delf1,delf2);
	ff[k]=gdf[k][IC(i-1,j)]+0.5*dksi2*delf;
        }
        else ff[k]=gdf[k][IC(i-1,j)]; 
	}
Sl=ff[5];

     NV1=(nlr[0][Ilr(i,j)]*ff[1]+nlr[1][Ilr(i,j)]*ff[2])/nlr[2][Ilr(i,j)];
     TV1=(tlr[0][Ilr(i,j)]*ff[1]+tlr[1][Ilr(i,j)]*ff[2])/nlr[2][Ilr(i,j)];
     rqp[0][0]=ff[0];
     rqp[1][0]=NV1;
     rqp[2][0]=ff[3];

// right values for riem

        dksi1=dksi2;
	dx = xc[IC(i+1,j)] - xc[IC(i,j)];
	dy = yc[IC(i+1,j)] - yc[IC(i,j)];
	dksi2 = sqrt(dx*dx+dy*dy);
	for( k=0; k < 6; k++)
	{
	delf1=(gdf[k][IC(i,j)]-gdf[k][IC(i-1,j)])/dksi1;
	delf2=(gdf[k][IC(i+1,j)]-gdf[k][IC(i,j)])/dksi2;
          if(i_appr)
           {
             delf=min_mod(delf1,delf2);
             ff[k]=gdf[k][IC(i,j)]-0.5*dksi1*delf;
           } 
          else ff[k]=gdf[k][IC(i,j)];
	}
Sr=ff[5];

     NV2=(nlr[0][Ilr(i,j)]*ff[1]+nlr[1][Ilr(i,j)]*ff[2])/nlr[2][Ilr(i,j)];
     TV2=(tlr[0][Ilr(i,j)]*ff[1]+tlr[1][Ilr(i,j)]*ff[2])/nlr[2][Ilr(i,j)];
     rqp[0][1]=ff[0];
     rqp[1][1]=NV2;
     rqp[2][1]=ff[3];

       for(l=0; l < 2; l++)
        for(k=0; k < 3; k++)
         par[k+3*l]=rqp[k][l]; 

       riemi(par,&alarm);
//       riemi_H2(par,&alarm);

        if(alarm )
          {
          printf("ALARM_RIEM FROM l_r : i=%d j=%d  \n",i,j);
          printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
          printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
          exit(-1);
          }

     r_o=par[6];  u_o=par[7] ;   p_o=par[8] ;
     s1 =par[9];  s2 =par[10];   s3 =par[11]; 

//     Sb=Riem_SC(Sl,Sr,u_o);

     v_o=TV1;
     if(s2 < 0.) v_o=TV2; 

     db=r_o; pb=p_o;
     ub=(u_o*nlr[0][Ilr(i,j)]+v_o*tlr[0][Ilr(i,j)])/nlr[2][Ilr(i,j)];
     vb=(u_o*nlr[1][Ilr(i,j)]+v_o*tlr[1][Ilr(i,j)])/nlr[2][Ilr(i,j)];
     e=pb/(GAM-1.)+0.5*db*(ub*ub+vb*vb);

// Eu_flux

//	w1 = wlr[0][Ilr(i,j)];
//	w2 = wlr[1][Ilr(i,j)];

	w1 = 0.;
	w2 = 0.;


     fllr[0][Ilr(i,j)]=db*(ub - w1)*nlr[0][Ilr(i,j)]+db*(vb-w2)*nlr[1][Ilr(i,j)];
     fllr[1][Ilr(i,j)]=(db*ub*(ub - w1)+pb)*nlr[0][Ilr(i,j)]+db*(ub - w1)*vb*nlr[1][Ilr(i,j)];
     fllr[2][Ilr(i,j)]=db*ub*(vb-w2)*nlr[0][Ilr(i,j)]+(db*vb*(vb-w2)+pb)*nlr[1][Ilr(i,j)];
     fllr[3][Ilr(i,j)]=(ub - w1)*(e+pb)*nlr[0][Ilr(i,j)]+(vb-w2)*(e+pb)*nlr[1][Ilr(i,j)];

// NS_flux
// !!!!!!!!!!!!!!!!!!!!

     if(Flag_turb == 1){
      muef = mu_ref * (pow((pb/db),0.76) + mu_turb_lr[Ilr(i,j)]);
     }
     else{
      muef=mu_ref*pow((pb/db),0.76);
     }
 

  t11=2.*(2.*dxflr[1][Ilr(i,j)]-dyflr[2][Ilr(i,j)])/3.;
  t12=(dyflr[1][Ilr(i,j)]+dxflr[2][Ilr(i,j)]);
  t22=2.*(2.*dyflr[2][Ilr(i,j)]-dxflr[1][Ilr(i,j)])/3.;

     kapef=muef/0.72;                                                                    
     kapef*=3.5;

     if(mu_max < muef) mu_max=muef;

     t11*=muef; 
     t12*=muef; 

//!!!!!!!!!!!!!!!!!!!!!
//if(i_Rey)
//  t12-=Sb*db;

     t21=t12;
     t22*=muef; 

n=Ilr(i,j);
//Full QGD
if(i_qgd > 1)
{
mum=mu_ref*pow((pb/db),0.76);
tau=mum/pb;
wzx=tau*(db*(ub*dxflr[1][n]+vb*dyflr[1][n])+dxflr[3][n]);
wzy=tau*(db*(ub*dxflr[2][n]+vb*dyflr[2][n])+dyflr[3][n]);
dv=dxflr[1][n]+dyflr[2][n];
rz=tau*(ub*dxflr[3][n]+vb*dyflr[3][n]+GAM*pb*dv);
t11+=ub*wzx+rz;
t12+=vb*wzy;
t21+=vb*wzx;
t22+=vb*wzy+rz;
}

if(coord > 0.)
{
t11-=2.*muef*vb/(3.*yc[IC(i,j)]);
t22-=2.*muef*vb/(3.*yc[IC(i,j)]);
}
     q1=kapef*dxflr[4][Ilr(i,j)];
     q2=kapef*dyflr[4][Ilr(i,j)];
     
     fllr[1][Ilr(i,j)]-=t11*nlr[0][Ilr(i,j)]+t12*nlr[1][Ilr(i,j)];
     fllr[2][Ilr(i,j)]-=t21*nlr[0][Ilr(i,j)]+t22*nlr[1][Ilr(i,j)];
     fllr[3][Ilr(i,j)]-=(ub*t11+vb*t12+q1)*nlr[0][Ilr(i,j)]+
                        (ub*t21+vb*t22+q2)*nlr[1][Ilr(i,j)];
// calc time-step
if(i > 2 && i < IT-2)
{
     dksi1=ndu[2][Idu(i-1,j)];
     dksi2=ndu[2][Idu(i,j)];
     dt1=1.e17;
     if( s1 < 0.)
     dt1=dksi1/(fabs(s1)+eforstep);
     dt2=1.e17;
     if( s3 > 0.)
     dt2=dksi2/(fabs(s3)+eforstep);
     dt1=min(dt1,dt2);
     dt2=0.25*dksi1*dksi2/muef;
     dt_max=min(dt1,dt2); 
     dt_lr=min(dt_max,dt_lr);
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EQUATION 5
//  fllr5[Ilr(i,j)]=(db*ub*Sb-dxflr5[0][Ilr(i,j)]-dxflr5[2][Ilr(i,j)]+muef*(vb*dxflr[1][n]+ub*dxflr[2][n]))*nlr[0][Ilr(i,j)]+
//                  (db*vb*Sb-dyflr5[1][Ilr(i,j)]-dyflr5[3][Ilr(i,j)]+muef*(vb*dyflr[1][n]+ub*dyflr[2][n]))*nlr[1][Ilr(i,j)];

   }

return dt_lr;
}


//*****************************************************

double min_mod(double a, double b)
{
 double tmp;
 tmp=0.;

 if( a*b > 0.)
  {
   tmp=a;
   if( fabs(a) > fabs(b) )  tmp=b;
  }
  return tmp;
}

/*
double min_mod(double a, double b)
// Van Albada
{
 double tmp;
 tmp=((b*b+eforstep)*a+(a*a+eforstep)*b)/(a*a+b*b+2.*eforstep);
  return tmp;
}
*/ 


//*****************************************************

double  flow_du(void)
{
 int i,j,k,l,lom,alarm;
 double dx,dy,db,ub,vb,pb,e,NV1,TV1,NV2,TV2;
 double t11,t12,t21,t22,q1,q2;
 double mum,mut,muef,kapm,kapt,kapef;
// static double dfx[6],dfy[6],ff[6];
 double dfx[6],dfy[6],ff[6];
 double xc3,yc3,xc4,yc4;
// static double par[32],rqp[3][2],r_o,u_o,v_o,p_o,om,s1,s2,s3;
 double par[32],rqp[3][2],r_o,u_o,v_o,p_o,om,s1,s2,s3;
 double deta1,deta2;
 double delf,delf1,delf2,s;
 double w1,w2;
 double dt_du,dt1,dt2,dt_max;
 double Sb,Sl,Sr;

int n;
double wzx,wzy,rz,dv; 
double pr1, pr2;

 lom=0;  om=0.;
 mu_max=0.;
 dt_du=1.e17;

 for(i=2; i < IT-2; i++)
  for(j=2; j <= JT-2; j++)
   {
// left(down)  values for riem
	dx=xc[IC(i,j-1)]-xc[IC(i,j-2)];
	dy=yc[IC(i,j-1)]-yc[IC(i,j-2)];
	deta1=sqrt(dx*dx+dy*dy);
	dx=xc[IC(i,j)]-xc[IC(i,j-1)];
	dy=yc[IC(i,j)]-yc[IC(i,j-1)];
	deta2=sqrt(dx*dx+dy*dy);

	for( k=0; k < 6; k++)
	{
	delf1=(gdf[k][IC(i,j-1)]-gdf[k][IC(i,j-2)])/deta1;
	delf2=(gdf[k][IC(i,j)]-gdf[k][IC(i,j-1)])/deta2;
	if(i_appr)
        {
	delf=min_mod(delf1,delf2);
	ff[k]=gdf[k][IC(i,j-1)]+0.5*deta2*delf;
        }
        else ff[k]=gdf[k][IC(i,j-1)]; 
	}
Sl=ff[5];

     NV1=(ndu[0][Idu(i,j)]*ff[1]+ndu[1][Idu(i,j)]*ff[2])/ndu[2][Idu(i,j)];
     TV1=(tdu[0][Idu(i,j)]*ff[1]+tdu[1][Idu(i,j)]*ff[2])/ndu[2][Idu(i,j)];
     rqp[0][0]=ff[0];
     rqp[1][0]=NV1;
     rqp[2][0]=ff[3];

// right (up) values for riem

	deta1=deta2;
	dx=xc[IC(i,j+1)]-xc[IC(i,j)];
	dy=yc[IC(i,j+1)]-yc[IC(i,j)];
	deta2=sqrt(dx*dx+dy*dy);
	for( k=0; k < 6; k++)
	{
	delf1=(gdf[k][IC(i,j)]-gdf[k][IC(i,j-1)])/deta1;
	delf2=(gdf[k][IC(i,j+1)]-gdf[k][IC(i,j)])/deta2;
        if(i_appr)
        {
	delf=min_mod(delf1,delf2);
        ff[k]=gdf[k][IC(i,j)]-0.5*deta1*delf;
        }
        else ff[k]=gdf[k][IC(i,j)];
	}
Sr=ff[5];

     NV2=(ndu[0][Idu(i,j)]*ff[1]+ndu[1][Idu(i,j)]*ff[2])/ndu[2][Idu(i,j)];
     TV2=(tdu[0][Idu(i,j)]*ff[1]+tdu[1][Idu(i,j)]*ff[2])/ndu[2][Idu(i,j)];
     rqp[0][1]=ff[0];
     rqp[1][1]=NV2;
     rqp[2][1]=ff[3];
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(j==2)
      {
      NV1=-NV2;
      TV1=TV2;
      rqp[0][0]=rqp[0][1];
      rqp[1][0]=-rqp[1][1];
      rqp[2][0]=rqp[2][1];
      }
/*
     if(j==JT-2)
      {
      NV2=-NV1;
      TV2=TV1;
      rqp[0][1]=rqp[0][0];
      rqp[1][1]=-rqp[1][0];
      rqp[2][1]=rqp[2][0];
      }
*/
       for(l=0; l < 2; l++)
        for(k=0; k < 3; k++)
         par[k+3*l]=rqp[k][l]; 

       riemi(par,&alarm);
//       riemi_H2(par,&alarm);

        if(alarm )
          {
          printf("ALARM_RIEM FROM l_r : i=%d j=%d  \n",i,j);
          printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
          printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
          exit(-1);
          }

     r_o=par[6];  u_o=par[7] ;   p_o=par[8] ;
     s1 =par[9];  s2 =par[10];   s3 =par[11]; 

//     Sb=Riem_SC(Sl,Sr,u_o);

     v_o=TV1;
     if(s2  < 0.) v_o=TV2; 
     db=r_o; pb=p_o;
     ub=(u_o*ndu[0][Idu(i,j)]+v_o*tdu[0][Idu(i,j)])/ndu[2][Idu(i,j)];
     vb=(u_o*ndu[1][Idu(i,j)]+v_o*tdu[1][Idu(i,j)])/ndu[2][Idu(i,j)];

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(j==2) ub=vb=0.;

     e=pb/(GAM-1.)+0.5*db*(ub*ub+vb*vb);
 

// Eu_flux
 w1=0.;
 w2=0.;

     fldu[0][Idu(i,j)]=db*(ub-w1)*ndu[0][Idu(i,j)]+db*(vb-w2)*ndu[1][Idu(i,j)];
     fldu[1][Idu(i,j)]=(db*ub*(ub-w1)+pb)*ndu[0][Idu(i,j)]+db*(ub-w1)*vb*ndu[1][Idu(i,j)];
     fldu[2][Idu(i,j)]=db*ub*(vb-w2)*ndu[0][Idu(i,j)]+(db*vb*(vb-w2)+pb)*ndu[1][Idu(i,j)];
     fldu[3][Idu(i,j)]=(ub-w1)*(e+pb)*ndu[0][Idu(i,j)]+(vb-w2)*(e+pb)*ndu[1][Idu(i,j)];



// NS_flux

     if(Flag_turb == 1){
      muef = mu_ref * (pow((pb/db),0.76) + mu_turb_du[Idu(i,j)]);
     }
     else{
      muef=mu_ref*pow((pb/db),0.76);
     }

  t11=2.*(2.*dxfdu[1][Idu(i,j)]-dyfdu[2][Idu(i,j)])/3.;
  t12=(dyfdu[1][Idu(i,j)]+dxfdu[2][Idu(i,j)]);   
  t22=2.*(2.*dyfdu[2][Idu(i,j)]-dxfdu[1][Idu(i,j)])/3.;

     kapef=muef/0.72; 
     kapef*=3.5;
     if(mu_max < muef) mu_max=muef;

     t11*=muef; 
     t12*=muef; 
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//if(i_Rey)
//     t12-=db*Sb;

     t21=t12;
     t22*=muef; 

n=Idu(i,j);
if(i_qgd > 1)
{
mum=mu_ref*pow((pb/db),0.76);
tau=mum/pb;
wzx=tau*(db*(ub*dxfdu[1][n]+vb*dyfdu[1][n])+dxfdu[3][n]);
wzy=tau*(db*(ub*dxfdu[2][n]+vb*dyfdu[2][n])+dyfdu[3][n]);
dv=dxfdu[1][n]+dyfdu[2][n];
rz=tau*(ub*dxfdu[3][n]+vb*dyfdu[3][n]+GAM*pb*dv);
t11+=ub*wzx+rz;
t12+=vb*wzy;
t21+=vb*wzx;
t22+=vb*wzy+rz;
}

if(coord > 0.)
{
t11-=2.*muef*vb/(3.*yc[IC(i,j)]);
t22-=2.*muef*vb/(3.*yc[IC(i,j)]);
}
     q1=kapef*dxfdu[4][Idu(i,j)];
     q2=kapef*dyfdu[4][Idu(i,j)];
     
     fldu[1][Idu(i,j)]-=t11*ndu[0][Idu(i,j)]+t12*ndu[1][Idu(i,j)];
     fldu[2][Idu(i,j)]-=t21*ndu[0][Idu(i,j)]+t22*ndu[1][Idu(i,j)];
     fldu[3][Idu(i,j)]-=(ub*t11+vb*t12+q1)*ndu[0][Idu(i,j)]+
                        (ub*t21+vb*t22+q2)*ndu[1][Idu(i,j)];

if( j > 2 && j < JT-2)
{
     deta1=nlr[2][Ilr(i,j-1)];
     deta2=nlr[2][Ilr(i,j)];
     dt1=1.e17;
     if(s1 < 0.)
     dt1=deta1/(fabs(s1)+eforstep);
     dt2=1.e17;
     if(s3 > 0.)
     dt2=deta2/(fabs(s3)+eforstep);
     dt1=min(dt1,dt2);
     dt2=0.25*deta1*deta2/muef;
     dt_max=min(dt1,dt2);
     dt_du=min(dt_max,dt_du);
}


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EQUATION5

//  fldu5[Idu(i,j)]=(db*ub*Sb-dxfdu5[0][Idu(i,j)]-dxfdu5[2][Idu(i,j)]+muef*(vb*dxfdu[1][n]+ub*dxfdu[2][n]))*ndu[0][Idu(i,j)]+
//                  (db*vb*Sb-dyfdu5[1][Idu(i,j)]-dyfdu5[3][Idu(i,j)]+muef*(vb*dyfdu[1][n]+ub*dyfdu[2][n]))*ndu[1][Idu(i,j)];

                                                  
   }
return dt_du;
}


//*****************************************************

void new_gdf(void)
{
 int i,j,n,k;
 double p,e;
 double rhs[5],a0[5],a1[5],a_[5],fl[5],dxfc[5],dyfc[5];
 double q1,q2,mum,c,t11,t12,t21,t22,t33,muef,kapef;
 double fact;
 double fl5,rhs5,a05,a_5,al5;

 double gdftmp1,gdftmp2,gdftmp3,tmp1,tmp2,tmp3;

 double ddt1,ddt2,d2dt2;

 rhs[0]=rhs[1]=rhs[2]=rhs[3]=rhs[4]=0.;
 RES=0.;

 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    n=IC(i,j);
    for(k=0; k < 5; k++)
     {
      dxfc[k]=0.25*(dxflr[k][Ilr(i,j)]+dxflr[k][Ilr(i+1,j)]+
                    dxfdu[k][Idu(i,j)]+dxfdu[k][Idu(i,j+1)]);
      dyfc[k]=0.25*(dyflr[k][Ilr(i,j)]+dyflr[k][Ilr(i+1,j)]+
                    dyfdu[k][Idu(i,j)]+dyfdu[k][Idu(i,j+1)]);
     }

if(coord > 0.)
{  
/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  to check                                  
   if(Flag_turb == 1)
    {
    muef=mum=mu_ref*(pow((gdf[3][n]/gdf[0][n]),0.76) + 0.25 * (mu_turb_lr[0][IC(i,j)]+mu_turb_lr[1][IC(i,j)]+ 
                                                               mu_turb_du[0][IC(i,j)]+mu_turb_du[1][IC(i,j)]));
    }
    else{muef=mum=mu_ref*pow((gdf[3][n]/gdf[0][n]),0.76);}
*/
    muef=mum=mu_ref*pow((gdf[3][n]/gdf[0][n]),0.76);
    kapef=3.5*muef/0.72;
    t11=2.*muef*(2.*dxfc[1]-dyfc[2]-gdf[2][n]/yc[n])/3.;
    t12=muef*(dyfc[1]+dxfc[2]);
    t21=t12;
    t22=2.*muef*(2.*dyfc[2]-dxfc[1]-gdf[2][n]/yc[n])/3.;
    t33=2.*muef*(-dxfc[1]-dyfc[2]+2.*gdf[2][n]/yc[n])/3.;
    q1=kapef*dxfc[4];
    q2=kapef*dyfc[4];
    rhs[0]=-gdf[0][n]*gdf[2][n]/yc[n];
    rhs[1]=-(gdf[0][n]*gdf[1][n]*gdf[2][n]-t12)/yc[n];
    rhs[2]=-(gdf[0][n]*gdf[2][n]*gdf[2][n]-t22+t33)/yc[n];
    e=gdf[3][n]/(GAM-1.)+0.5*gdf[0][n]*
      (gdf[1][n]*gdf[1][n]+gdf[2][n]*gdf[2][n]);
    rhs[3]=-(gdf[2][n]*(e+gdf[3][n])+q2-gdf[1][n]*t12-gdf[2][n]*t22)/yc[n];
}

muef=mum=mu_ref*pow((gdf[3][n]/gdf[0][n]),0.76);

if(i_qgd) tau=muef/gdf[3][n];
else tau=0.;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!????????????????????????????????
//tau*=2.;
          
    a0[0]=gdf[0][n];
    a0[1]=gdf[0][n]*gdf[1][n];
    a0[2]=gdf[0][n]*gdf[2][n];
    e=gdf[3][n]/(GAM-1.)+0.5*gdf[0][n]*
      (gdf[1][n]*gdf[1][n]+gdf[2][n]*gdf[2][n]);
    a0[3]=e;

    a_[0]=gdf_[0][n];
    a_[1]=gdf_[0][n]*gdf_[1][n];
    a_[2]=gdf_[0][n]*gdf_[2][n];
    e=gdf_[3][n]/(GAM-1.)+0.5*gdf_[0][n]*
      (gdf_[1][n]*gdf_[1][n]+gdf_[2][n]*gdf_[2][n]);
    a_[3]=e;

     for(k=0; k < 4; k++)
     {
     fl[k]=(fllr[k][Ilr(i+1,j)]-fllr[k][Ilr(i,j)]+fldu[k][Idu(i,j+1)]-fldu[k][Idu(i,j)]);

/*       
if(i_qgd)
  {
// 3 level t_app2
//     a1[k]=a_[k]*(1.-2.*tau/dt)+4.*tau*a0[k]/dt-2.*dt*fl[k]/vol[n]+2.*dt*rhs[k];
//     a1[k]=a1[k]/(1.+2.*tau/dt);

// 2 level t_app1
     a1[k]=a_[k]*(-tau/dt)+a0[k]*(1.+2*tau/dt)-dt*fl[k]/vol[n]+dt*rhs[k];
     a1[k]=a1[k]/(1.+tau/dt);
  }
else
*/

     a1[k]=a0[k]-dt*fl[k]/vol[n]+dt*rhs[k];
  
     }
    RES=max(RES,fabs(fl[0]));

//     ddt1=(gdf[0][n]*gdf[1][n]*gdf[2][n]-gdf_[0][n]*gdf_[1][n]*gdf_[2][n])/dt;

//     a05=gdf[5][n]*gdf[0][n];  a_5=gdf_[5][n]*gdf_[0][n]; 
     

for(k=0; k < 4; k++)
     gdf_[k][n]=gdf[k][n];
     
    gdf[0][n]=a1[0];
    gdf[1][n]=a1[1]/a1[0];
    gdf[2][n]=a1[2]/a1[0];
    p=(GAM-1.)*(a1[3]-0.5*(a1[1]*gdf[1][n]+a1[2]*gdf[2][n]));
    gdf[3][n]=p;
    gdf[4][n]=p/gdf[0][n];

/*
//EQUATION UV

c=sqrt(GAM*gdf[3][n]/gdf[0][n]);
//tau=muef/gdf[3][n]+sqrt(lxmin*lymin)/c;
tau=muef/gdf[3][n];
//tau=0.;

     fl5=(fllr5[Ilr(i+1,j)]-fllr5[Ilr(i,j)]+fldu5[Idu(i,j+1)]-fldu5[Idu(i,j)]);

     rhs5=0.;
     ddt2=(gdf[0][n]*gdf[1][n]*gdf[2][n]-gdf_[0][n]*gdf_[1][n]*gdf_[2][n])/dt;
     d2dt2=(ddt2-ddt1)/dt;
     rhs5-=0.5*tau*d2dt2;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!16.08.21
     rhs5+=rhs5_1[n];

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2) 01.09.21
//!!!!!!!!!!!!!!!!!?????????????????????????????????120122
     rhs5-=gdf_[0][n]*gdf[5][n]*(dxfc[1]+dyfc[2]);
//     rhs5+=gdf_[0][n]*gdf[5][n]*(dxfc[1]+dyfc[2]);


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!04.09.21
rhs5+=2.*muef*(dxfc[1]*dxfc[2]+dyfc[1]*dyfc[2]);
     
//     al5=a_5*(-tau/dt)+a05*(1.+2*tau/dt)-dt*fl5/vol[n]+dt*rhs5;
//     al5=al5/(1.+tau/dt);

     al5=a05-dt*fl5/vol[n]+dt*rhs5;
     gdf_[5][n]=gdf[5][n];
     gdf[5][n]=al5/a1[0];
*/
     
   }

}

//*****************************************************

void turb_lr(void) // Функция, заполняющая массив mu_turb_lr вязкостями;
{  
 int i, j;
 int Flag = 0;
int JTH;
 double eps = 1.e-20;
 double mu_To, mu_Ti;
 double k = 0.4, k_caps = 0.0168, C_cp = 26., C_kleb = 0.3, C_wk = 0.25, A = 26.;
 double F_wake, Gamma, D;
 double F_max = 0., Y_max = 0, F_y, t_w, nu_;
 double U_max = -9.9e+9;
 double U_dif;
 double nu, omega, mu_m;
 double mu_To_max = 0., mu_Ti_max = 0., bug_min = 9.9e+9;
 double gdf_0, gdf_1, gdf_3, y_,YUP;

//JTH=(JT-3-2)/2;


 for(i=2; i < IT-2; i++){

//  YUP=y[IN(i,JT-2)];

//DOWN => UP

  F_max = 0.;
  Y_max = 0.;
  U_max = 0.;
   for(j=2; j < JT-2; j++){

    y_ = (yc[IC(i,j)] + yc[IC(i-1,j)]) / 2.;
    gdf_0 = (gdf[0][IC(i,j)] + gdf[0][IC(i-1,j)]) / 2.;
    gdf_1 = (gdf[1][IC(i,j)] + gdf[1][IC(i-1,j)]) / 2.;
    gdf_3 = (gdf[3][IC(i,j)] + gdf[3][IC(i-1,j)]) / 2.;

    mu_m=mu_ref*pow((gdf_3/gdf_0),0.76);
    nu = mu_m / gdf_0;// динамическая вязкость
if(j==2)
{
    t_w = nu * dyflr[1][Ilr(i,j)];// напряжение трения на стенке
    nu_ = sqrt(fabs(t_w / gdf_0));// динамическая скорость
}
    omega = dxflr[2][Ilr(i,j)] - dyflr[1][Ilr(i,j)];
    F_y = y_ * fabs(omega) * (1 - exp((-y_ * nu_) / (A * nu))); // F(y)
    // Определяем Y_max, F_max
    if(F_max < F_y){
     F_max = F_y;
     Y_max = y_;
    }
    // Определим U_dif = | U_max(x) - U_min(x) |
    U_max = max(U_max, gdf_1);
   }
  U_dif = fabs(U_max); //  - U_min
  F_wake = min(Y_max * F_max, C_wk * Y_max * pow(U_dif, 2.) / F_max );

  Flag=0;
  for(j=2; j < JT-2; j++){
   y_ = (yc[IC(i,j)] + yc[IC(i-1,j)]) / 2.;
   gdf_0 = (gdf[0][IC(i,j)] + gdf[0][IC(i-1,j)]) / 2.;
   gdf_1 = (gdf[1][IC(i,j)] + gdf[1][IC(i-1,j)]) / 2.;
   gdf_3 = (gdf[3][IC(i,j)] + gdf[3][IC(i-1,j)]) / 2.;

   Gamma = pow((1 + 5.5 * pow((C_kleb * y_ / Y_max), 6.)), -1. );
   mu_m=mu_ref*pow((gdf_3/gdf_0),0.76);
   nu = mu_m / gdf_0;// динамическая вязкость
if(j==2)
{
   t_w = nu * dyflr[1][Ilr(i,j)];// напряжение трения на стенке
   nu_ = sqrt(fabs(t_w / gdf_0));// динамическая скорость
}
   omega = dxflr[2][Ilr(i,j)] - dyflr[1][Ilr(i,j)];
   D = pow((1 - exp((-y_ * nu_) / (A * nu))), 2.);

   
  mu_Ti = pow((k * y_), 2) * D * fabs(omega); // Turb in
   mu_To = k_caps * C_cp * F_wake * Gamma; // Turb out
//printf("lr:  %d %d %e %e \n",i,j,mu_Ti,mu_To);

   //Отбор вязкости
   if(mu_To  > mu_Ti && Flag == 0)
    {
      mu_turb_lr[Ilr(i,j)] =  mu_Ti;
    }
   else
    {
      mu_turb_lr[Ilr(i,j)] = mu_To;
      Flag=1;
    }
   mu_turb_lr[Ilr(i,j)] = mu_turb_lr[Ilr(i,j)] * gdf_0;
// printf("lr:  %d %d turb=%e \n",i,j,mu_turb_lr[Ilr(i,j)]);
  }

    }

}

//*****************************************************


void turb_du(void) // Функция, заполняющая массив mu_turb_du вязкостями;
{  
 int i, j;
 int Flag = 0;
int JTH;

 double eps = 1.e-20;
 double mu_To, mu_Ti;
 double k = 0.4, k_caps = 0.0168, C_cp = 26., C_kleb = 0.3, C_wk = 0.25, A = 26.;
 double F_wake, Gamma, D;
 double F_max = 0, Y_max = 0, F_y, t_w, nu_;
 double U_max = -9.9e+9;
 double U_dif;
 double nu, omega, mu_m;
 double mu_To_max = 0., mu_Ti_max = 0., bug_min = 9.9e+9;
 double gdf_0, gdf_1, gdf_3, y_,YUP;

//JTH=(JT-3-2+1)/2;

 for(i=2; i < IT-2; i++){

//YUP=y[IN(i,JT-2)];

//DOWN => UP
     F_max = 0.;
     Y_max = 0.;
     U_max = 0.;
   for(j=2; j <=JT-2; j++){

    y_ = y[IN(i,j)];
    gdf_0 = (gdf[0][IC(i,j)] + gdf[0][IC(i,j-1)]) / 2.;
    gdf_1 = (gdf[1][IC(i,j)] + gdf[1][IC(i,j-1)]) / 2.;
    gdf_3 = (gdf[3][IC(i,j)] + gdf[3][IC(i,j-1)]) / 2.;

    mu_m=mu_ref*pow((gdf_3/gdf_0),0.76);
    nu = mu_m / gdf_0;// динамическая вязкость
if(j==2)
{
    t_w = nu * dyfdu[1][Idu(i,j)];// напряжение трения на стенке
    nu_ = sqrt(fabs(t_w / gdf_0));// динамическая скорость
}
    omega = dxfdu[2][Idu(i,j)] - dyfdu[1][Idu(i,j)];
    F_y = y_ * fabs(omega) * (1. - exp((-y_ * nu_) / (A * nu))); // F(y)
    // Определяем Y_max, F_max
// 
//printf("du:  %d %d  Y=%e  F=%e  Om=%e \n",i,j,y_,F_y,omega);

    if(F_max < F_y){
     F_max = F_y;
     Y_max = y_;
    }
    // Определим U_dif = | U_max(x) - U_min(x) |
    U_max = max(U_max, gdf_1);
   }

  U_dif = fabs(U_max); // 
  if(Y_max < y[IN(i,3)])
  {
  Y_max = y[IN(i,3)];
  F_wake=0.;
  }
  else
  F_wake = min(Y_max * F_max, C_wk * Y_max * pow(U_dif, 2.) / F_max );

  Flag=0;

  for(j=2; j <= JT-2; j++){
   y_ = y[IN(i,j)];
   gdf_0 = (gdf[0][IC(i,j)] + gdf[0][IC(i,j-1)]) / 2.;
   gdf_1 = (gdf[1][IC(i,j)] + gdf[1][IC(i,j-1)]) / 2.;
   gdf_3 = (gdf[3][IC(i,j)] + gdf[3][IC(i,j-1)]) / 2.;

   Gamma = pow((1 + 5.5 * pow((C_kleb * y_ / Y_max), 6.)), -1. );
   mu_m=mu_ref*pow((gdf_3/gdf_0),0.76);
   nu = mu_m / gdf_0;// динамическая вязкость
if(j==2)
{
   t_w = nu * dyfdu[1][Idu(i,j)];// напряжение трения на стенке
   nu_ = sqrt(fabs(t_w / gdf_0));// динамическая скорость
}
   omega = dxfdu[2][Idu(i,j)] - dyfdu[1][Idu(i,j)];
   D = pow((1 - exp((-y_ * nu_) / (A * nu))), 2);
   
   mu_Ti = pow((k * y_), 2) * D * fabs(omega); // Turb in
   mu_To = k_caps * C_cp * F_wake * Gamma; // Turb out
//printf("du:  %d %d %e %e %e %e %e  \n",i,j,y_,D,omega,F_wake,Gamma);

   
   //Отбор вязкости
   if(mu_Ti < mu_To && Flag == 0)
    {
      mu_turb_du[Idu(i,j)] =  mu_Ti;
    }
   else
    {
      mu_turb_du[Idu(i,j)] = mu_To;
      Flag=1;
    }                       
// printf("du:  %d %d turb_i=%e turb_o=%e \n",i,j,mu_Ti,mu_To);

   mu_turb_du[Idu(i,j)] = mu_turb_du[Idu(i,j)] * gdf_0;
// printf("du:  %d %d turb=%e \n",i,j,mu_turb_du[Idu(i,j)]);
// printf("du:  %d %d turb_i=%e turb_o=%e \n",i,j,mu_Ti,mu_To);

  }

  }
}
