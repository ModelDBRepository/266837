#include <stdio.h>
#include <math.h>  /*per sin i cos*/
#include <stdlib.h> /*pel malloc i calloc*/
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

int neur=50;
double imax;
double imaxa;
double dgap;
int ngap;
double eps;

int func (double t, const double x[], double dx[], void *params);

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);

double inaglut1(double naea, double glut);

double max(a,b){
   double max;
   max=a;
   if (b>a){
    max=b;
   }
   return max;
}

double min(a,b){
   double min;
   min=a;
   if (b<a){
    min=b;
   }
   return min;
}


int mod (int a, int b){
   int ret = a % b;
   if(ret < 0){
     ret+=b;
   }
   return ret;
}

int main(int argc, char * argv[]){
    double *x,*dx;
    int ndim=11*neur;;
    int i;
    FILE * output;
    double tol=1.0e-8;
    double *xa;
    double mu;
    
    //const gsl_odeiv_step_type * T  = gsl_odeiv_step_rk8pd;  // Runge-Kutta Prince-Dormand 8-9.
    //const gsl_odeiv_step_type * T  = gsl_odeiv_step_rk2imp; // Runge-Kutta 2 implicit.
    const gsl_odeiv_step_type * T  = gsl_odeiv_step_rk4imp; // Runge-Kutta 4 implicit.
    //const gsl_odeiv_step_type * T  = gsl_odeiv_step_bsimp; // Burlirsch-Stoer implicit.. Need jacobian!
    //const gsl_odeiv_step_type * T  = gsl_odeiv_step_gear1; // Gear 1 implicit.
    //const gsl_odeiv_step_type * T  = gsl_odeiv_step_gear2; // Gear 2 implicit.

//    const gsl_odeiv_step_type * T=gsl_odeiv_step_rkf45;

    gsl_odeiv_step * s  = gsl_odeiv_step_alloc (T, ndim);
    gsl_odeiv_control * c  = gsl_odeiv_control_y_new (tol, 0.0);
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (ndim);

    gsl_odeiv_system sys = {func, jac, ndim, &mu};
    
    double t = 0.0; 
    double tf=160000;//150000.;
    double h = 0.05;
    int status;
    int nite;
    int ne=0;
    
    char name1[50];
    int number1;

    if(argc != 6){
     printf("usage: %s dgap ngap eps imax imaxa\n", argv[0]);
     abort();
    }

    sscanf(argv[1], "%lg", &dgap);
    sscanf(argv[2], "%d", &ngap);
    sscanf(argv[3], "%lg", &eps);
    sscanf(argv[4], "%lg", &imax);
    sscanf(argv[5], "%lg", &imaxa);
    
    sprintf(name1, "a2_int_epsoff_dgap%lg_ngap%d_eps%lg_imax%lg_imaxa%lg.dat", dgap,ngap,eps,imax,imaxa);  

   /*open the data file*/
    output=fopen(name1, "w");
    if (output==NULL){
    printf("Error in opening data file\n");
    exit(1);
    }

    /*memory allocation*/
    x=(double *)calloc(ndim, sizeof(double));
    xa=(double *)calloc(ndim, sizeof(double));
    dx=(double *)calloc(ndim, sizeof(double));

   /*Initial point*/
    for (i=0;i<neur;i++){
      x[i]=-70.; //v
      xa[i]=x[i];
      x[i+neur]=.9751; //hp
      x[i+neur*2]=.25512; //n
      x[i+neur*3]=3.5;//ke
   /*   if (i<5){
        x[i+neur*3]=20.;//90.; //ke 
      }
      if (i>4){
        x[i+neur*3]=3.5; //ke
      }*/
      x[i+neur*4]=138.116; //ki
      x[i+neur*5]=133.574; //nae
      x[i+neur*6]=4.297; //nai  
      x[i+neur*7]=-75.1757; //va
      x[i+neur*8]=124.; //kia
      x[i+neur*9]=6.6; //naia
      x[i+neur*10]=0.;//sna
    }


   fprintf(output,"%.5f",t);
    for (i=0;i<neur;i++){
     fprintf(output," %.5f ", x[i]);
     } 
    fprintf(output,"\n");
    
    
   /*integrate*/

    t=0.;
    nite=0.;

    while (t<tf){

//turn epsilon off
  
   /*   if (t>10000){
       eps=0;
     }*/

     status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tf, &h, x);
   
    if (status != GSL_SUCCESS)
     break;

   if (((x[23]+30)*(xa[23]+30)<0) && x[23]>-30){
     ne++;
   }
   if (((x[24]+30)*(xa[24]+30)<0) && x[24]>-30){
     ne++;
   }
   if (((x[25]+30)*(xa[25]+30)<0) && x[25]>-30){
     ne++;
   }
   
   if (((x[26]+30)*(xa[26]+30)<0) && x[26]>-30){
     ne++;
   }

   if (ne==4){
     eps=0;
   }

   // if (nite==1){
    fprintf(output,"%.5f",t);
    for (i=0;i<neur;i++){
//     fprintf(output," %.5f %.5f", x[i], x[i+neur*3]);
     xa[i]=x[i];  
     fprintf(output," %.5f ", x[i]);
     } 
//    fprintf(output,"%.5f %d\n", eps, ne);
    fprintf(output,"\n");
    printf("%.2f %.4f\n", t, eps);
     }

    free(x);
    free(dx);
    free(xa);

    fclose(output);

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);
  
}

/* Vector field we want to integrate*/


// Fast sodium
double ina(double v, double n, double vna){ 
   double gna=3.;//50
   double thm=-34.;
   double sigm=5.;
   double phih=.05;
   double minf;
   double ina;

   minf=1./(1.+exp(-(v-thm)/sigm));
   ina=gna*pow(minf,3)*(1.-n)*(v-vna);

   return ina;
  
}

//NaP
double inap(double v, double hp, double vna){
   double gnap=.4;//0.8
   double thmp=-40.,sigmp=6.;
   double minfp, inap;

   minfp=1./(1.+exp(-(v-thmp)/sigmp));
   inap=gnap*minfp*hp*(v-vna);

   return inap;
}

double hinfp(double v){
   double thhp=-48.,sighp=-6.;
   double hinfp;
   hinfp=1./(1.+exp(-(v-thhp)/sighp));
   return hinfp;
}

double tauhp(double v){
   double tauhp;
   double vt=-49.,sig=6.;
   double taubar=10000.;
   tauhp=taubar/cosh((v-vt)/(2.*sig));
   return tauhp;
}

// IK

double ik(double v, double n, double vk){
   double gk=5.;
   double ik;
   ik=gk*(pow(n,4))*(v-vk);

   return ik;
}

double taun(double v){
   double taun0=.05,taun1=.27,thnt=-40.,sn=-12.;
   double taun;
   taun=taun0+taun1/(1.+exp(-(v-thnt)/sn));
   return taun;
}


double ninf(double v){
   double thn=-55.,sgn=14.;
   double ninf;
   ninf=1./(1.+exp(-(v-thn)/sgn));
   return ninf;
}

// Ileak
double il(double v){
   double gl=.3;
   double vl=-70.;
   double il;
   il=gl*(v-vl);
   return il;
}




double ika(double v, double frt, double ke, double ki){
  double gamma=0.2;
  double pka=6.e-6;
  double F=96485;
  double ika;
  double phia;
  phia=(1./frt)*v;
  if (fabs(phia)>1.0e-10){
  ika=(1.-gamma)*pka*F*phia*(ke*exp(-phia)-ki)/(exp(-phia)-1.);
  }
  else{
  ika=(1.-gamma)*pka*F*(ki-ke);
  }
  return ika;
}

double inaa(double v, double frt, double nae,double nai){
  double pnaa=0.015e-6;
  double F=96485;
  double inaa; 
  double phia;
  phia=(1./frt)*v;
  if (fabs(phia)>1.0e-10){
  inaa=pnaa*F*phia*(nae*exp(-phia)-nai)/(exp(-phia)-1.);
  }
  else{
  inaa=pnaa*F*(nai-nae);
  }
  return inaa;
}


// NMDA
double inmda(double v, double sn){
   double thetat=-10.,trise=2.;
   double gnmda=0;
   double vnmda=0;
   double binf;
   double inmda;
   
   binf=1./(1.+exp(-(v-thetat)/16.13));
   inmda=gnmda*sn*binf*(v-vnmda);
   return inmda;
}
   

/*double sinfa(double v){
  double theta=0.;
  double kappa=1.;
  double sinfa;
  sinfa=1./(1.+exp(-((v-theta)/kappa)));
  return sinfa;
}  

double inaglut(double naea, double glut){

  double tg=.01, kg=.0001;
  double thc=40., kc=2.;
  double sinfn, sinfg;
  double inaglut;
  double gnaa=8.; //=0.; //8.;

  sinfg=1./(1.+exp(-((glut-tg)/kg)));
  sinfn=1./(1.+exp(-(naea-thc)/kc));
  inaglut=gnaa*sinfn*sinfg;
 
  return inaglut;
}*/

double gap(double x, double y, double a, double b, double frt){
  double phigap;
  double F=96485.;
  double gap;
  phigap=(x-y)/frt;
  if (fabs(phigap)>1.0e-10){
  gap=F*phigap*((b*exp(-phigap)-a)/(exp(-phigap)-1.));
  }
  else{
    gap=F*(a-b);
  }
 // printf("%.10f %.10f %.10f %.10f %.10f %.10f\n",x,y,a,b,phigap, gap);
 // exit(1);
  return gap;
}




int func (double t, const double x[], double dx[], void *params)
{
   //double mu = *(double *)params;
   
   double R=8310;
   double Temp=310.;
   double F=96485;
   double frt=R*Temp/F;
   double su=922;
   double voln=2160;
   double delta=.1;
   double sa=1600;
   double vola=2000;
   double vole=delta*(voln+vola);
   double c1=10.*su/(F*voln);
   double c2=10.*su/(F*vole);
   double c3a=10*sa/(F*vola);
   double c3e=10*sa/(F*vole);
   double iapp=0.;
   double phin=0.8; 
   double phih=0.05;
//   double imax=5;
   double glut=0.;

  
   double tdecay=1.,alphan=.5; 

   int i,j,k;

   double v[neur],hp[neur],n[neur],ke[neur],ki[neur],nae[neur],nai[neur],kia[neur],naia[neur],va[neur],sna[neur];
   double vna[neur], vk[neur];
   double ipump[neur], inapump[neur], ikpump[neur];
   double kdiff[neur],nadiff[neur];


 
    for (i=0;i<neur;i++){
      v[i]=x[i]; //v
      hp[i]=x[i+neur]; //hp
      n[i]=x[i+neur*2]; //n
      ke[i]=x[i+neur*3]; //ke
      ki[i]=x[i+neur*4]; //ki
      nae[i]=x[i+neur*5]; //nae
      nai[i]=x[i+neur*6]; //nai
      va[i]=x[i+neur*7]; //va
      kia[i]=x[i+neur*8]; //kia
      naia[i]=x[i+neur*9]; //naia
      sna[i]=x[i+neur*10]; //sna
      //glut[i]=x[i+neur*13]; //glut
    }

    double dk=0.002;
    double dna=0.00133;

    for (i=0;i<neur;i++){
      ipump[i]=imax/(pow(1.+2./ke[i],2)*pow(1.+7.7/nai[i],3));
      inapump[i]=3.*ipump[i];
      ikpump[i]=-2*ipump[i];
      vk[i]=frt*log(ke[i]/ki[i]);
      vna[i]=frt*log(nae[i]/nai[i]);
     }

   double ke0=3.5;
   double ken1=3.5;
   double nae0=135;
   double naen1=135;

//kdiff[0]=dk*(ke[1]-ke[0]);
kdiff[0]=dk*(ke0+ke[1]-2.*ke[0]);
for (i=1;i<(neur-1);i++){
 kdiff[i]=dk*(ke[i+1]+ke[i-1]-2.*ke[i]);
}
kdiff[neur-1]=dk*(ken1+ke[neur-2]-2.*ke[neur-1]);
//kdiff[neur-1]=dk*(ke[neur-2]-ke[neur-1]);


//nadiff[0]=dna*(nae[1]-nae[0]);
nadiff[0]=dna*(nae0+nae[1]-2.*nae[0]);
for (i=1;i<(neur-1);i++){
  nadiff[i]=dna*(nae[i+1]+nae[i-1]-2.*nae[i]);
}
nadiff[neur-1]=dna*(naen1+nae[neur-2]-2.*nae[neur-1]);
//nadiff[neur-1]=dna*(nae[neur-2]-nae[neur-1]);

   
//Astrocyte

//double dgap=0.;
double gkir=50;
//double imaxa=5.;
double cma=1.;
double iappa=0.;
double pkgap,pnagap;
double pka=6.0e-6;
double pnaa=0.015e-6;

//double buff=1; //vs 0 propagation
//double vbuff=-70;

double ipumpa[neur], igap[neur],inagap[neur],ikgap[neur];
//int ngap=5;

pkgap=dgap*pka;
pnagap=0.8*pkgap;


for (i=0; i<neur;i++){
 ikgap[i]=0;
 inagap[i]=0;
 for (j=max(i-ngap,0);j<min(i+ngap+1,neur);j++){
   if (j!=i){
//   printf("(i,j) %d %d\n",i,j);
   k=j;
   //printf("%d %d %d\n",i,j,k);
   ikgap[i]+=pkgap*gap(va[i],va[k],kia[i],kia[k],frt);
   inagap[i]+=pnagap*gap(va[i],va[k],naia[i],naia[k],frt);
   }
 }
}

//exit(1);

for (i=0; i<neur;i++){
 ipumpa[i]=imaxa/(pow(1.+2./ke[i],2)*pow(1.+10./naia[i],3));
 igap[i]=inagap[i]+ikgap[i];
 //va[i]=buff*vbuff+(1-buff)*vb[i];
}


double dkx=0.;//.0002;
double dnax=0.;//.0002;
double kx=3.5;
double nax=140.;

/*printf("%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",dx[i],dx[i+neur],
dx[i+2*neur],dx[i+3*neur],dx[i+4*neur],dx[i+5*neur],dx[i+6*neur],dx[i+7*neur],dx[i+8*neur],dx[i+9*neur],dx[i+10*neur]);*/

//All the equations
for (i=0; i<neur;i++){
  dx[i]=-(ina(v[i],n[i],vna[i])+inap(v[i],hp[i],vna[i])+ik(v[i],n[i],vk[i])+ipump[i]+il(v[i])+inmda(v[i],sna[i]))+iapp; //v
  dx[i+neur]=phih*(hinfp(v[i])-hp[i])/tauhp(v[i]); //hp
  dx[i+neur*2]= phin*(ninf(v[i])-n[i])/taun(v[i]); //n
//  if (i<12.5 || i>15.5){
  if ((neur/2.-2.5)<i && i<(neur/2.+1.5)){
//    printf("%d %.3f\n ",i,eps);
    dx[i+neur*3]=c2*(ik(v[i],n[i],vk[i])+ikpump[i])+kdiff[i]+c3e*(ika(va[i],frt,ke[i],kia[i])-2.*ipumpa[i])+eps+dkx*(kx-ke[i]); //ke  
      }
  else{
//  if (i>12.5 && i<15.5){
    dx[i+neur*3]=c2*(ik(v[i],n[i],vk[i])+ikpump[i])+kdiff[i]+c3e*(ika(va[i],frt,ke[i],kia[i])-2.*ipumpa[i])+dkx*(kx-ke[i]); //ke
  }
  //dx[i+neur*3]=c2*(ik(v[i],n[i],vk[i])+ikpump[i])+kdiff[i]+c3e*(ika(va[i],frt,ke[i],kia[i])-2.*ipumpa[i]); //ke
  dx[i+neur*4]=-c1*(ik(v[i],n[i],vk[i])+ikpump[i]); //ki
  dx[i+neur*5]=c2*(ina(v[i],n[i],vna[i])+inap(v[i],hp[i],vna[i])+inapump[i])+nadiff[i]+c3e*(inaa(va[i],frt,nae[i],naia[i])+3.*ipumpa[i])+dnax*(nax-nae[i]); //nae
  dx[i+neur*6]=-c1*(ina(v[i],n[i],vna[i])+inap(v[i],hp[i],vna[i])+inapump[i]); //nai
  dx[i+neur*7]=-(ika(va[i],frt,ke[i],kia[i])+inaa(va[i],frt,nae[i],naia[i])+ipumpa[i]+igap[i]-iappa)/cma; //va
  dx[i+neur*8]=-c3a*(ika(va[i],frt,ke[i],kia[i])-2*ipumpa[i]+ikgap[i]); //kia
  dx[i+neur*9]=-c3a*(inaa(va[i],frt,nae[i],naia[i])+3*ipumpa[i]+inagap[i]); //Naia
  dx[i+neur*10]=-sna[i]/tdecay+alphan*glut*(1.-sna[i]); //sna
  
}

   return GSL_SUCCESS;
}

/* The jacobian
 * In case we wanna use implict Burslish-Stoer methods
 */
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    double mu = *(double *)params;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;

    // If we use Burlish-Stoer we need to provide the Jacobian
    /*gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
    gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));*/

    // Otherwise, set the Jacobian to 0
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 0, 0.0);
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
  
    return GSL_SUCCESS;
}



