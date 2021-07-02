#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// Change this struct according to colloid.h
typedef struct {

  int index;
  int rebuild;
  int nbonds;
  int nangles;
  int isfixedr;
  int isfixedv;
  int isfixedw;
  int isfixeds;
  int type;
  int bond[2];
  int rng;
  int intpad[20];
  double a0;
  double ah;
  double r[3];
  double v[3];
  double w[3];
  double s[3];
  double m[3];
  double b1;
  double b2;
  double c;
  double h;
  double dr[3];
  double deltaphi;
  double q0;            /* magnitude charge 0 */
  double q1;            /* magnitude charge 1 */
  double epsilon;       /* permittivity */
  double deltaq0;       /* surplus/deficit of charge 0 at change of shape */
  double deltaq1;       /* surplus/deficit of charge 1 at change of shape */
  double sa;            /* surface area (finite difference) */
  double saf;           /* surface area to fluid (finite difference grid) */
  double spare1;
  double spare2[3];
  double tumbletheta;
  double tumblephi;
  double jangle1;
  double jangle2;
  double mu1;
  double mu2;
  double muangle;
  double alpha1;
  double alpha2;
  double alphaangle;
  double dpad[2];

} colloid_state_type;


int main (int argc, char **argv) {

	int ic_g,jc_g,kc_g, index_l,index_g;
	int ETIME, STIME, TSTEP;

  int i,j,k,index;
	struct stat st = {0};

  int it,ib;
  char fname[128];
  FILE *fin, *fout0, *fout, *fout2;

  int Ncoll, Ncoll2;
  colloid_state_type *colloid;

  double *rx, *ry, *rz, *z, *size;
  double *ux, *uy, *uz, *ut;
  double *wx, *wy, *wz;
  double *mx, *my, *mz;
  double *mmx, *mmy, *mmz;
  double *sx, *sy, *sz;
  double *xold,*yold, *zold, *zoldy;
  double angulo, theta, phi;
  double temp;
  double pi = acos(0.0)*2.0;

  if(argc!=4) {
    fprintf(stderr,"Error! Usage is: %s <init time> <end time> <time step> \n",argv[0]);
    exit(2);
  }

  STIME = atoi(argv[1]);
  ETIME = atoi(argv[2]);
  TSTEP = atoi(argv[3]);

  // If folder does not exist, create it!
  if (stat("cortesC", &st) == -1) {
    mkdir("cortesC", 0777);
  }

  //First, I look for how many particles there are in the first file
  sprintf(fname,"config.cds%.8d.001-001",STIME);
  fin = fopen(fname,"r");
  if(fin==NULL) fprintf(stderr,"Error! File %s not found!\n",fname);
  fread(&Ncoll,sizeof(int),1,fin);
  fclose(fin);

	//Then, I get the RAM needed
  ut = (double *)malloc(sizeof(double)*Ncoll);
  z = (double *)malloc(sizeof(double)*Ncoll);
  ut = (double *)malloc(sizeof(double)*Ncoll);
  mmx = (double *)malloc(sizeof(double)*Ncoll);
  mmy = (double *)malloc(sizeof(double)*Ncoll);
  mmz = (double *)malloc(sizeof(double)*Ncoll);
  xold = (double *)malloc(sizeof(double)*Ncoll);
  yold = (double *)malloc(sizeof(double)*Ncoll);
  zold = (double *)malloc(sizeof(double)*Ncoll);
  zoldy = (double *)malloc(sizeof(double)*Ncoll);

  // Positions
  rx = (double *)malloc(sizeof(double)*Ncoll);
  ry = (double *)malloc(sizeof(double)*Ncoll);
  rz = (double *)malloc(sizeof(double)*Ncoll);
  size = (double *)malloc(sizeof(double)*Ncoll);

  // Velocity
  ux = (double *)malloc(sizeof(double)*Ncoll);
  uy = (double *)malloc(sizeof(double)*Ncoll);
  uz = (double *)malloc(sizeof(double)*Ncoll);

  // Velocity (angular)
  wx = (double *)malloc(sizeof(double)*Ncoll);
  wy = (double *)malloc(sizeof(double)*Ncoll);
  wz = (double *)malloc(sizeof(double)*Ncoll);

  // Direction
  mx = (double *)malloc(sizeof(double)*Ncoll);
  my = (double *)malloc(sizeof(double)*Ncoll);
  mz = (double *)malloc(sizeof(double)*Ncoll);

  // s vector. (?)
  sx = (double *)malloc(sizeof(double)*Ncoll);
  sy = (double *)malloc(sizeof(double)*Ncoll);
  sz = (double *)malloc(sizeof(double)*Ncoll);

  for(it=STIME;it<=ETIME;it+=TSTEP) {

    sprintf(fname,"config.cds%.8d.001-001",it);
    fin = fopen(fname,"r");
    if(fin==NULL) fprintf(stderr,"Error! File %s not found!\n",fname);
    fread(&Ncoll2,sizeof(int),1,fin);

    if (Ncoll2!=Ncoll){
      printf("WARNING! We expected %d particles, but we found %d at step %d\n",Ncoll, Ncoll2, it);
    }

    colloid=(colloid_state_type *)malloc(sizeof(colloid_state_type)*Ncoll);
    fread(colloid,sizeof(colloid_state_type),Ncoll,fin);
    fclose(fin);

    for(ib=0;ib<Ncoll;ib++) {
      if (it==STIME){
        printf("Wetting for %d particle. H %f	C	%f\n", ib, colloid[ib].h, colloid[ib].c);
      }

      if (it==STIME){
        xold[ib]=0.0;
        yold[ib]=0.0;
        zold[ib]=0.0;
        zoldy[ib] = 0.0;
      }else{
        xold[ib]=colloid[ib].r[0]-rx[ib];
        yold[ib]=colloid[ib].r[1]-ry[ib];
        zold[ib]=colloid[ib].r[2]/colloid[ib].ah-rz[ib]/colloid[ib].ah;
        zoldy[ib]=acos(colloid[ib].m[2])/pi-acos(mz[ib])/pi;
        //zoldy[ib] = zold[ib];

        //Normaliza

        //xold[ib]/=sqrt(xold[ib]*xold[ib]+zold[ib]*zold[ib]);
        //yold[ib]/=sqrt(yold[ib]*yold[ib]+zold[ib]*zold[ib]);
        //zold[ib]/=sqrt(zold[ib]*zold[ib]+zoldy[ib]*zoldy[ib]);
        //zoldy[ib]/=sqrt(zold[ib]*zold[ib]+zoldy[ib]*zoldy[ib]);
      }

      rx[ib] = colloid[ib].r[0];
      ry[ib] = colloid[ib].r[1];
      rz[ib] = colloid[ib].r[2];

      ux[ib] = colloid[ib].v[0];
      uy[ib] = colloid[ib].v[1];
      uz[ib] = colloid[ib].v[2];

      wx[ib] = colloid[ib].w[0];
      wy[ib] = colloid[ib].w[1];
      wz[ib] = colloid[ib].w[2];

      mx[ib] = colloid[ib].m[0];
      my[ib] = colloid[ib].m[1];
      mz[ib] = colloid[ib].m[2];

      sx[ib] = colloid[ib].s[0];
      sy[ib] = colloid[ib].s[1];
      sz[ib] = colloid[ib].s[2];

      size[ib] = colloid[ib].a0;

      if (it==STIME){
        ut[ib] = 0.0;
        z[ib] = 0.0;
        mmx[ib] = 0.0;
        mmy[ib] = 0.0;
        mmz[ib] = 0.0;
      }

    }

    for(ib=0;ib<Ncoll;ib++) {
      sprintf(fname,"cortesC/pos+vel_%d.dat",ib);					//LUCAS
      fout = fopen(fname,"a+");

      sprintf(fname,"cortesC/pos+vel_grafica_%d.dat",ib);					//LUCAS
      fout2 = fopen(fname,"a+");

      if (ux[ib]>0.1||uy[ib]>0.1||uz[ib]>0.1) {
        printf("ERROR!!!! velocidad m�s alta!!! %d %e %e %e", ib, ux[ib],uy[ib],uz[ib]);
      }else{

        angulo = 0.0;

        if ((mz[ib]>0)&&(mx[ib]>=0)) {
          angulo = 180 - acos(mz[ib])*180/pi;
        } else if ((mz[ib]>=0)&&(mx[ib]<0)) {
          angulo = 180 + acos(mz[ib])*180/pi;
        } else if ((mz[ib]<=0)&&(mx[ib]>0)) {
          angulo = acos(-1*mz[ib])*180/pi;
        } else if ((mz[ib]<0)&&(mx[ib]<=0)) {
          angulo =  -1*acos(-1*mz[ib])*180/pi;
        } else {
          printf("ERROR\n");
        }

        //theta
        theta = atan2(hypot(mx[ib],my[ib]), mz[ib])*180/M_PI;

        //phi
        phi = atan2(my[ib],mx[ib])*180/M_PI;

        fprintf(fout," %d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",ib, it, rx[ib],ry[ib],rz[ib],
	                                                ux[ib],uy[ib],uz[ib],
                                                        wx[ib],wy[ib],wz[ib],
                                                        mx[ib],my[ib],mz[ib],
                                                        angulo, phi, theta);
        if (it>30000) {
          ut[ib] += sqrt(pow(ux[ib],2)+pow(uy[ib],2)+pow(uz[ib],2));
        }

        if (it>30000) {
          z[ib] += rz[ib];
          mmx[ib] += mx[ib];
          mmy[ib] += my[ib];
          mmz[ib] += mz[ib];
        }

        temp = colloid[ib].mu2*colloid[ib].alpha2/(0.0625*2);

        fprintf(fout2,"%d %d %e %e %e %e %e\n",ib, it, rz[ib]/colloid[ib].ah, acos(mz[ib])/pi, uy[ib]/temp,zold[ib],zoldy[ib]);

      }
      fclose(fout);
      fclose(fout2);
    }

  	free(colloid);

  } /* time loop */

  return 0;

}
