#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#define val 800000  /* this can be changed */

typedef struct {
	double x,y,z;
} trivec;

int main (int argc, char **argv) {
	int NX,NY,NZ,NTOT;
	int PEx,PEy,PEz;
	int lnodes,nox,noy,noz;
	int ic_g,jc_g,kc_g, index_l,index_g;
	int nlocal[3],zone;
	int ETIME, STIME, TSTEP;

	int i,j,k,index;
	struct stat st = {0};

	int it;
	char fname[128];
	FILE *fin, *fout;

	static double vx[val], vy[val], vz[val];	// Change in case you need more memory!
	trivec *vel;

	// parameters needed:
	// SIZE BOX, DECOMPOSITION, TIME, SLICE
	// NX NY NZ PEx PEy PEz STIME ETIME TSTEP SLICE
	NX = atoi(argv[1]);
	NY = atoi(argv[2]);
	NZ = atoi(argv[3]);
	PEx=atoi(argv[4]);
	PEy=atoi(argv[5]);
	PEz=atoi(argv[6]);
	STIME = atoi(argv[7]);
	ETIME = atoi(argv[8]);
	TSTEP = atoi(argv[9]);

	nlocal[0] = NX/PEx;
	nlocal[1] = NY/PEy;
	nlocal[2] = NZ/PEz;
 	NTOT = NX*NY*NZ;
	printf("Num of nodes = %d \n", NTOT);
	lnodes=nlocal[0]*nlocal[1]*nlocal[2];

	// We will print for Y slice given at argv[10] and save inside a particular folder
	// if that folder does not exist, create it!
	sprintf(fname,"cortes%d",atoi(argv[10]));
	if (stat(fname, &st) == -1) {
		mkdir(fname, 0777);
	}

	vel=(trivec *)malloc(sizeof(trivec)*NTOT);

	// Main loop for each velocity file printed
	for(it=STIME;it<=ETIME;it+=TSTEP) {
		// Check that the file exists first
		sprintf(fname,"vel-%.8d.001-001",it);	// it can happen that you need less zeros... change the 8 then!
		fin=fopen(fname,"rb");
		if(fin==NULL) {
			fprintf(stderr,"File %s not found!\n",fname);
		} else {
			fread(vel,sizeof(trivec),NTOT,fin);
			fclose(fin);
		} /* condition on finding file */

		// Extract data from the binary field following the parallelization from Ludwig.
		zone=0;
		for(nox=0;nox<PEx;nox++){
			for(noy=0;noy<PEy;noy++){
				for(noz=0;noz<PEz;noz++){
					zone+=1;
					for(i=0;i<nlocal[0];i++) {
						ic_g = nox*nlocal[0] + i;
						for(j=0;j<nlocal[1];j++) {
							jc_g = noy*nlocal[1] + j;
							for(k=0;k<nlocal[2];k++) {
								kc_g = noz*nlocal[2] + k;
								index_l = k + j*nlocal[2] + i*nlocal[1]*nlocal[2];
								index_g = kc_g + jc_g*NZ + ic_g*NY*NZ;
								index = index_l+(zone-1)*lnodes;
								vx[index_g]=vel[index].x;
								vy[index_g]=vel[index].y;
								vz[index_g]=vel[index].z;
							}
						}
					}
				}
			}
		}

		// Print only XZ components of velocity to a file inside folder cortesYSLICENUMBER
		sprintf(fname,"cortes%d/vel%d",atoi(argv[10]),it);		//LUCAS
		fout=fopen(fname,"w");

		j=atoi(argv[10]);
		for(i=0;i<NX;i++) {
			for(k=0;k<NZ;k++) {
				index = k + j*NZ + i*NY*NZ;
				fprintf(fout," %d %d %e %e \n",i,k,vx[index],vz[index]);
			}
 		}

		fclose(fout);

	} /* time loop */

	return 0;

}
