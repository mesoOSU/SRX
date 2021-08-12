#include <stdio.h>
#include <mpi.h>
#include "evp.h"

void PrintTensor(ten2nd m)
{
	int i,j;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++)
			printf("%e\t",m[i][j]);
		printf("\n");
	}
	return;
}/*end PrintTensor()*/

void PrintTensorNorm(ten2nd m)
{
	real norm = 0.;
	T2_loop{
		norm += m[mi][mj]*m[mi][mj];
	}
	norm = sqrt(norm);
	printf("pyz: Tensor norm = %e\n",norm);

}/*end PrintTensorNorm()*/

void PrintVoigt(voigt m)
{
	int i;
	for(i=0;i<6;i++){
		printf("%.3e\t",m[i]);
	}
	printf("\n");
	return;
}/*end PrintVoigt()*/

void PrintVoigt66(voigt66 m)
{
	int i,j;
	for(i=0;i<6;i++){
		for(j=0;j<6;j++)
			printf("%e\t",m[i][j]);
		printf("\n");
	}
	return;
}/*end PrintVoigt66()*/

void PrintB(void)
{
	int m;
	for(m=0;m<6;m++){
		printf("B[%d]:\n",m+1);
		PrintTensor(B[m]);
	}

	return;
}/*end PrintB()*/

void WriteEtaMPI(char *s)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;

	AllocMem(tmpVector, lsize, real);
	
	local_loop{
		tmpVector[pIDX] = (real)grain_f[pIDX];
	}

	sprintf(fname, "%s.iout", s);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, mpirank*lsize*sizeof(double), MPI_DOUBLE,
			MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}/*end WriteEtaMPI()*/

void PrintSchmidXt(void)
{
	int i,j;
	for(j=0;j<N_modes;j++){
		for(i=0;i<nsmx;i++){
			printf("pyz: Schmid vector of mode#%d:",i);
			printf("%.3f %.3f %.3f %.3f %.3f\n",
					Schm_xt[j][i][0], Schm_xt[j][i][1], Schm_xt[j][i][2], Schm_xt[j][i][3], Schm_xt[j][i][4]);
		}
	}
	return;
}/*end PrintSchmidXt()*/

static void GenerateLamellae(char *s, real vol, real t1[],real ph[],real t2[])
{
	FILE *fp;
	char fname[80];
	int i,j,k;
	int jph,jgr;
	int N=64;
	real half_d; 

	half_d = sqrt(2)/2.*N*(1.-sqrt(1-vol));
	sprintf(fname,"%s.ms",s);
	fp=fopen(fname,"w");
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			for(k=0;k<N;k++){
				if((k>=N-sqrt(2)*half_d-j)&&(k<N+sqrt(2)*half_d-j)){
					jph = 2;
					jgr = 1;
					fprintf(fp,"%2f\t%2f\t%2f\t%d\t%d\t%d\t%d\t%d\n",
							t1[jph-1],ph[jph-1],t2[jph-1],i+1,j+1,k+1,jgr,jph);
				}
				else{
					jph = 1;
					jgr = 0;
					fprintf(fp,"%2f\t%2f\t%2f\t%d\t%d\t%d\t%d\t%d\n",
							t1[jph-1],ph[jph-1],t2[jph-1],i+1,j+1,k+1,jgr,jph);
				}
			}
		}
	}
	fclose(fp);

	return;
}/*end GenerateLamellae()*/

void Ti_alpha_beta(void)
{
	/* calculate the Euler angle for seven samples of
	   single alpha/beta colony */

	ten2nd e_xt;	// xtal directions in computational basis (sample axes)
					// for hcp, we use [2 -1 -1 0] (x-axis), [0 1 -1 0] (y-axis), and [0 0 0 1] (z-axis)
	ten2nd tran_m;	// sample -> xtal(alpha)
	ten2nd tran_BurgerOR;	// xtal(alph) -> xtal(beta)
	ten2nd tran_tot;
	real t1, ph, t2;
	real et1[2],eph[2],et2[2];
	int jph;
	ten2nd e_sa = {{1.0,0.,0.},{0.,1.0,0.},{0.,0.,1.}};	// sample axes
	int i,j,k;
	real theta;

	/* transformation matrix */
	theta = atan(1./sqrt(2));
	e_xt[0][0]=cos(theta)/sqrt(2); e_xt[0][1]=-1.*sin(theta)/sqrt(2); e_xt[0][2]=1./sqrt(2);
	e_xt[1][0]=-1.*cos(theta)/sqrt(2); e_xt[1][1]=sin(theta)/sqrt(2); e_xt[1][2]=1./sqrt(2);
	e_xt[2][0]=-1.*sin(theta); e_xt[2][1]=-1.*cos(theta); e_xt[2][2]=0.0;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_BurgerOR[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}


	/* prismatic a1 */
	e_xt[0][0]=0.0; e_xt[0][1]=1./sqrt(2); e_xt[0][2]=-1./sqrt(2);
	e_xt[1][0]=0.0; e_xt[1][1]=1./sqrt(2); e_xt[1][2]=1./sqrt(2);
	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=0.0;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Prismatic a1 of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Prismatic a1 of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(mpirank==0){
		GenerateLamellae("Prism_1_64x64x64",0.12,et1,eph,et2);
	}

	/* prismatic a2 */
	theta = 75.*PI/180.;
	e_xt[0][0]=0.0; e_xt[0][1]=cos(theta); e_xt[0][2]=-1.*sin(theta);
	e_xt[1][0]=0.0; e_xt[1][1]=sin(theta); e_xt[1][2]=cos(theta);
	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=0.0;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Prismatic a2 of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Prismatic a2 of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(mpirank==0){
		GenerateLamellae("Prism_2_64x64x64",0.12,et1,eph,et2);
	}

	/* prismatic a3 */
	theta = 75.*PI/180.;
	e_xt[0][0]=0.0; e_xt[0][2]=-1.*cos(theta); e_xt[0][1]=sin(theta);
	e_xt[1][0]=0.0; e_xt[1][2]=sin(theta); e_xt[1][1]=cos(theta);
	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=0.0;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Prismatic a3 of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Prismatic a3 of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(mpirank==0){
		GenerateLamellae("Prism_3_64x64x64",0.12,et1,eph,et2);
	}

	/* basal a1 */
	e_xt[0][0]=0.0; e_xt[0][1]=1.0/sqrt(2); e_xt[0][2]=-1.0/sqrt(2);
	e_xt[1][0]=1.0; e_xt[1][1]=0.; e_xt[1][2]=0.;
	e_xt[2][0]=0.0; e_xt[2][1]=1.0/sqrt(2); e_xt[2][2]=1.0/sqrt(2);
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Basal a1 of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Basal a1 of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(mpirank==0){
		GenerateLamellae("Basal_1_64x64x64",0.12,et1,eph,et2);
	}

	/* do a rotation about xtal c axis to obtain proper coord.
	   of basal a2 and a3 */
	/* basal a2 */
	e_xt[0][0]=0.866025; e_xt[0][1]=-0.353553; e_xt[0][2]=0.353553;
	e_xt[1][0]=0.5; e_xt[1][1]=0.612372; e_xt[1][2]=-0.612372;
	e_xt[2][0]=0.0; e_xt[2][1]=1.0/sqrt(2); e_xt[2][2]=1.0/sqrt(2);
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Basal a2 of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Basal a2 of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(mpirank==0){
		GenerateLamellae("Basal_2_64x64x64",0.12,et1,eph,et2);
	}

	/* basal a3 */
	e_xt[0][0]=-0.866025; e_xt[0][1]=-0.353553; e_xt[0][2]=0.353553;
	e_xt[1][0]=0.5; e_xt[1][1]=-0.612372; e_xt[1][2]=0.612372;
	e_xt[2][0]=0.0; e_xt[2][1]=1.0/sqrt(2); e_xt[2][2]=1.0/sqrt(2);
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Basal a3 of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Basal a3 of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(mpirank==0){
		GenerateLamellae("Basal_3_64x64x64",0.12,et1,eph,et2);
	}

	/* pyramidal c+a */
	e_xt[0][0]=0.; e_xt[0][1]=1.; e_xt[0][2]=0.;
	e_xt[1][0]=-1.; e_xt[1][1]=0.; e_xt[1][2]=0.;
	e_xt[2][0]=0.; e_xt[2][1]=0.; e_xt[2][2]=1.;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	jph=1;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Pyramidal c+a of alpha:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_tot[i][j] = 0.0;
			for(k=0;k<3;k++)
				tran_tot[i][j] += tran_BurgerOR[i][k]*tran_m[k][j];
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_tot,1);
	jph=2;
	et1[jph-1] = t1; eph[jph-1]=ph; et2[jph-1]=t2;
	if(mpirank==0){
		printf("Pyramidal c+a of beta:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(mpirank==0){
		GenerateLamellae("Pyramidal_64x64x64",0.12,et1,eph,et2);
	}

	return;
}/*end Ti_alpha_beta()*/

static void GenerateAlSingleXtal(char *s, real t1, real ph, real t2, int Nx, int Ny, int Nz)
{
	FILE *fp;
	char fname[80];
	int i,j,k;
	int jph,jgr;

	sprintf(fname,"%s.ms",s);
	fp=fopen(fname,"w");
	jph=1, jgr=0;
	for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			for(k=0;k<Nz;k++){
				fprintf(fp,"%2f\t%2f\t%2f\t%d\t%d\t%d\t%d\t%d\n",
						t1,ph,t2,i+1,j+1,k+1,jgr,jph);
			}
		}
	}
	fclose(fp);

	return;
}/*end GenerateAlSingleXtal()*/

void SingleXtal_Al(void)
{
	/* Generate a Al single crystal with <110> axis parallel to the
	   compression axis (z-axis in our simulation) */

	ten2nd e_xt;	// xtal directions in computational basis (sample axes)
					// for fcc, we use [0 0 1] (x-axis), [1 -1 0] (y-axis), and [1 1 0] (z-axis)
	ten2nd tran_m;	// sample -> xtal(alpha)
	real t1, ph, t2;
	ten2nd e_sa = {{1.0,0.,0.},{0.,1.0,0.},{0.,0.,1.}};	// sample axes
	int i,j,k;

	/* transformation matrix */
	e_xt[0][0]=0.0; e_xt[0][1]=0.0; e_xt[0][2]=1.0;
	e_xt[1][0]=1./sqrt(2); e_xt[1][1]=-1./sqrt(2); e_xt[1][2]=0.0;
	e_xt[2][0]=1./sqrt(2); e_xt[2][1]=1./sqrt(2); e_xt[2][2]=0.0;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1,&ph,&t2,tran_m,1);
	//t1=0.0; ph=45.0; t2=0.0;
	t1=0.0; ph=0.0; t2=0.0;
	if(mpirank==0){
		printf("Al single crytal:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1,ph,t2);
	}
	if(mpirank==0){
		//GenerateAlSingleXtal("AlSingleXtal_32x32x32",t1,ph,t2,32,32,32);
		//GenerateAlSingleXtal("AlSingleXtal_16x32x32",t1,ph,t2,16,32,32);
		GenerateAlSingleXtal("SphericalInclusion_128",t1,ph,t2,128,128,128);
	}

	return;
}/*end SingleXtal_Al()*/

void Bicrystal_Ti_2(int nx, int ny, int nz)
{
	/* Generate a BCC-Ti bicrystal, using Exp. Euler angle*/

	//ten2nd e_xt;	// xtal directions in computational basis (sample axes)
	//ten2nd tran_m;	// sample -> xtal(alpha)
	real t1[2], ph[2], t2[2];
	//ten2nd e_sa = {{1.0,0.,0.},{0.,1.0,0.},{0.,0.,1.}};	// sample axes
	int i,j,k;
	//real norm;
	FILE *fp;

//	/* Grain#1 */
//	/* transformation matrix */
//	e_xt[0][0]=0.; e_xt[0][1]=1.0; e_xt[0][2]=0.;
//	e_xt[1][0]=-1.; e_xt[1][1]=0.; e_xt[1][2]=0.;
//	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=1.0;
//	for(i=0;i<3;i++){
//		norm=0.0;
//		for(j=0;j<3;j++){
//			norm += e_xt[i][j]*e_xt[i][j];
//		}
//		norm = sqrt(norm);
//		for(j=0;j<3;j++){
//			e_xt[i][j] /= norm;
//		}
//	}
//	for(i=0;i<3;i++){
//		for(j=0;j<3;j++){
//			tran_m[i][j]=0.0;
//			for(k=0;k<3;k++){
//				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
//			}
//		}
//	}
//	EulerToTransMatrix(&t1[0],&ph[0],&t2[0],tran_m,1);
	t1[0] = 17.2; ph[0] = 130.8; t2[0] = 35.2;
	if(mpirank==0){
		printf("BCC-Ti, Grain#1:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1[0],ph[0],t2[0]);
	}

//	/* Grain#2 */
//	/* transformation matrix */
//	e_xt[0][0]=0.; e_xt[0][1]=1.0; e_xt[0][2]=0.;
//	e_xt[1][0]=-1.; e_xt[1][1]=0.; e_xt[1][2]=0.;
//	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=1.0;
//	for(i=0;i<3;i++){
//		norm=0.0;
//		for(j=0;j<3;j++){
//			norm += e_xt[i][j]*e_xt[i][j];
//		}
//		norm = sqrt(norm);
//		for(j=0;j<3;j++){
//			e_xt[i][j] /= norm;
//		}
//	}
//	for(i=0;i<3;i++){
//		for(j=0;j<3;j++){
//			tran_m[i][j]=0.0;
//			for(k=0;k<3;k++){
//				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
//			}
//		}
//	}
//	EulerToTransMatrix(&t1[1],&ph[1],&t2[1],tran_m,1);
	t1[1] = 228.8; ph[1] = 83.9; t2[1] = 90.1;
	if(mpirank==0){
		printf("BCC-Ti, Grain#2:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1[1],ph[1],t2[1]);
	}


	/* Create ms file */
	fp = fopen("BiXtal_Ti.ms","w");
	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			for(k=0;k<nz;k++){
				if(k<nz/2){
					fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
							t1[0],ph[0],t2[0],i+1,j+1,k+1,0,1);
				}else{
					fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
							t1[1],ph[1],t2[1],i+1,j+1,k+1,1,1);
				}
			}
		}
	}
	fclose(fp);

	return;
}/*end Bicrystal_Ti_2()*/

void Bicrystal_Ti_3(int nx, int ny, int nz)
{
	/* Generate a BCC-Ti bicrystal, symmetric twist, theta = 10.52*/

	ten2nd e_xt;	// xtal directions in computational basis (sample axes)
	ten2nd tran_m;	// sample -> xtal(alpha)
	real t1[2], ph[2], t2[2];
	ten2nd e_sa = {{1.0,0.,0.},{0.,1.0,0.},{0.,0.,1.}};	// sample axes
	int i,j,k;
	real norm;
	FILE *fp, *fp1;

	/* Grain#1 */
	/* transformation matrix */
	e_xt[0][0]=0.0; e_xt[0][1]=1.0; e_xt[0][2]=0;
	e_xt[1][0]=-0.904114; e_xt[1][1]=0.0; e_xt[1][2]=1.08746;
	e_xt[2][0]=1.08746; e_xt[2][1]=0.0; e_xt[2][2]=0.904114;
	for(i=0;i<3;i++){
		norm=0.0;
		for(j=0;j<3;j++){
			norm += e_xt[i][j]*e_xt[i][j];
		}
		norm = sqrt(norm);
		for(j=0;j<3;j++){
			e_xt[i][j] /= norm;
		}
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1[0],&ph[0],&t2[0],tran_m,1);
	if(mpirank==0){
		printf("BCC-Ti, Grain#1:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1[0],ph[0],t2[0]);
	}

	/* Grain#2 */
	/* transformation matrix */
	e_xt[0][0]=0.0; e_xt[0][1]=1.0; e_xt[0][2]=0;
	e_xt[1][0]=-1.08746; e_xt[1][1]=0.0; e_xt[1][2]=0.904114;
	e_xt[2][0]=0.904114; e_xt[2][1]=0.0; e_xt[2][2]=1.08746;
	for(i=0;i<3;i++){
		norm=0.0;
		for(j=0;j<3;j++){
			norm += e_xt[i][j]*e_xt[i][j];
		}
		norm = sqrt(norm);
		for(j=0;j<3;j++){
			e_xt[i][j] /= norm;
		}
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1[1],&ph[1],&t2[1],tran_m,1);
	if(mpirank==0){
		printf("BCC-Ti, Grain#2:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1[1],ph[1],t2[1]);
	}


	/* Create ms file */
	fp = fopen("BiXtal_Ti_Twist.ms","w");
	fp1 = fopen("BiXtal_Ti_Twist.out","w");
	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			for(k=0;k<nz;k++){
				if((i<nx/4)||(i>=(nx-nx/4))){
					fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
							t1[0],ph[0],t2[0],i+1,j+1,k+1,0,1);
					fprintf(fp1,"%d\n",0);
				}else{
					fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
							t1[1],ph[1],t2[1],i+1,j+1,k+1,1,1);
					fprintf(fp1,"%d\n",1);
				}
			}
		}
	}
	fclose(fp);
	fclose(fp1);

	return;
}/*end Bicrystal_Ti_3()*/

void Bicrystal_Ti_diffuse(int nx, int ny, int nz)
{
	/* The GB is treated as a second phase with the average properties */
	/* Generate a BCC-Ti bicrystal, symmetric tilt, theta = 10.52*/

	ten2nd e_xt;	// xtal directions in computational basis (sample axes)
	ten2nd tran_m;	// sample -> xtal(alpha)
	real t1[2], ph[2], t2[2];
	ten2nd e_sa = {{1.0,0.,0.},{0.,1.0,0.},{0.,0.,1.}};	// sample axes
	int i,j,k;
	real norm;
	FILE *fp, *fp1;
	int GB_width = 4;
	real t1_gb, ph_gb, t2_gb;

	/* Grain#1 */
	/* transformation matrix */
	e_xt[0][0]=0.0648243; e_xt[0][1]=0.995789; e_xt[0][2]=-0.0648243;
	e_xt[1][0]=-0.995789; e_xt[1][1]=0.129649; e_xt[1][2]=0.995789;
	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=1.0;
	for(i=0;i<3;i++){
		norm=0.0;
		for(j=0;j<3;j++){
			norm += e_xt[i][j]*e_xt[i][j];
		}
		norm = sqrt(norm);
		for(j=0;j<3;j++){
			e_xt[i][j] /= norm;
		}
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1[0],&ph[0],&t2[0],tran_m,1);
	if(mpirank==0){
		printf("BCC-Ti, Grain#1:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1[0],ph[0],t2[0]);
	}

	/* Grain boundary */
	/* transformation matrix */
	e_xt[0][0]=0.; e_xt[0][1]=1.0; e_xt[0][2]=0.;
	e_xt[1][0]=-1.; e_xt[1][1]=0.; e_xt[1][2]=1.;
	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=1.0;
	for(i=0;i<3;i++){
		norm=0.0;
		for(j=0;j<3;j++){
			norm += e_xt[i][j]*e_xt[i][j];
		}
		norm = sqrt(norm);
		for(j=0;j<3;j++){
			e_xt[i][j] /= norm;
		}
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1_gb,&ph_gb,&t2_gb,tran_m,1);
	if(mpirank==0){
		printf("BCC-Ti, Grain boundary:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1_gb,ph_gb,t2_gb);
	}

	/* Grain#2 */
	/* transformation matrix */
	e_xt[0][0]=-0.0648243; e_xt[0][1]=0.995789; e_xt[0][2]=0.0648243;
	e_xt[1][0]=-0.995789; e_xt[1][1]=-0.129649; e_xt[1][2]=0.995789;
	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=1.0;
	for(i=0;i<3;i++){
		norm=0.0;
		for(j=0;j<3;j++){
			norm += e_xt[i][j]*e_xt[i][j];
		}
		norm = sqrt(norm);
		for(j=0;j<3;j++){
			e_xt[i][j] /= norm;
		}
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1[1],&ph[1],&t2[1],tran_m,1);
	if(mpirank==0){
		printf("BCC-Ti, Grain#2:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1[1],ph[1],t2[1]);
	}

	t1[0] = 90.0; ph[0] = 45.0; t2[0] = -50.0;
	t1[1] = 90.0; ph[1] = 45.0; t2[1] = 50.0;
	

	/* Create ms file */
	fp = fopen("BiXtal_Ti_diffuse.ms","w");
	fp1 = fopen("BiXtal_Ti_diffuse.out","w");
	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			for(k=0;k<nz;k++){
				if(((i<GB_width/2)||((i>=(nx/2-GB_width/2))&&(i<(nx/2+GB_width/2))))||(i>=nx-GB_width/2)){
					fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
							t1_gb,ph_gb,t2_gb,i+1,j+1,k+1,3,2);
					fprintf(fp1,"%d\n",3);
				}
				else if(i<nx/2){
					fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
							t1[0],ph[0],t2[0],i+1,j+1,k+1,0,1);
					fprintf(fp1,"%d\n",0);
				}
				else{
					fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
							t1[1],ph[1],t2[1],i+1,j+1,k+1,1,1);
					fprintf(fp1,"%d\n",1);
				}
			}
		}
	}
	fclose(fp);
	fclose(fp1);

	return;
}/*end Bicrystal_Ti_diffuse()*/
void Bicrystal_Ti(int nx, int ny, int nz)
{
	/* Generate a BCC-Ti bicrystal, symmetric tilt, theta = 10.52*/

	ten2nd e_xt;	// xtal directions in computational basis (sample axes)
	ten2nd tran_m;	// sample -> xtal(alpha)
	real t1[2], ph[2], t2[2];
	ten2nd e_sa = {{1.0,0.,0.},{0.,1.0,0.},{0.,0.,1.}};	// sample axes
	int i,j,k;
	real norm;
	FILE *fp, *fp1;

	/* Grain#1 */
	/* transformation matrix */
//	e_xt[0][0]=0.0648243; e_xt[0][1]=0.995789; e_xt[0][2]=-0.0648243;
//	e_xt[1][0]=-0.995789; e_xt[1][1]=0.129649; e_xt[1][2]=0.995789;
//	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=1.0;
	e_xt[0][0]=0.2418; e_xt[0][1]=0.9397; e_xt[0][2]=-0.2418;
	e_xt[1][0]=-0.9397; e_xt[1][1]=0.4837; e_xt[1][2]=0.9397;
	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=1.0;
	for(i=0;i<3;i++){
		norm=0.0;
		for(j=0;j<3;j++){
			norm += e_xt[i][j]*e_xt[i][j];
		}
		norm = sqrt(norm);
		for(j=0;j<3;j++){
			e_xt[i][j] /= norm;
		}
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1[0],&ph[0],&t2[0],tran_m,1);
	if(mpirank==0){
		printf("BCC-Ti, Grain#1:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1[0],ph[0],t2[0]);
	}

	/* Grain#2 */
	/* transformation matrix */
//	e_xt[0][0]=-0.0648243; e_xt[0][1]=0.995789; e_xt[0][2]=0.0648243;
//	e_xt[1][0]=-0.995789; e_xt[1][1]=-0.129649; e_xt[1][2]=0.995789;
//	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=1.0;
	e_xt[0][0]=-0.2418; e_xt[0][1]=0.9397; e_xt[0][2]=0.2418;
	e_xt[1][0]=-0.9397; e_xt[1][1]=-0.4837; e_xt[1][2]=0.9397;
	e_xt[2][0]=1.0; e_xt[2][1]=0.0; e_xt[2][2]=1.0;
	for(i=0;i<3;i++){
		norm=0.0;
		for(j=0;j<3;j++){
			norm += e_xt[i][j]*e_xt[i][j];
		}
		norm = sqrt(norm);
		for(j=0;j<3;j++){
			e_xt[i][j] /= norm;
		}
	}
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_m[i][j]=0.0;
			for(k=0;k<3;k++){
				tran_m[i][j] += e_xt[i][k]*e_sa[j][k];
			}
		}
	}
	EulerToTransMatrix(&t1[1],&ph[1],&t2[1],tran_m,1);
	if(mpirank==0){
		printf("BCC-Ti, Grain#2:\n");
		printf("Euler angle = (%e,%e,%e)\n",t1[1],ph[1],t2[1]);
	}


//	t1[0] = 0.0; ph[0] = 0.0; t2[0] = 0.0;
	/* Create ms file */
	fp = fopen("BiXtal_Ti.ms","w");
	fp1 = fopen("BiXtal_Ti.out","w");
	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			for(k=0;k<nz;k++){
//				if((i<nx/4)||(i>=(nx-nx/4))){
//					fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
//							t1[0],ph[0],t2[0],i+1,j+1,k+1,0,1);
//					fprintf(fp1,"%d\n",0);
//				}else{
//					fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
//							t1[1],ph[1],t2[1],i+1,j+1,k+1,1,1);
//					fprintf(fp1,"%d\n",1);
//				}
				if(i<nx/2){
					fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
							t1[0],ph[0],t2[0],i+1,j+1,k+1,0,1);
					fprintf(fp1,"%d\n",0);
				}else{
					fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
							t1[1],ph[1],t2[1],i+1,j+1,k+1,1,1);
					fprintf(fp1,"%d\n",1);
				}
			}
		}
	}
	fclose(fp);
	fclose(fp1);

	return;
}/*end Bicrystal_Ti()*/

void Gamma_GammaPrime(int nx, int ny, int nz, real VolFrac)
{
	int i, j, k;
	int phase;
	int channel_xs, channel_xe;
	int channel_ys, channel_ye;
	int channel_zs, channel_ze;
  real WithFrac;  // the ratio between gamma-pirme cube edge and whole simulateion box edge
	FILE *fp;
	FILE *fp1;
  
  WithFrac = 1.0-2.0*pow(VolFrac/8.0, 1.0/3);

	channel_xs = nx/2 - (int)(WithFrac*nx)/2;
	channel_xe = nx/2 + (int)(WithFrac*nx)/2;
	channel_ys = nx/2 - (int)(WithFrac*ny)/2;
	channel_ye = nx/2 + (int)(WithFrac*ny)/2;
	channel_zs = nx/2 - (int)(WithFrac*nz)/2;
	channel_ze = nx/2 + (int)(WithFrac*nz)/2;

	fp = fopen("NiAlloy_test2.ms","w");
	fp1 = fopen("NiAlloyVTK.out","w");
	for(i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			for(k=0;k<nz;k++){
				if((i>=channel_xs&&i<channel_xe)||(j>=channel_ys&&j<channel_ye)||(k>=channel_zs&&k<channel_ze)){
					phase = 1;	// gamma matrix
				}else{
					phase = 2;	// gamma prime
				}
				fprintf(fp,"%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
						0.,45.,0.,i+1,j+1,k+1,0,phase);
				fprintf(fp1,"%d\n",phase);
			}
		}
	}
	fclose(fp);
	
	return;
}/*end Gamma_GammaPrime()*/

#ifdef PF_DRX
void SphericalInclusion()
{
	real strain = 0.01;
	real factor;

	real r0 = 0.1;	// radius
	real h0 = 0.5;
	real rx = CellDim[0]*r0*2*h0;	// the box is scaled with 2*h0;
	real ry = CellDim[1]*r0*2*h0;
	real rz = CellDim[2]*r0*2*h0;

	real x0 = 0.5*CellDim[0]*2*h0;
	real y0 = 0.5*CellDim[1]*2*h0;
	real z0 = 0.5*CellDim[2]*2*h0;

	// isotropic elasticity
	real mu = 1000.0;
	real nu = 0.3;
	real lambda = 2.0*nu/(1.0-2.0*nu)*mu;

	ten4th c0ijkl = {0.0};

	ten2nd epsavg = {{1.0/3.0, 1.0/2.0, 0.0},{1.0/2.0, 1.0/3.0, 0.0}, {0.0, 0.0, 1.0/3.0}};

	real check=0.0;
	C4_loop{
		c0ijkl[mi][mj][mk][ml] = lambda*(mi==mj)*(mk==ml)+mu*((mi==mk)*(mj==ml)+(mi==ml)*(mj==mk));
		check += fabs(c0ijkl[mi][mj][mk][ml]-C0[mi][mj][mk][ml]);
	}
	printf("pyz:the difference between c0ijkl and C0 is %lf\n",check);

	T2_loop{
		epsavg[mi][mj] *= strain;
	}

	local_loop{
		factor = pow((px+lxs-x0),2.0)/rx/rx + pow((py-y0),2.0)/ry/ry + pow((pz-z0),2.0)/rz/rz;
		factor = 0.5+tanh((1.0-factor)*1000)/2.0;
		T2_loop{
			Eps[pIDX][mi][mj] = factor*epsavg[mi][mj];
		}
	}

	T2_loop{
		epsavg[mi][mj] = 0.0;
	}
	WriteEpsMPI("eps_initial",0);

	//HomogeneousStressSolver(c0ijkl, Eps, epsavg,
	//	DisGrad, Sig);
	InhomogeneousStressSolver(C_gr, Eps, epsavg,
			c0ijkl, Eps0_last,
			Tau_0, Tau,
			DisGrad, Sig);

	WriteSigMPI("sig_final",0);

	return;
}/*end SphericalInclusion()*/

void SphericalVoid()
{
	real factor;

	real r0 = 0.1;	// radius
	real h0 = 0.5;
	real rx = CellDim[0]*r0*2*h0;	// the box is scaled with 2*h0;
	real ry = CellDim[1]*r0*2*h0;
	real rz = CellDim[2]*r0*2*h0;

	real x0 = 0.5*CellDim[0]*2*h0;
	real y0 = 0.5*CellDim[1]*2*h0;
	real z0 = 0.5*CellDim[2]*2*h0;

	// isotropic elasticity
	real mu = 1000.0;
	real nu = 0.3;
	real lambda = 2.0*nu/(1.0-2.0*nu)*mu;

	ten4th c0ijkl = {0.0};
	C4_loop{
		c0ijkl[mi][mj][mk][ml] = lambda*(mi==mj)*(mk==ml)+mu*((mi==mk)*(mj==ml)+(mi==ml)*(mj==mk));
	}

	ten2nd epsavg = {{-nu, 0.0, 0.0},{0.0, -nu, 0.0}, {0.0, 0.0, 1.0}};
	T2_loop{
		epsavg[mi][mj] *= 1.0/2./mu/(1.0+nu);
	}

	local_loop{
		factor = pow((px+lxs-x0),2.0)/rx/rx + pow((py-y0),2.0)/ry/ry + pow((pz-z0),2.0)/rz/rz;
		factor = 0.5+tanh((1.0-factor)*1000)/2.0;
		C6_loop{
			C_gr[pIDX][mi][mj] *= (1.0-factor);
		}
		T2_loop{
			Eps[pIDX][mi][mj] = 0.0;
		}
	}

	InhomogeneousStressSolver(C_gr, Eps, epsavg,
			c0ijkl, Eps0_last,
			Tau_0, Tau,
			DisGrad, Sig);

	WriteSigMPI("sig_final",0);

	return;
}/*end SphericalVoid()*/
#endif
