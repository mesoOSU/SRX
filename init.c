#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fftw_mpi.h>
#include <assert.h>
#include "evp.h"

#define CC(i,j,k,l) (C_ijkl[i-1][j-1][k-1][l-1])

static void ElastStiffnessMatrix(real *ElConst, ten4th cijkl)
{
	real c11, c12, c44;
	real C_ijkl[3][3][3][3] = {0.0};
	int i,j,k,l;

	c11 = *ElConst;
	c12 = *(ElConst+1);
	c44 = *(ElConst+2);

	CC(1,1,1,1) = c11;
	CC(2,2,2,2) = c11;
	CC(3,3,3,3) = c11;
	CC(1,1,2,2) = c12;
	CC(2,2,3,3) = c12;
	CC(3,3,1,1) = c12;
	CC(1,2,1,2) = c44;
	CC(2,3,2,3) = c44;
	CC(3,1,3,1) = c44;
	
	CC(2,2,1,1) = c12;
	CC(3,3,2,2) = c12;
	CC(1,1,3,3) = c12;
	CC(2,1,2,1) = c44;
	CC(1,2,2,1) = c44;
	CC(2,1,1,2) = c44;
	CC(3,2,3,2) = c44;
	CC(3,2,2,3) = c44;
	CC(2,3,3,2) = c44;
	CC(1,3,1,3) = c44;
	CC(1,3,3,1) = c44;
	CC(3,1,1,3) = c44;

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			for(k=0;k<3;k++)
				for(l=0;l<3;l++)
					cijkl[i][j][k][l] = CC(i+1,j+1,k+1,l+1);

	return;
}/*end ElastStiffnessMatrix()*/

static void PlasticInit(char *s, int phase)
{
	FILE *fp;
	char buffer[80] = {0};
	int i, j, mx, m_count;
	int matcheck;
	real norm;
	voigt5 aux5;
	ten2nd aux33;
	voigt66 aux66;
	ten4th aux3333;


	fp = fopen(s, "r");

	fgets(buffer, 80, fp);

	fgets(buffer, 80, fp);
	sscanf(buffer,"%s",XtalSys);

	fgets(buffer, 80, fp);
	fgets(buffer, 80, fp);

	fgets(buffer, 80, fp);
	sscanf(buffer,"%lf %lf %lf",&XtalAxis[0], &XtalAxis[1], &XtalAxis[2]);
	
	fgets(buffer, 80, fp);
	fgets(buffer, 80, fp);
	
	fgets(buffer, 80, fp);
	sscanf(buffer,"%d", &N_modes_max);

	fgets(buffer, 80, fp);
	fgets(buffer, 80, fp);
	
	fgets(buffer, 80, fp);
	sscanf(buffer,"%d", &N_modes);
	AllocMem(iMode,N_modes,int);

	fgets(buffer, 80, fp);
	fgets(buffer, 80, fp);
	for(i=0;i<N_modes;i++){
		fgets(buffer, 80, fp);
		sscanf(buffer,"%d",&iMode[i]);
	}


	fgets(buffer, 80, fp);

	m_count = 0;
	for(i=0;i<N_modes_max;i++){
		fgets(buffer, 80, fp);
		sscanf(buffer, "%*s %*s %d %*s",&mx);
		matcheck = 0;
		j = 0;
		while(j<N_modes){
			if(mx==iMode[j]){
				matcheck = 1;
				m_count++;
				break;
			}
			j++;
		}
			
			fgets(buffer,80,fp);
			fgets(buffer,80,fp);
			sscanf(buffer,"%d\t\t%lf\t\t%lf\t\t%lf\t\t%d", &nsmx, &nrsx, &gamd0x, &twshx, &isectwx);

			fgets(buffer,80,fp);
			fgets(buffer,80,fp);
			sscanf(buffer,"%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf", &tau0xf, &tau0xb, &tau1x, &thet0, &thet1);
			fgets(buffer,80,fp);
			fgets(buffer,80,fp);
			sscanf(buffer,"%lf\t\t%lf", &hselfx, &hlatex);
			fgets(buffer,80,fp);
		if(matcheck==1){
			assert(nsmx<=NSYSMX);
			nSYS[phase] = nsmx;

			if(thet0<thet1){
				PError("Initial hardening rate is lower than final hardening rate !!",230);
			}
			if(tau1x<1.E-6){	// linear hardening, independent of tau0
				tau1x =  1.E-6;	// avoid division by zero
				thet0 = thet1;
			}

			for(j=0;j<nsmx;j++){
				// define strain rate sensitivity and crss
				nSRS[phase][j] = nrsx;
				gam0[phase][j] = gamd0x;
				tau[phase][j][0] = tau0xf;
				tau[phase][j][1] = tau0xb;
				tau[phase][j][2] = tau1x;
				thet[phase][j][0] = thet0;
				thet[phase][j][1] = thet1;

				fgets(buffer,80,fp);
				sscanf(buffer,"%lf %lf %lf\t\t%lf %lf %lf",
						&bn[phase][j][0], &bn[phase][j][1], &bn[phase][j][2],
						&bb[phase][j][0], &bb[phase][j][1], &bb[phase][j][2]);
			}
		//	fgets(buffer,80,fp);
		}
		else{
			for(j=0;j<nsmx;j++)
				fgets(buffer,80,fp);
		}
	}

	fclose(fp);


	for(i=0;i<N_modes;i++){
	/* normalize bn and bb */
	if(mpirank==0){
		printf("Slip mode#%d:\n", iMode[i]);
		printf("Maximum slip systems = %d\n",nSYS[phase]);
	}
	for(j=0;j<nSYS[phase];j++){
		norm = sqrt(pow(bn[phase][j][0],2.0) + pow(bn[phase][j][1],2.0) + pow(bn[phase][j][2],2.0));
		for(mx=0;mx<3;mx++){
			bn[phase][j][mx] /= norm;
		}
		if(mpirank==0)printf("bn[%d][%d] = (%lf, %lf, %lf)\n",i,j,bn[phase][j][0],bn[phase][j][1],bn[phase][j][2]);

		norm = sqrt(pow(bb[phase][j][0],2.0) + pow(bb[phase][j][1],2.0) + pow(bb[phase][j][2],2.0));
#ifdef DD_BASED_FLAG
#ifdef DD_BCC
		bb_B[phase][j] = a0[phase]*sqrt(3.0)/2.0;	// nm
#else
		bb_B[phase][j] = a0[phase]/sqrt(2.0);	// nm
#endif
#endif
		for(mx=0;mx<3;mx++){
			bb[phase][j][mx] /= norm;
		}
		if(mpirank==0)printf("bb[%d][%d] = (%lf, %lf, %lf)\n",i,j,bb[phase][j][0],bb[phase][j][1],bb[phase][j][2]);

		norm = 0.0;	// used to check orthogonality
		for(mx=0;mx<3;mx++){
			norm += bn[phase][j][mx]*bb[phase][j][mx];
		}
		if(norm>1E-3){
			PError("Slip systems is NOT orthogonal !!",7);
		}
#ifdef DD_BASED_FLAG
		// sense vector: p = b X n
		bp[phase][j][0] = bb[phase][j][1]*bn[phase][j][2]-bb[phase][j][2]*bn[phase][j][1];
		bp[phase][j][1] = bb[phase][j][2]*bn[phase][j][0]-bb[phase][j][0]*bn[phase][j][2];
		bp[phase][j][2] = bb[phase][j][0]*bn[phase][j][1]-bb[phase][j][1]*bn[phase][j][0];
#endif

		// Schmid tensors
		T2_loop{
			aux33[mi][mj] = (bn[phase][j][mi]*bb[phase][j][mj]+bn[phase][j][mj]*bb[phase][j][mi])/2.0;
		}
		chg_basis5(aux5,aux33,aux66,aux3333,2);
		for(mx=0;mx<5;mx++){
			Schm_xt[phase][j][mx] = aux5[mx];
		}
	}
	}

	// initialize self and latent hardening coeff. for each system of the phase
	// notice here we assue latent hardeing is the same for all the off diangonal
	// elements, which is of course not always the case.
	for(i=0;i<nSYS[phase];i++){
		for(j=0;j<nSYS[phase];j++){
			if(i==j){
				Hard[phase][i][j] = hselfx;
			}
			else{
				Hard[phase][i][j] = hlatex;
			}
		}
	}

	return;
}/*end PlasticInit()*/

static void InitMicroStruct(char *s)
{
	ten4th caux3333;
	voigt aux6;
	ten2nd aux33;
	//voigt66 caux66;
	voigt66 c066_local = {0.0};
	real ph, th, om;
	int jgr, jph;
	int ii,jj,kk;
	int nph1, nph1_all;
	FILE *fp;
	char buffer[80] = {0};
	int i, idx;
	int EmptySteps = mpirank*Nxyz/NumPE;
	ten2nd sa2xt;	// transform matrix (sample -> xtal)
	std::vector<G_Info>::iterator it, end;
	std::vector<G_Info> local_gID_list;
	int tmp_flag;

	fp = fopen(s,"r");
	/*empty reading to go to the corresponding slabbed region*/
	for(i=0;i<EmptySteps;i++){
		fgets(buffer,80,fp);
	}

	nph1 = 0;
	local_loop{
		fgets(buffer,80,fp);
		sscanf(buffer,"%lf %lf %lf %d %d %d %d %d", &ph, &th, &om, &ii, &jj, &kk, &jgr, &jph);
		idx = ((ii-1-lxs)*CellDim[1]+jj-1)*CellDim[2]+kk-1; 
		if(idx!=pIDX){
			PError("Error in reading microstructure data file (inconsistent index system)!!", 1139);
		}

		if(jph==1) nph1++;
		grain_f[pIDX] = jgr; 
		phase_f[pIDX] = jph; 

		if(!Type_phases[jph-1]){	// NOT gas!!
			ph *= PI/180.;	th *= PI/180.;	om *= PI/180.;
			EulerToTransMatrix(&ph, &th, &om, sa2xt, 2);
			T2_loop{
				TranMat_xt2sa[pIDX][mi][mj] = sa2xt[mj][mi];
			}
			Ten4thTransform(caux3333,sa2xt,Cijkl[jph-1],2);	// Cijkl[jph] is in xtal ref. Do inverse transform
			chg_basis(aux6, aux33, caux66,caux3333,4);
			C6_loop{
				C_gr[pIDX][mi][mj] = caux66[mi][mj];
				c066_local[mi][mj] += caux66[mi][mj];
			}
			
			/* construct local_gID_list */
			end=local_gID_list.end();
			tmp_flag=0;
			for(it=local_gID_list.begin();it!=end;it++){
				if(grain_f[pIDX]==(*it).ID){
					tmp_flag=1;
					break;
				}
			}
			if(tmp_flag==0){
				struct G_Info tmp_ID={grain_f[pIDX],ph,th,om};
				local_gID_list.push_back(tmp_ID);
			}

		}
	}
	fclose(fp);

	/* construct global gID_list */
	int local_size = local_gID_list.size();
	int global_size;
	MPI_Allreduce(&local_size, &global_size, 1,
			MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	std::vector<G_Info> global_gID_list(global_size);
	std::vector<int> nsize(NumPE);
	MPI_Allgather(&local_size,1,MPI_INT,&nsize[0],1,MPI_INT,MPI_COMM_WORLD);
	std::vector<int> disps(NumPE);
	for(int i=0;i<NumPE;i++){
		disps[i]=(i>0)?(disps[i-1]+nsize[i-1]):0;
	}

	// define G_ID_Entry
	struct G_Info value = {0};
	MPI_Datatype G_ID_Entry;
	int count=4;
	int block_lens[4]={1,1,1,1};
	MPI_Aint indices[4];
	MPI_Datatype old_types[4]={MPI_INT, MPI_real, MPI_real, MPI_real};
	MPI_Address(&value,&indices[0]);
	MPI_Address(&value.t1,&indices[1]);
	MPI_Address(&value.Phi,&indices[2]);
	MPI_Address(&value.t2,&indices[3]);
	// make relative
	for(int i=count-1;i>=0;i--)
		indices[i] -= indices[0];
	MPI_Type_struct(count,block_lens,indices,old_types,&G_ID_Entry);
	MPI_Type_commit(&G_ID_Entry);

	MPI_Allgatherv(&local_gID_list[0],local_size,G_ID_Entry,
			&global_gID_list[0],&nsize[0],&disps[0],G_ID_Entry,MPI_COMM_WORLD);
	std::vector<G_Info>::iterator it_g, end_g;
	end_g=global_gID_list.end();
	for(it_g=global_gID_list.begin();it_g!=end_g;it_g++){
		end=gID_list.end();
		int check=0;
		for(it=gID_list.begin();it!=end;it++){
			if((*it).ID==(*it_g).ID){
				check=1;
				break;
			}
		}
		if(check==0){
			gID_list.push_back((*it_g));
		}
	}
	/* sort gID_list */
	sort(gID_list.begin(),gID_list.end(),G_Info::before);
	init_gid=gID_list.back().ID;
	if(mpirank==0){
		printf("The initial microstructure contains %d grains.\n",gID_list.size());
		printf("The initial largest grain ID is %d\n",gID_list.back().ID);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&nph1, &nph1_all, 1, MPI_INT,
			MPI_SUM, MPI_COMM_WORLD);
	Wgt_ph1 = (1.0*(real)nph1_all)/Nxyz;

	C6_loop{
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&c066_local[mi][mj], &C066[mi][mj], 1, MPI_real,
				MPI_SUM, MPI_COMM_WORLD);
		C066[mi][mj] *= WGT;
	}
	printf("Value of C066[1][1] is %le\n",C066[0][0]);

	C6_loop{
		S066[mi][mj] = C066[mi][mj];
	}

	LU_inv_66(S066);

	chg_basis(aux6,aux33,C066,C0,3);
	chg_basis(aux6,aux33,S066,S0,3);

//	/* test LU_inv_66 */
//	voigt66 s0 = {0.0};
//	s0[0][1] = 1.0;
//	s0[1][0] = 2.0; s0[1][2] = 2.0;
//	s0[2][1] = 3.0; s0[2][3] = 1.0;
//	s0[3][2] = 1.0; s0[3][4] = 2.0;
//	s0[4][3] = 3.0; s0[4][5] = 1.0;
//	s0[5][4] = 2.0;
//	printf("pyz: Before inverse:\n");
//	PrintVoigt66(s0);
//	LU_inv_66(s0);
//	printf("pyz: After inverse:\n");
//	PrintVoigt66(s0);


	return;
}/*end InitMicroStruct()*/

#ifdef DRX_ELASTIC_TREAT
static void InitFSField()
{
	ten4th caux3333;
	voigt aux6;
	ten2nd aux33;
	voigt66 saux66;
	voigt66 taux66;
	real dum;
	int itmp;
	int jph;

	local_loop{
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			C6_loop{
				saux66[mi][mj] = C_gr[pIDX][mi][mj];
			}
			LU_inv_66(saux66);
		}
		else{
			C6_loop{
				saux66[mi][mj] = 0.0;
			}
		}
		C6_loop{
			dum = 0.0;
			for(itmp=0;itmp<6;itmp++){
				dum += C066[mi][itmp]*saux66[itmp][mj];
			}
			taux66[mi][mj] = (mi==mj)+dum;
		}
		LU_inv_66(taux66);
		C6_loop{
			dum = 0.0;
			for(itmp=0;itmp<6;itmp++){
				dum += saux66[mi][itmp]*taux66[itmp][mj];
			}
			FS_gr[pIDX][mi][mj] = dum;
		}
	}

	return;
}/* InitFSField()*/
#endif

static void InitFields(void)
{
	int jph;
	voigt aux6;
	ten2nd aux33;
	voigt66 cgaux66;
	ten4th cg;
  char fname[80]={0};

	// Macro strain
	T2_loop{
		DisGradAvg[mi][mj] = Udot[mi][mj]*TimeStep;
		dDisGradAvg[mi][mj] = 0.0;
		dDisGradAvg_acum[mi][mj] = 0.0;
	}

	local_loop{
		//jgr = grain_f[pIDX];	// there could be thermo strain associated with grains. Here we ignore it for the time being.
		// strain, strain rate, disp. grad.
		T2_loop{
			Edot[pIDX][mi][mj] = 0.0;
			Eps[pIDX][mi][mj] = 0.0;
			fluctuation[pIDX][mi][mj] = 0.0;
			DisGrad[pIDX][mi][mj] = 0.0;
		}

		// inital guess for the stress, assume elasticity
		jph = phase_f[pIDX];
		if(!Type_phases[jph-1]){
			C6_loop{
				cgaux66[mi][mj] = C_gr[pIDX][mi][mj];
			}
			chg_basis(aux6,aux33,cgaux66,cg,3);
			T2_loop{
				Sig[pIDX][mi][mj] = 0.0;
				T2p_loop{
					Sig[pIDX][mi][mj] += cg[mi][mj][mip][mjp]*DisGradAvg[mip][mjp];
				}
				if(CREEP_FLAG==1){
					Sig[pIDX][mi][mj] += Scauchy[mi][mj];
				}
			}
			//PrintTensorNorm(Sig[pIDX]);
		}
		else{	// gas
			T2_loop{
				Sig[pIDX][mi][mj] = 0.0;
			}
		}
	}

	e_vm = d_vm*TimeStep;

	if(mpirank==0){
		fprintf(fp_vm,"TIME\tEVM\tEVMP\tDVM\tDVMP\tSVM\tSVMSVM_1\n");
#ifdef DD_BASED_FLAG
		fprintf(fp_dd,"TIME\tEVM\tDD_F\tDD_P\tDD_s\tDD_m\tDD_g\n");
#ifdef PF_DRX
		fprintf(fp_drx,"TIME\tEVM\tNEW_IDs\n");
#endif
#endif
		fprintf(fp_err,"IT\tERRE\tERRS\tSVM\n");
	}

  /* 09/16/15 --- PYZ */
  /* statistic of DRX nucleation types */
#ifdef PF_DRX
  //sprintf(fname,"drx_nucl_stat_PE%02d.out",mpirank);
  //fp_stat = fopen(fname,"w+");
#endif

	return;
}/*end InitFields()*/

static void M3Inverse(real a[3][3], real ia[3][3])
{
	int i, j;
	real deta = 0.0;

	deta = a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]) \
		   -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]) \
		   +a[0][2]*(a[1][0]*a[2][1]-a[2][0]*a[1][1]);
	ia[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
	ia[0][1] = a[0][2]*a[2][1] - a[0][1]*a[2][2];
	ia[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
	ia[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
	ia[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
	ia[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
	ia[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
	ia[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1];
	ia[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			ia[i][j] /=deta;

	return;
}/*end M3Inverse()*/

static void CalcGAMMA(void)
{
	int i,j,k,l;
	real nvector[3];
	real iomega[3][3], omega[3][3];

	local_loop{
		C4_loop{
			GAMMA[pIDX][mi][mj][mk][ml] = 0.0;
		}
	}

	local_loop{
		nvector[0] = g[pIDX].x;
		nvector[1] = g[pIDX].y;
		nvector[2] = g[pIDX].z;
		if(sIDX==0){		//k=0, boundary condition
			T2_loop{
				T2p_loop{
					if(ElastBC==0){
						GAMMA[pIDX][mi][mj][mip][mjp] = -1.0*S0[mi][mj][mip][mjp];
					}
					else{
						GAMMA[pIDX][mi][mj][mip][mjp] = 0.0;
					}
				}
			}
		}
		else{
			// Calculate Omega
			for(i=0;i<3;i++){
				for(k=0;k<3;k++){
					iomega[i][k] = 0.0;
					for(j=0;j<3;j++)
						for(l=0;l<3;l++)
							iomega[i][k] += C0[i][j][k][l]*nvector[j]*nvector[l];
				}
			}

			M3Inverse(iomega, omega);

			T2_loop{
				T2p_loop{
					GAMMA[pIDX][mi][mj][mip][mjp] = -1.0*omega[mi][mip]*nvector[mj]*nvector[mjp];
				}
			}
		}
	}

	return;
}/*end CalcGAMMA()*/

static void InteractionStrengthMatrix(char *type, real self, real coplane,
		real cross, real glissile, real hirth, real lomer, real **chi)
{
	/* The assignment of interaction strength matrix dependes on the label
	   of the slip system. Currently only FCC is considered according to
	  "Arsenlis and Parks, JMPS, 50(2002), 1979" */
	int i,j;
	if(strcmp(type, "FCC")==0||strcmp(type, "fcc")==0){
		chi[0][0] = self;
		chi[0][1] = coplane;
		chi[0][2] = coplane;
		chi[0][3] = hirth;
		chi[0][4] = glissile;
		chi[0][5] = lomer;
		chi[0][6] = cross;
		chi[0][7] = glissile;
		chi[0][8] = glissile;
		chi[0][9] = hirth;
		chi[0][10] = lomer;
		chi[0][11] = glissile;

		chi[1][1] = self;
		chi[1][2] = coplane;
		chi[1][3] = glissile;
		chi[1][4] = cross;
		chi[1][5] = glissile;
		chi[1][6] = glissile;
		chi[1][7] = hirth;
		chi[1][8] = lomer;
		chi[1][9] = lomer;
		chi[1][10] = hirth;
		chi[1][11] = glissile;

		chi[2][2] = self;
		chi[2][3] = lomer;
		chi[2][4] = glissile;
		chi[2][5] = hirth;
		chi[2][6] = glissile;
		chi[2][7] = lomer;
		chi[2][8] = hirth;
		chi[2][9] = glissile;
		chi[2][10] = glissile;
		chi[2][11] = cross;

		chi[3][3] = self;
		chi[3][4] = coplane;
		chi[3][5] = coplane;
		chi[3][6] = hirth;
		chi[3][7] = lomer;
		chi[3][8] = glissile;
		chi[3][9] = cross;
		chi[3][10] = glissile;
		chi[3][11] = glissile;

		chi[4][4] = self;
		chi[4][5] = coplane;
		chi[4][6] = lomer;
		chi[4][7] = hirth;
		chi[4][8] = glissile;
		chi[4][9] = glissile;
		chi[4][10] = hirth;
		chi[4][11] = lomer;

		chi[5][5] = self;
		chi[5][6] = glissile;
		chi[5][7] = glissile;
		chi[5][8] = cross;
		chi[5][9] = glissile;
		chi[5][10] = lomer;
		chi[5][11] = hirth;

		chi[6][6] = self;
		chi[6][7] = coplane;
		chi[6][8] = coplane;
		chi[6][9] = hirth;
		chi[6][10] = glissile;
		chi[6][11] = lomer;

		chi[7][7] = self;
		chi[7][8] = coplane;
		chi[7][9] = glissile;
		chi[7][10] = cross;
		chi[7][11] = glissile;

		chi[8][8] = self;
		chi[8][9] = lomer;
		chi[8][10] = glissile;
		chi[8][11] = hirth;

		chi[9][9] = self;
		chi[9][10] = coplane;
		chi[9][11] = coplane;

		chi[10][10] = self;
		chi[10][11] = coplane;

		chi[11][11] = self;

		/* symmetric matrix */
		for(i=1;i<12;i++){
			for(j=0;j<i;j++){
				chi[i][j] = chi[j][i];
			}
		}
	}
	else{
		PError("Interaction strength matrix is only supported for FCC currently!!", 1110);	
	}
	return;
}/*end InteractionStrengthMatrix()*/

#ifdef DD_BASED_FLAG
#ifdef DD_AL
static void T_DependentElastConst(real *elastConst, real T)
{
	/* Balasubramanian and Anand, JMPS, 2002;50;101 */
	/* The anisotropic elastic tensor of aluminium */
	real k11, k12, k44;

	// GPa
	k11 = 123.323 + (6.7008E-8)*pow(T,3.0) - (1.1342E-4)*pow(T,2.0) - (7.8788E-3)*T;
	k12 = 70.6512 + (4.4105E-8)*pow(T,3.0) - (7.5498E-5)*pow(T,2.0) - (3.9992E-3)*T;
	k44 = 31.2071 + (7.0477E-9)*pow(T,3.0) - (1.2136E-5)*pow(T,2.0) - (8.3274E-3)*T;
	
	// MPa
	elastConst[0] = k11*1000;
	elastConst[1] = k12*1000;
	elastConst[2] = k44*1000;

	return;
}/*end T_DependentElastConst()*/
#endif
#endif

void SetupJob(void)
{
	int i,j,m;
	int jph;
	real sense_vector[3];
	voigt aux6;
	ten2nd aux33;
	voigt66 aux66;
	ten4th aux3333;
	ten2nd dev_Udot_s;

	if(mpirank==0){
		printf("\n========================================\n");
		printf("------Start setting up the job......------\n");
	}

	// total number of grid points
	Nxyz = CellDim[0]*CellDim[1]*CellDim[2];
	WGT = 1./Nxyz;

	assert(N_phases<=NPHMAX);

	/************************
	  fft_mpi initialization 
	 ************************/
	plan = fftw3d_mpi_create_plan(MPI_COMM_WORLD, CellDim[0], CellDim[1], CellDim[2],
				FFTW_FORWARD, FFTW_ESTIMATE);
	plan_pf = fftw3d_mpi_create_plan(MPI_COMM_WORLD, CellDim_pf[0], CellDim_pf[1], CellDim_pf[2],
				FFTW_FORWARD, FFTW_ESTIMATE);
	iplan = fftw3d_mpi_create_plan(MPI_COMM_WORLD, CellDim[0], CellDim[1], CellDim[2],
				FFTW_BACKWARD, FFTW_ESTIMATE);
	/* slab decomposition 
	   all PEs work coorporatively
	   to perform FFT*/
	fftwnd_mpi_local_sizes(plan, &lnx, &lxs, &lnyt, &lyst, &lsize);
	fftwnd_mpi_local_sizes(plan_pf, &lnx_pf, &lxs_pf, &lnyt_pf, &lyst_pf, &lsize_pf);
//	exit(0);
	fft_data = (fftw_complex*)fftw_malloc(lsize*sizeof(fftw_complex));
	fft_work = (fftw_complex*)fftw_malloc(lsize*sizeof(fftw_complex));

	AllocMem(kSig_r,lsize,ten2nd);	// stress in k-space, real part
	AllocMem(kSig_i,lsize,ten2nd);	// stress in k-space, imaginary part

	if(mpirank==0) printf("\nfft_mpi initialized\n");

	/************************
	 k-space setup
	 ************************/
	AllocMem(g, lsize, Vec3R);
	AllocMem(g2, lsize, real);
	AllocMem(GAMMA, lsize, ten4th);
	real tmp[3];
	local_loop{
		tmp[0] = (real)(px+lxs) * (real)(1.0/CellDim[0]);
		tmp[1] = (real)(py) * (real)(1.0/CellDim[1]);
		tmp[2] = (real)(pz) * (real)(1.0/CellDim[2]);
		for(m=0; m<3; m++){
			tmp[m] -= round(tmp[m]);
			tmp[m] *= TWO_PI;
		}

		VSet(g[pIDX], tmp[0], tmp[1], tmp[2]);
		g2[pIDX] = VSecNorm(g[pIDX]);

		if((px+lxs+py+pz) == 0){
			g2[pIDX] = 1.0;
			VZero(g[pIDX]);
		}
		VScale(g[pIDX], 1./g2[pIDX]);
	}
	if(mpirank==0) printf("\nk-space initialized\n");

	/*******************************
	  Material properties initialization
	  ******************************/
	chg_basis(aux6, aux33, aux66, aux3333, 0);	// calculate B[6]
	if(mpirank==0){
	//	PrintB();
	}
	AllocMem(nSRS,N_phases,real*);
	AllocMem(nSYS,N_phases,int);
	for(i=0;i<N_phases;i++){
		nSYS[i] = 0;
		AllocMem(nSRS[i],NSYSMX,real);
	}
	AllocMem(C_gr,lsize,voigt66);
	AllocMem(TranMat_xt2sa,lsize,ten2nd);
	AllocMem(Schm_gr,lsize,voigtch);
	if(!Type_phases[0]){
#ifdef DD_BASED_FLAG
#ifdef DD_AL
		T_DependentElastConst(ElastConst_PhaseI, T__K);
#endif

#ifdef DD_CU
		ElastConst_PhaseI[0] = (1.818-0.0004045*T__K)*100*1000;	// MPa
		ElastConst_PhaseI[1] = (1.285-0.0002*T__K)*100*1000;	// MPa
		ElastConst_PhaseI[2] = (0.8396-0.0002711*T__K)*100*1000;	// MPa
#endif
		Shear_G[0] = 0.1*(ElastConst_PhaseI[0] - ElastConst_PhaseI[1] + 3*ElastConst_PhaseI[2]) + (0.5)*(5*(ElastConst_PhaseI[0] - ElastConst_PhaseI[1])*ElastConst_PhaseI[2])/(4*ElastConst_PhaseI[2] + 3*(ElastConst_PhaseI[0] - ElastConst_PhaseI[1]));
#endif
		if(mpirank==0){
			printf("Phase I: elastic stiffness(c11,c12,c44)=(%lf,%lf,%lf)\n",
					ElastConst_PhaseI[0],ElastConst_PhaseI[1],ElastConst_PhaseI[2]);
		}
		ElastStiffnessMatrix(ElastConst_PhaseI, Cijkl[0]);
		PlasticInit(Slip_PhaseI, 0);
	}
	if(N_phases==2&&(!Type_phases[1])){
#ifdef DD_BASED_FLAG
#ifdef DD_AL
		T_DependentElastConst(ElastConst_PhaseII, T__K);
#endif
#ifdef DD_CU
		ElastConst_PhaseI[0] = (1.818-0.0004045*T__K)*100*1000;	// MPa
		ElastConst_PhaseI[1] = (1.285-0.0002*T__K)*100*1000;	// MPa
		ElastConst_PhaseI[2] = (0.8396-0.0002711*T__K)*100*1000;	// MPa
#endif
		Shear_G[1] = 0.1*(ElastConst_PhaseII[0] - ElastConst_PhaseII[1] + 3*ElastConst_PhaseII[2]) + (0.5)*(5*(ElastConst_PhaseII[0] - ElastConst_PhaseII[1])*ElastConst_PhaseII[2])/(4*ElastConst_PhaseII[2] + 3*(ElastConst_PhaseII[0] - ElastConst_PhaseII[1]));
#endif
		if(mpirank==0){
			printf("Phase II: elastic stiffness(c11,c12,c44)=(%lf,%lf,%lf)\n",
					ElastConst_PhaseII[0],ElastConst_PhaseII[1],ElastConst_PhaseII[2]);
		}
		ElastStiffnessMatrix(ElastConst_PhaseII, Cijkl[1]);
		PlasticInit(Slip_PhaseII, 1);
	}
	AllocMem(grain_f,lsize,int);
        AllocMem(grain_f_new,lsize,int);
	AllocMem(phase_f,lsize,int);
	InitMicroStruct(initial_ms);
//	PrintSchmidXt();

	/************************
	  Boundary conditions
	 ************************/
	// strain rate
	if((VelGrad_BC_Flag[0]+VelGrad_BC_Flag[1]+VelGrad_BC_Flag[2])==2){
		PError("Cannot enforce only two deviatoric components (Check input VelGrad_BC_Flag)",116);
	}
	
	for(i=0;i<6;i++){
		VelGrad_BC[i] *= VelGrad_BC_Flag[i];
	}
	VoigtToFull(VelGrad_BC, Udot);
	if(mpirank==0){
		printf("Applied velocity gradient:\n");
	//	PrintTensor(Udot);
	}
	SymAntDecompose(Udot, Udot_s, Udot_a);
	chg_basis(D_bar6,Udot_s,aux66,aux3333,2);
	for(i=0;i<5;i++) D_bar5[i] = D_bar6[i];
	chg_basis5(D_bar5,dev_Udot_s,aux66,aux3333,1);
	d_vm = 0.0;
	T2_loop{
		d_vm += pow(dev_Udot_s[mi][mj],2.0);
	}
	d_vm = sqrt(2./3.*d_vm);

	// Cauchy stress
	for(i=0;i<6;i++){	// check BC validity
		if((Stress_BC_Flag[i]*VelGrad_BC_Flag[i]!=0)||((Stress_BC_Flag[i]+VelGrad_BC_Flag[i])!=1))
			PError("Check boundary conditions on strain rate and stress!!",114);
	}
	CREEP_FLAG = 0;
	for(i=0;i<6;i++){	// creep?
		CREEP_FLAG += Stress_BC_Flag[i];
	}
	CREEP_FLAG /= 6;
	VoigtToFull(Stress_BC, Scauchy);
	if(mpirank==0){
		printf("Applied (Cauchy) stress:\n");
	//	PrintTensor(Scauchy);
	}

	/************************
	  constitutive law related
	 ************************/
#if  defined(DD_BASED_FLAG)
	AllocMem(rho_m, lsize,real*);
	AllocMem(rho_s, lsize,real*);
	AllocMem(rho_dipole, lsize,real*);
	AllocMem(rho_forest, lsize,real*);
	AllocMem(rho_inplane, lsize,real*);
	AllocMem(trial_rho_s, lsize,real*);
	AllocMem(rho_g1, lsize,real*);
	AllocMem(rho_g2, lsize,real*);
	AllocMem(rho_g3, lsize,real*);
	AllocMem(rho_dot_s, lsize,real*);
	AllocMem(rho_dot_g1, lsize,real*);
	AllocMem(rho_dot_g2, lsize,real*);
	AllocMem(gadot_grad, lsize,Vec3R*);
	AllocMem(ga_grad, lsize,Vec3R*);
	AllocMem(local_pxp_send, CellDim[1]*CellDim[2], real);
	AllocMem(local_pxm_send, CellDim[1]*CellDim[2], real);
	AllocMem(local_pxp_recv, CellDim[1]*CellDim[2], real);
	AllocMem(local_pxm_recv, CellDim[1]*CellDim[2], real);
	AllocMem(local_pxpG_send, CellDim[1]*CellDim[2], int);
	AllocMem(local_pxmG_send, CellDim[1]*CellDim[2], int);
	AllocMem(local_pxpG_recv, CellDim[1]*CellDim[2], int);
	AllocMem(local_pxmG_recv, CellDim[1]*CellDim[2], int);
	AllocMem(GB_checkX, lsize, int*);
	AllocMem(GB_checkY, lsize, int*);
	AllocMem(GB_checkZ, lsize, int*);
	AllocMem(gam_acum, lsize,real*);
	AllocMem(trial_gam_acum, lsize,real*);
	AllocMem(trial_rho_g1, lsize,real*);
	AllocMem(trial_rho_g2, lsize,real*);
	AllocMem(trial_rho_g3, lsize,real*);
	AllocMem(rho_F, lsize,real*);
	AllocMem(rho_P, lsize,real*);
	local_loop{
		AllocMem(rho_m[pIDX],NSYSMX,real);
		AllocMem(rho_s[pIDX],NSYSMX,real);
		AllocMem(rho_dipole[pIDX],NSYSMX,real);
		AllocMem(rho_forest[pIDX],NSYSMX,real);
		AllocMem(rho_inplane[pIDX],NSYSMX,real);
		AllocMem(trial_rho_s[pIDX],NSYSMX,real);
		AllocMem(rho_g1[pIDX],NSYSMX,real);
		AllocMem(rho_g2[pIDX],NSYSMX,real);
		AllocMem(rho_g3[pIDX],NSYSMX,real);
	  AllocMem(rho_dot_s[pIDX], NSYSMX,real);
	  AllocMem(rho_dot_g1[pIDX], NSYSMX,real);
	  AllocMem(rho_dot_g2[pIDX], NSYSMX,real);
		AllocMem(gadot_grad[pIDX],NSYSMX,Vec3R);
		AllocMem(ga_grad[pIDX],NSYSMX,Vec3R);
		AllocMem(gam_acum[pIDX], NSYSMX,real);
		AllocMem(trial_gam_acum[pIDX], NSYSMX,real);
		AllocMem(trial_rho_g1[pIDX],NSYSMX,real);
		AllocMem(trial_rho_g2[pIDX],NSYSMX,real);
		AllocMem(trial_rho_g3[pIDX],NSYSMX,real);
		AllocMem(rho_F[pIDX],NSYSMX,real);
		AllocMem(rho_P[pIDX],NSYSMX,real);
		AllocMem(GB_checkX[pIDX],2,int);
		AllocMem(GB_checkY[pIDX],2,int);
		AllocMem(GB_checkZ[pIDX],2,int);
	}
	AllocMem(rho_tot,lsize,real);
	//AllocMem(dd_average,lsize,real);
#ifdef PF_DRX
	assert(CellDim_pf[0]%NumPE==0);
	if(CellDim[0]!=CellDim_pf[0]||CellDim[1]!=CellDim_pf[1]||CellDim[2]!=CellDim_pf[2]){
    DRX_Interpolation_Flag=1;
  }else{
    DRX_Interpolation_Flag=0;
  }
	lx_pf = CellDim_pf[0]/NumPE; ly_pf = CellDim_pf[1]; lz_pf = CellDim_pf[2];
	lsize_pf = lx_pf*ly_pf*lz_pf;
  L00 = L0;
  TimeStep_DRX = pf_time_step;  // [s] time step in phase-field simulation
	if(DRX_Interpolation_Flag==1){
		if(mpirank==0){
			printf("FFT and PF have different grid size:\nAllocate grain_f_drx, rho_tot_drx, gID_new_drx, and grex_new_drx\n");
		}
  	AllocMem(rho_tot_drx,lsize_pf,real);
  //	exit(0);
  	
  	AllocMem(dd_average_drx,lsize,real);
  	AllocMem(grain_f_drx,lsize_pf,int);
  //	AllocMem(gID_rex_drx,lsize_pf,int);
  	AllocMem(gID_new_drx,lsize_pf,int);
  //	AllocMem(grex_new_drx,lsize_pf,int);
  	AllocMem(first_pf_drx,lsize_pf,int);
  }
  //it(0);
  /* The following are on FFT grid */
	AllocMem(diff_rho,lsize,real);
	AllocMem(gID_rex,lsize,int); 
	AllocMem(gID_new,lsize,int);
	AllocMem(nuclei,lsize,int);
	AllocMem(nucleation_count,lsize,int);
        AllocMem(growth,lsize,int);
	AllocMem(swap_check,lsize,int);
//	AllocMem(grex_new,lsize,int);
	AllocMem(first_pf,lsize,int);
	AllocMem(GB_indicator,lsize,int);
	AllocMem(GB_type,lsize,int);
	AllocMem(t_last,lsize,real);
	if(Nucl_Static_Flag==1){
		AllocMem(kappa_drx,lsize,real);	// static nucleation strength pre-determined
		AllocMem(k_c_field,lsize,real);	// the critical nucleation strength varies with the local GB migration activation energy
	}
//exit(0);
 srand(time(NULL));
	gsl_rng_default_seed = rand();
	RandType = gsl_rng_default;
	RandInstance = gsl_rng_alloc(RandType);
//exit(0);
  local_loop{
    diff_rho[pIDX] = 0.0;
    gID_rex[pIDX] = 0;
    nuclei[pIDX] = 0;
    nucleation_count[pIDX] = 0;
    growth[pIDX] = 0;
     swap_check[pIDX] = 0;
    gID_new[pIDX] = 0;
 //   grex_new[pIDX] = 0;
    first_pf[pIDX] = 0;
    GB_indicator[pIDX] = 0;
    GB_type[pIDX] = 0;
    t_last[pIDX] = -1E3*tau_DRX*TimeStep;
  }

	if(DRX_Interpolation_Flag==1){
		for(int i=0;i<lx_pf;i++){
			for(int j=0;j<ly_pf;j++){
				for(int k=0;k<lz_pf;k++){
					int idx = (i*ly_pf+j)*lz_pf+k;
				//exit(0);
          rho_tot_drx[idx]=0.0;
        //  dd_average_drx[idx]=0.0;
         // first_pf_drx[idx]=0;
          grain_f_drx[idx]=0;
      //    gID_rex_drx[idx]=0;
          gID_new_drx[idx]=0;
        //exit(0);
     //     grex_new_drx[idx]=0;
        }
      }
      
    }
   
  }
  
	if(Nucl_Static_Flag==1){
    local_loop{
			kappa_drx[pIDX]=gsl_ran_weibull(RandInstance,k_c,alpha);
      k_c_field[pIDX]=k_c;  //initially all the same
    }
	}
	
	
#ifdef DRX_ELASTIC_TREAT
	AllocMem(Eps0_last,lsize,ten2nd);
	AllocMem(Tau_0,lsize, ten2nd);
	AllocMem(Tau,lsize, ten2nd);
	AllocMem(FS_gr,lsize,voigt66);
	InitFSField();
#endif
#endif
	local_loop{
		for(i=0;i<NSYSMX;i++){
			rho_m[pIDX][i] = 0.0;
		//	jph1 = phase_f[pIDX];
		
			rho_s[pIDX][i] = rho_SSD_initial;	// 1E12/m^2
			trial_rho_s[pIDX][i] = rho_SSD_initial;	// 1E12/m^2
		
	
			rho_g1[pIDX][i] = 0.0;	// GND is initially zero
			rho_g2[pIDX][i] = 0.0;	// GND is initially zero
			rho_g3[pIDX][i] = 0.0;	// GND is initially zero
      rho_dot_s[pIDX][i] = 0.0;
      rho_dot_g1[pIDX][i] = 0.0;
      rho_dot_g2[pIDX][i] = 0.0;
			VZero(gadot_grad[pIDX][i]);
			VZero(ga_grad[pIDX][i]);
			gam_acum[pIDX][i] = 0.0;	// accumulated shear on each slip system
			trial_gam_acum[pIDX][i] = 0.0;	// accumulated shear on each slip system, trial version
			trial_rho_g1[pIDX][i] = 0.0;	// GND is initially zero
			trial_rho_g2[pIDX][i] = 0.0;	// GND is initially zero
			trial_rho_g3[pIDX][i] = 0.0;	// GND is initially zero
			rho_F[pIDX][i] = 0.0;
			rho_P[pIDX][i] = 0.0;
			rho_dipole[pIDX][i] = 0.0;
			rho_forest[pIDX][i] = 0.0;
			rho_inplane[pIDX][i] = 0.0;
		}
		rho_tot[pIDX] = 0.0;
		//dd_average[pIDX]=0.0;
	}
	// interaction strength matrix
	for(jph=0;jph<NPHMAX;jph++){
		AllocMem(Chi[jph],NSYSMX,real*);
		for(i=0;i<NSYSMX;i++){
			AllocMem(Chi[jph][i],NSYSMX,real);
		}
		AllocMem(Chi_1[jph],NSYSMX,real*);
		for(i=0;i<NSYSMX;i++){
			AllocMem(Chi_1[jph][i],NSYSMX,real);
		}
		InteractionStrengthMatrix("FCC",Selfinter0[jph],Coplanar0[jph],
			0.0,0.0,0.0,0.0,
				Chi[jph]);
				InteractionStrengthMatrix("FCC",0.0,0.0,
				CrossSlip0[jph], GlissileJunction0[jph],
				HirthLock0[jph], LomerCottrellLock0[jph],
				Chi_1[jph]);
	}
	// sin and cos values for projection
	for(jph=0;jph<NPHMAX;jph++){
		for(i=0;i<nSYS[jph];i++){
			for(j=0;j<nSYS[jph];j++){
				sense_vector[0] = bn[jph][j][1]*bb[jph][j][2]-bn[jph][j][2]*bb[jph][j][1];
				sense_vector[1] = bn[jph][j][2]*bb[jph][j][0]-bn[jph][j][0]*bb[jph][j][2];
				sense_vector[2] = bn[jph][j][0]*bb[jph][j][1]-bn[jph][j][1]*bb[jph][j][0];

				cos_theta_s[jph][i][j] = fabs(bn[jph][i][0]*sense_vector[0]+
					bn[jph][i][1]*sense_vector[1]+bn[jph][i][2]*sense_vector[2]);
				sin_theta_s[jph][i][j] = fabs(sqrt(fabs(1.0-cos_theta_s[jph][i][j]*cos_theta_s[jph][i][j])));

				cos_theta_g1[jph][i][j] = fabs(bn[jph][i][0]*bb[jph][j][0]+
					bn[jph][i][1]*bb[jph][j][1]+bn[jph][i][2]*bb[jph][j][2]);
				sin_theta_g1[jph][i][j] = fabs(sqrt(fabs(1.0-cos_theta_g1[jph][i][j]*cos_theta_g1[jph][i][j])));

				cos_theta_g2[jph][i][j] = cos_theta_s[jph][i][j];
				sin_theta_g2[jph][i][j] = sin_theta_s[jph][i][j];

				cos_theta_g3[jph][i][j] = fabs(bn[jph][i][0]*bn[jph][j][0]+
					bn[jph][i][1]*bn[jph][j][1]+bn[jph][i][2]*bn[jph][j][2]);
				sin_theta_g3[jph][i][j] = fabs(sqrt(fabs(1.0-cos_theta_g3[jph][i][j]*cos_theta_g3[jph][i][j])));
			}
		}
	}
#ifdef DD_GND
	Gradient_ExchangeGrainID();
#endif

#else
	// Initialize CRSS and accum shear for grains
	AllocMem(gamacum,lsize,real);
	AllocMem(xkin,lsize,real*);
	AllocMem(trial_tau,lsize,real**);
	AllocMem(crss,lsize,real**);
	local_loop{
		AllocMem(xkin[pIDX],NSYSMX,real);
		AllocMem(trial_tau[pIDX],NSYSMX,real*);
		AllocMem(crss[pIDX],NSYSMX,real*);
		for(i=0;i<NSYSMX;i++){
			AllocMem(trial_tau[pIDX][i],2,real);
			AllocMem(crss[pIDX][i],2,real);
		}
	}
	local_loop{
		gamacum[pIDX] = 0.0;
		jph = phase_f[pIDX];
		if(Type_phases[jph-1]==0){
			for(i=0;i<nSYS[jph-1];i++){
				crss[pIDX][i][0] = tau[jph-1][i][0];
				crss[pIDX][i][1] = tau[jph-1][i][1];
				trial_tau[pIDX][i][0] = tau[jph-1][i][0];
				trial_tau[pIDX][i][1] = tau[jph-1][i][1];
				xkin[pIDX][i] = 0.0;
			}
		}
	}
#endif
				


	/*******************************
	  global arrays for stress/strain...
	  arrays allocation and initialization
	  ******************************/
	AllocMem(Sig, lsize, ten2nd);
	AllocMem(fluctuation,lsize,ten2nd);
	AllocMem(Eps, lsize, ten2nd);
	AllocMem(Edot, lsize, ten2nd);
	AllocMem(DisGrad, lsize, ten2nd);
	AllocMem(VelGrad, lsize, ten2nd);
	AllocMem(gamdot,lsize,real*);
	AllocMem(new_position,lsize,real*);
	AllocMem(displacement_fluct,lsize,real*);
	local_loop{
		AllocMem(gamdot[pIDX],NSYSMX,real);
		AllocMem(new_position[pIDX],3,real);
		AllocMem(displacement_fluct[pIDX],3,real);
	}


	/*******************************
	  Initialize the stress/strain,rates...
	  ******************************/
	InitFields();

	/*******************************
	  Green's operator in k-space
	  ******************************/
	CalcGAMMA();	// calculate Green operator in k-space




	if(mpirank==0){
		printf("\n-----------Job set up-----------------\n");
		printf("========================================\n");
	}
	return;
}/*end SetupJob()*/

void DestroyJob(void)
{
	int i;

	if(mpirank==0){
		printf("\n========================================\n");
		printf("--------Start destroying job......----------\n");
	}

	/*******************************
	  Material properties initialization
	  ******************************/
	for(i=0;i<N_phases;i++)
		free(nSRS[i]);
	free(nSRS);
	free(nSYS);
	free(C_gr);
	free(TranMat_xt2sa);
	free(grain_f);
        free(grain_f_new);
	free(phase_f);
	free(Schm_gr);

	free(iMode);

	/************************
	  crss, hardening,...
	 ************************/
#ifdef DD_BASED_FLAG
	int jph;
	local_loop{
		free(rho_m[pIDX]);
		free(rho_s[pIDX]);
		free(rho_dipole[pIDX]);
		free(rho_forest[pIDX]);
		free(rho_inplane[pIDX]);
		free(trial_rho_s[pIDX]);
		free(rho_g1[pIDX]);
		free(rho_g2[pIDX]);
		free(rho_g3[pIDX]);
		free(rho_dot_s[pIDX]);
		free(rho_dot_g1[pIDX]);
		free(rho_dot_g2[pIDX]);
		free(trial_rho_g1[pIDX]);
		free(trial_rho_g2[pIDX]);
		free(trial_rho_g3[pIDX]);
		free(rho_F[pIDX]);
		free(rho_P[pIDX]);
		free(gadot_grad[pIDX]);
		free(gam_acum[pIDX]);
		free(trial_gam_acum[pIDX]);
		free(ga_grad[pIDX]);
		free(GB_checkX[pIDX]);
		free(GB_checkY[pIDX]);
		free(GB_checkZ[pIDX]);
	}
	free(rho_m);
	free(rho_s);
	free(rho_dipole);
	free(rho_forest);
	free(rho_inplane);
	free(trial_rho_s);
	free(rho_g1);
	free(rho_g2);
	free(rho_g3);
	free(rho_dot_s);
	free(rho_dot_g1);
	free(rho_dot_g2);
	free(trial_rho_g1);
	free(trial_rho_g2);
	free(trial_rho_g3);
	free(rho_F);
	free(rho_P);
	free(rho_tot);
	//free(dd_average);
#ifdef PF_DRX
	free(diff_rho);
	if(DRX_Interpolation_Flag==1){
    free(rho_tot_drx);
    free(grain_f_drx);
//	  free(gID_rex_drx);
	  free(gID_new_drx);
	//  free(grex_new_drx);
	  free(first_pf_drx);
  }
	free(gID_rex);
	free(nuclei);
	free(nucleation_count);
	free(growth);
    free(swap_check);
	free(gID_new);
//	free(grex_new);
	free(first_pf);
	free(GB_indicator);
	free(GB_type);
	free(t_last);
	if(Nucl_Static_Flag==1){
		free(kappa_drx);
		free(k_c_field);
	}
	gsl_rng_free(RandInstance);
#ifdef DRX_ELASTIC_TREAT
	free(Eps0_last);
	free(Tau_0);
	free(Tau);
	free(FS_gr);
#endif
#endif
	free(gadot_grad);
	free(gam_acum);
	free(trial_gam_acum);
	free(ga_grad);
	free(local_pxp_send);
	free(local_pxm_send);
	free(local_pxp_recv);
	free(local_pxm_recv);
	free(local_pxpG_send);
	free(local_pxmG_send);
	free(local_pxpG_recv);
	free(local_pxmG_recv);
	free(GB_checkX);
	free(GB_checkY);
	free(GB_checkZ);
	for(jph=0;jph<NPHMAX;jph++){
		for(i=0;i<NSYSMX;i++){
			free(Chi[jph][i]);
		}
		free(Chi[jph]);
			for(i=0;i<NSYSMX;i++){
			free(Chi_1[jph][i]);
		}
		free(Chi_1[jph]);
	}
#else
	free(gamacum);
	local_loop{
		for(i=0;i<NSYSMX;i++){
			free(trial_tau[pIDX][i]);
			free(crss[pIDX][i]);
		}
		free(xkin[pIDX]);
		free(trial_tau[pIDX]);
		free(crss[pIDX]);
	}
	free(xkin);
	free(trial_tau);
	free(crss);
#endif

	/************************
	  fft_mpi finalization
	 ************************/
	fftwnd_mpi_destroy_plan(plan);
		fftwnd_mpi_destroy_plan(plan_pf);
	fftwnd_mpi_destroy_plan(iplan);
	fftw_free(fft_data);
	fftw_free(fft_work);
	free(kSig_r);
	free(kSig_i);
	if(mpirank==0) printf("\nfft_mpi finalized\n");

	/************************
	 free k-space
	 ************************/
	free(g);
	free(g2);
	free(GAMMA);
	if(mpirank==0) printf("\nk-space finalized\n");


	/************************
	 free global arrays (stress/strain)
	 ************************/
	free(Sig);
	free(fluctuation);
	free(Eps);
	free(Edot);
	free(DisGrad);
	free(VelGrad);
	local_loop{
		free(gamdot[pIDX]);
		free(new_position[pIDX]);
		free(displacement_fluct[pIDX]);
	}
	free(gamdot);
	free(new_position);
	free(displacement_fluct);

  /* 09/16/15 --- PYZ */
#ifdef PF_DRX
  //fclose(fp_stat);
#endif

	if(mpirank==0){
		printf("\n------------Job destroyed-----------\n");
		printf("========================================\n");
	}
	return;
}/*end DestroyJob()*/

