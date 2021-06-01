#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "evp.h"

//typedef enum{N_I, N_R} VType;
typedef enum{N_I, N_R, N_S} VType;

typedef struct{
	char *vName;
	void *vPtr;
	VType vType;
	int vLen;
	int vStatus;
}NameList;

#define NameI(x) {#x, &x, N_I, sizeof(x)/sizeof(int)}
#define NameR(x) {#x, &x, N_R, sizeof(x)/sizeof(real)}
#define NameS(x) {#x, &x, N_S, 1}
#define NP_I ((int *)(nameList[k].vPtr) + j)
#define NP_R ((real *)(nameList[k].vPtr) + j)
#define NP_S ((char *)(nameList[k].vPtr) + j)

NameList nameList[] = {
	NameI(CellDim),
	NameI(N_phases),
	NameI(Type_phases),
	NameR(TimeStep),
	NameI(N_steps),
	NameI(Update_Flag),
	NameI(Hard_Flag),
	NameI(Tex_Flag),
	NameI(IterMax),
	NameI(PrintControl),
	NameR(Err),
	NameI(VelGrad_BC_Flag),
	NameR(VelGrad_BC),
	NameI(Stress_BC_Flag),
	NameR(Stress_BC),
	NameI(ElastBC),
	NameR(ElastConst_PhaseI),
	NameR(ElastConst_PhaseII),
	NameS(Slip_PhaseI),
	NameS(Slip_PhaseII),
	NameS(initial_ms),
#ifdef DD_BASED_FLAG
	NameR(a0),
	NameR(L0),
	NameR(GND_ScaleRatio),
	NameR(T__K),
	NameR(Q_slip),
	NameR(Q_bulk),
	NameR(C_1),
	NameR(C_2),
	NameR(C_3),
	NameR(C_4),
	NameR(C_5),
	NameR(C_6),
	NameR(C_7),
	NameR(C_8),
	NameR(Selfinter0),
	NameR(Coplanar0),
	NameR(CrossSlip0),
	NameR(GlissileJunction0),
	NameR(HirthLock0),
	NameR(LomerCottrellLock0),
	NameR(rho_SSD_initial),
	NameR(D_0),
	NameI(nucleus_radius),
#ifdef PF_DRX
	NameI(CellDim_pf),
	NameI(RandSeed),
	NameR(E_gb),
	NameR(M_gb),
//	NameR(M_bar),
	NameR(Scale_Fdeform),
	NameR(k_c),
	NameI(Nucl_Static_Flag),
	NameR(alpha),
  NameR(zeta_kc),
  NameR(age_drx),
	NameR(ssdScaler_DRX),
	NameR(gndScaler_DRX),
	NameR(sigScaler_DRX),
	NameR(tau_DRX),
	NameR(DeltaQ_DRX),
	NameR(DeltaTau_DRX),
	NameR(pf_length_scale),
	NameR(pf_time_step),
	
#endif
#endif
};

int GetNameList(char **argv){
	int j, k, match, ok;
	char buff[200], *token;
	FILE *fp;
	strcpy(buff, argv[1]);
	if((fp = fopen(buff, "r")) == 0){
		printf("PE#%d: Input file not found!\n", mpirank);
		return(1);
	}
	for(k=0; k<sizeof(nameList)/sizeof(NameList); k++)
		nameList[k].vStatus = 0;
	ok = 1;
	while(1){
		fgets(buff, 200, fp);
		if(feof(fp))
			break;
		if(buff[0]=='#'||buff[0]=='\n')
			continue;
		token = strtok(buff, " \t\n");
		if(!token)
			break;
		match = 0;
		for(k = 0; k<sizeof(nameList)/sizeof(NameList); k++){
			if(strcmp(token, nameList[k].vName) == 0){
					match = 1;
					if(nameList[k].vStatus == 0){
						nameList[k].vStatus = 1;
						for(j=0; j<nameList[k].vLen; j++){
							token = strtok(NULL, ", \t\n");
							if(token){
								switch(nameList[k].vType){
									case N_I:
										*NP_I = atol(token);
										break;
									case N_R:
										*NP_R = atof(token);
										break;
									case N_S:
										strcpy(NP_S,token);
										break;
								}
							}
							else{
								nameList[k].vStatus = 2;
								ok = 0;
							}
						}
						token = strtok(NULL, ", \t\n");
						if(token){
							nameList[k].vStatus = 3;
							ok = 0;
						}
						break;
					}
					else{
						nameList[k].vStatus = 4;
						ok = 0;
					}
			}
		}
		if(!match) ok = 0;
	}
	fclose(fp);
	for(k=0; k<sizeof(nameList)/sizeof(NameList); k++){
		if(nameList[k].vStatus != 1)
			ok = 1;
	}
	return(ok);
}/* end GetNameList() */

void PrintNameList(){
	int j, k;

	printf("===================================================\n");
	printf("-----------NameList for input data-----------------\n\n");
	for(k=0; k<sizeof(nameList)/sizeof(NameList); k++){
		printf("%s\t", nameList[k].vName);
		if(strlen(nameList[k].vName) < 8) printf("\t");
		if(nameList[k].vStatus > 0){
			for(j=0; j<nameList[k].vLen; j++){
				switch(nameList[k].vType){
					case N_I:
						printf("%d ", *NP_I);
						break;
					case N_R:
						printf("%#e ", *NP_R);
						break;
					case N_S:
						printf("%s ", NP_S);
						break;
				}
			}
		}
		switch(nameList[k].vStatus){
			case 0:
				printf("** no data");
				break;
			case 1:
				break;
			case 2:
				printf("** missing data");
				break;
			case 3:
				printf("** extra data");
				break;
			case 4:
				printf("** multiply defined");
				break;
		}
		printf("\n");
	}
	printf("\n-------------------End of NameList-----------------\n");
	printf("===================================================\n\n");

	return;
}/* end PrintNameList() */

void OpenFiles(void)
{
	if(mpirank!=0){
		PError("ERROR: Only mpirank#0 is allowed to modify global files",28);
	}
	fp_vm = fopen("vm.out","w+");
#ifdef DD_BASED_FLAG
	fp_dd = fopen("dd.out","w+");
#ifdef PF_DRX
	fp_drx = fopen("drx.out","w+");
#endif
#endif
	fp_err = fopen("err.out","w+");

	return;
}/*end OpenFiles()*/

void CloseFiles(void)
{
	if(mpirank!=0){
		PError("ERROR: Only mpirank#0 is allowed to modify global files",28);
	}
	fclose(fp_vm);
#ifdef DD_BASED_FLAG
	fclose(fp_dd);
#ifdef PF_DRX
	fclose(fp_drx);
#endif
#endif
	fclose(fp_err);

	return;
}/*end CloseFiles()*/

void WriteGrainMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	int i;

	sprintf(fname, "%s_S%06d.iout", s ,step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, mpirank*lsize*sizeof(int), MPI_INT,
			MPI_INT, "native", MPI_INFO_NULL);
	MPI_File_write(fp, grain_f, lsize, MPI_INT, &status);
	MPI_File_close(&fp);

	return;
}/*end WriteGrainMPI()*/

#ifdef PF_DRX
void WriteGrainPFMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	int i;

	sprintf(fname, "%s_S%06d.iout", s ,step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, mpirank*lsize_pf*sizeof(int), MPI_INT,
			MPI_INT, "native", MPI_INFO_NULL);
	MPI_File_write(fp, gID_rex, lsize_pf, MPI_INT, &status);
	MPI_File_close(&fp);

	return;
}/*end WriteGrainPFMPI()*/
#endif

void WriteDisgradMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	
	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = DisGrad[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = DisGrad[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = DisGrad[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = DisGrad[pIDX][0][1];
			}
		}
	
		sprintf(fname, "%s_%01d_S%06d.iout", s, i+1 ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, mpirank*lsize*sizeof(real), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
	}

	free(tmpVector);

	return;
}/*end WriteDisgradMPI()*/

void WriteElsMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	
	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = DisGrad[pIDX][i][i] - Eps[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = (DisGrad[pIDX][1][2]+DisGrad[pIDX][2][1])/2.0 - Eps[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = (DisGrad[pIDX][0][2]+DisGrad[pIDX][2][0])/2.0 - Eps[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = (DisGrad[pIDX][0][1]+DisGrad[pIDX][1][0])/2.0 - Eps[pIDX][0][1];
			}
		}
	
		sprintf(fname, "%s_%01d_S%06d.iout", s, i+1 ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, mpirank*lsize*sizeof(real), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
	}

	free(tmpVector);

  return;
}/* WriteElsMPI()*/

void WriteEpsMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	
	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = Eps[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = Eps[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = Eps[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = Eps[pIDX][0][1];
			}
		}
	
		sprintf(fname, "%s_%01d_S%06d.iout", s, i+1 ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, mpirank*lsize*sizeof(real), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
	}

	free(tmpVector);

	return;
}/*end WriteEpsMPI()*/

void WriteSigMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	
	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = Sig[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = Sig[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = Sig[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = Sig[pIDX][0][1];
			}
		}
	
		sprintf(fname, "%s_%01d_S%06d.iout", s, i+1 ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, mpirank*lsize*sizeof(real), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
	}

	free(tmpVector);

	return;
}/*end WriteSigMPI()*/

void WriteEdotMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	
	for(i=0;i<6;i++){
		local_loop{
			if(i<3){
				tmpVector[pIDX] = Edot[pIDX][i][i];
			}
			else if(i==3){
				tmpVector[pIDX] = Edot[pIDX][1][2];
			}
			else if(i==4){
				tmpVector[pIDX] = Edot[pIDX][0][2];
			}
			else if(i==5){
				tmpVector[pIDX] = Edot[pIDX][0][1];
			}
		}
	
		sprintf(fname, "%s_%01d_S%06d.iout", s, i+1 ,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, mpirank*lsize*sizeof(real), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
	}

	free(tmpVector);

	return;
}/*end WriteEdotMPI()*/

void WriteTextureMPI(char *s, int step)
{
	typedef struct{
		real angle[3];	// Euler angles
		int coord[3];		// coordinates
		int jgr;			// grain type
		int jph;			//phase type
		int buff;			// buffer zone
	} EntryType;
	EntryType value = {0};
	MPI_Datatype TexEntry;

	EntryType *tmpVector;
	real t1,ph,t2;
	ten2nd sa2xt;
	char fname[100] = {0};
	MPI_File fp;
	int i;

	int count = 5;	
	int block_lens[5] = {3,3,1,1,1};
	MPI_Aint indices[5];
	MPI_Datatype old_types[5] = {MPI_real,MPI_INT,MPI_INT,MPI_INT,MPI_INT};
	
	MPI_Address(&value, &indices[0]);
	MPI_Address(value.coord, &indices[1]);
	MPI_Address(&value.jgr, &indices[2]);
	MPI_Address(&value.jph, &indices[3]);
	MPI_Address(&value.buff, &indices[4]);
	// make relative
	for(i=count-1;i>=0;i--)
		indices[i] -= indices[0];
	MPI_Type_struct(count,block_lens,indices,old_types,&TexEntry);
	MPI_Type_commit(&TexEntry);

	AllocMem(tmpVector, lsize, EntryType);
	local_loop{
		T2_loop{
			sa2xt[mi][mj] = TranMat_xt2sa[pIDX][mj][mi];	// transpose
		}
		EulerToTransMatrix(&t1,&ph,&t2,sa2xt,1);
		tmpVector[pIDX].angle[0] = t1;
		tmpVector[pIDX].angle[1] = ph;
		tmpVector[pIDX].angle[2] = t2;
		tmpVector[pIDX].coord[0] = px+lxs+1;
		tmpVector[pIDX].coord[1] = py+1;
		tmpVector[pIDX].coord[2] = pz+1;
		tmpVector[pIDX].jgr = grain_f[pIDX];
		tmpVector[pIDX].jph = phase_f[pIDX];
		tmpVector[pIDX].buff = 0;
	}

	sprintf(fname, "%s_S%06d.iout", s, step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, mpirank*lsize*sizeof(EntryType), TexEntry,
			TexEntry, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, TexEntry, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}/*end WriteTextureMPI()*/

/* void WriteNewPositionMPI(char *s, int step)
{
	typedef struct{
		//real angle[3];	// Euler angles
		//int coord[3];		// coordinates
	//	int jgr;			// grain type
		//int jph;			//phase type
		//real shear;
        real SS1,SS2,SS3,SS4,SS5,SS6,SS7,SS8,SS9;
                // buffer zone
	} EntryType;
	EntryType value = {0};
	MPI_Datatype TexEntry;

	EntryType *tmpVector;
	//real t1,ph,t2;
	//ten2nd sa2xt;
	char fname[200] = {0};
	MPI_File fp;
	int i;

	int count = 9;	
	int block_lens[9] = {1,1,1,1,1,1,1,1,1};
	MPI_Aint indices[9];
MPI_Datatype old_types[9] = {MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real};
	
	//MPI_Address(&value.coord, &indices[0]);
	//MPI_Address(value.coord, &indices[1]);
//MPI_Address(&value.jgr, &indices[0]);
	//MPI_Address(&value.jph, &indices[3]);
        MPI_Address(&value.SS1, &indices[0]);
        MPI_Address(&value.SS2, &indices[1]);
        MPI_Address(&value.SS3, &indices[2]);
        MPI_Address(&value.SS4, &indices[3]);
        MPI_Address(&value.SS5, &indices[4]);
        MPI_Address(&value.SS6, &indices[5]);
        MPI_Address(&value.SS7, &indices[6]);
        MPI_Address(&value.SS8, &indices[7]);
        MPI_Address(&value.SS9, &indices[8]);
       
       
	// make relative
	for(i=count-1;i>=0;i--)
		indices[i] -= indices[0];
	MPI_Type_struct(count,block_lens,indices,old_types,&TexEntry);
	MPI_Type_commit(&TexEntry);

	AllocMem(tmpVector, lsize, EntryType);
	local_loop{
		//T2_loop{
		//	sa2xt[mi][mj] = TranMat_xt2sa[pIDX][mj][mi];	// transpose
		//}
		//EulerToTransMatrix(&t1,&ph,&t2,sa2xt,1);
		//tmpVector[pIDX].angle[0] = t1;
		//tmpVector[pIDX].angle[1] = ph;
		//tmpVector[pIDX].angle[2] = t2;
		//tmpVector[pIDX].coord[0] = px+lxs+1;
		//tmpVector[pIDX].coord[1] = py+1;
		//tmpVector[pIDX].coord[2] = pz+1;
	//	tmpVector[pIDX].jgr = grain_f[pIDX];
		//tmpVector[pIDX].jph = phase_f[pIDX];
                tmpVector[pIDX].SS1 = DisGrad[pIDX][0][0];
                tmpVector[pIDX].SS2 = DisGrad[pIDX][0][1];
                tmpVector[pIDX].SS3 = DisGrad[pIDX][0][2];
                tmpVector[pIDX].SS4 =  DisGrad[pIDX][1][0];
                tmpVector[pIDX].SS5 =  DisGrad[pIDX][1][1];
                tmpVector[pIDX].SS6 =  DisGrad[pIDX][1][2];
                tmpVector[pIDX].SS7 =  DisGrad[pIDX][2][0];
                tmpVector[pIDX].SS8 =  DisGrad[pIDX][2][1];
                tmpVector[pIDX].SS9 =  DisGrad[pIDX][2][2];
               
               
	}

	sprintf(fname, "%s_S%04d.iout", s, step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, mpirank*lsize*sizeof(EntryType), TexEntry,
			TexEntry, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, TexEntry, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}
*/
void WriteSLIPMPI(char *s, int step)
{
	typedef struct{
		//real angle[3];	// Euler angles
		//int coord[3];		// coordinates
		int jgr;			// grain type
		//int jph;			//phase type
		//real shear;
        real SS1,SS2,SS3,SS4,SS5,SS6,SS7,SS8,SS9,SS10,SS11,SS12;
                // buffer zone
	} EntryType;
	EntryType value = {0};
	MPI_Datatype TexEntry;

	EntryType *tmpVector;
	//real t1,ph,t2;
	//ten2nd sa2xt;
	char fname[200] = {0};
	MPI_File fp;
	int i;

	int count = 13;	
	int block_lens[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
	MPI_Aint indices[13];
	MPI_Datatype old_types[13] = {MPI_INT,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real};
	
	//MPI_Address(&value.coord, &indices[0]);
	//MPI_Address(value.coord, &indices[1]);
	MPI_Address(&value.jgr, &indices[0]);
	//MPI_Address(&value.jph, &indices[3]);
        MPI_Address(&value.SS1, &indices[1]);
        MPI_Address(&value.SS2, &indices[2]);
        MPI_Address(&value.SS3, &indices[3]);
        MPI_Address(&value.SS4, &indices[4]);
        MPI_Address(&value.SS5, &indices[5]);
        MPI_Address(&value.SS6, &indices[6]);
        MPI_Address(&value.SS7, &indices[7]);
        MPI_Address(&value.SS8, &indices[8]);
        MPI_Address(&value.SS9, &indices[9]);
        MPI_Address(&value.SS10, &indices[10]);
        MPI_Address(&value.SS11, &indices[11]);
        MPI_Address(&value.SS12, &indices[12]);
       
	// make relative
	for(i=count-1;i>=0;i--)
		indices[i] -= indices[0];
	MPI_Type_struct(count,block_lens,indices,old_types,&TexEntry);
	MPI_Type_commit(&TexEntry);

	AllocMem(tmpVector, lsize, EntryType);
	local_loop{
		//T2_loop{
		//	sa2xt[mi][mj] = TranMat_xt2sa[pIDX][mj][mi];	// transpose
		//}
		//EulerToTransMatrix(&t1,&ph,&t2,sa2xt,1);
		//tmpVector[pIDX].angle[0] = t1;
		//tmpVector[pIDX].angle[1] = ph;
		//tmpVector[pIDX].angle[2] = t2;
		//tmpVector[pIDX].coord[0] = px+lxs+1;
		//tmpVector[pIDX].coord[1] = py+1;
		//tmpVector[pIDX].coord[2] = pz+1;
		tmpVector[pIDX].jgr = grain_f[pIDX];
		//tmpVector[pIDX].jph = phase_f[pIDX];
                tmpVector[pIDX].SS1 = trial_gam_acum[pIDX][0];
                tmpVector[pIDX].SS2 =  trial_gam_acum[pIDX][1];
                tmpVector[pIDX].SS3 =  trial_gam_acum[pIDX][2];
                tmpVector[pIDX].SS4 =  trial_gam_acum[pIDX][3];
                tmpVector[pIDX].SS5 =  trial_gam_acum[pIDX][4];
                tmpVector[pIDX].SS6 =  trial_gam_acum[pIDX][5];
                tmpVector[pIDX].SS7 =  trial_gam_acum[pIDX][6];
                tmpVector[pIDX].SS8 =  trial_gam_acum[pIDX][7];
                tmpVector[pIDX].SS9 =  trial_gam_acum[pIDX][8];
                tmpVector[pIDX].SS10 =  trial_gam_acum[pIDX][9];
                tmpVector[pIDX].SS11 =  trial_gam_acum[pIDX][10];
                tmpVector[pIDX].SS12 =  trial_gam_acum[pIDX][11];
               
	}

	sprintf(fname, "%s_S%04d.iout", s, step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, mpirank*lsize*sizeof(EntryType), TexEntry,
			TexEntry, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, TexEntry, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}

/*end WriteSLIPMPI()*/

void WriteNewPositionMPI(char *s, int step)
{
	typedef struct{
		//real angle[3];	// Euler angles
		//int coord[3];		// coordinates
	//	int jgr;			// grain type
		//int jph;			//phase type
		//real shear;
        real SS1,SS2,SS3;
                // buffer zone
	} EntryType;
	EntryType value = {0};
	MPI_Datatype TexEntry;

	EntryType *tmpVector;
	//real t1,ph,t2;
	//ten2nd sa2xt;
	char fname[200] = {0};
	MPI_File fp;
	int i;

	int count = 3;	
	int block_lens[3] = {1,1,1};
	MPI_Aint indices[3];
	MPI_Datatype old_types[3] = {MPI_real,MPI_real,MPI_real};
	
	//MPI_Address(&value.coord, &indices[0]);
	//MPI_Address(value.coord, &indices[1]);
//MPI_Address(&value.jgr, &indices[0]);
	//MPI_Address(&value.jph, &indices[3]);
        MPI_Address(&value.SS1, &indices[0]);
        MPI_Address(&value.SS2, &indices[1]);
        MPI_Address(&value.SS3, &indices[2]);
       
       
	// make relative
	for(i=count-1;i>=0;i--)
		indices[i] -= indices[0];
	MPI_Type_struct(count,block_lens,indices,old_types,&TexEntry);
	MPI_Type_commit(&TexEntry);

	AllocMem(tmpVector, lsize, EntryType);
	local_loop{
		//T2_loop{
		//	sa2xt[mi][mj] = TranMat_xt2sa[pIDX][mj][mi];	// transpose
		//}
		//EulerToTransMatrix(&t1,&ph,&t2,sa2xt,1);
		//tmpVector[pIDX].angle[0] = t1;
		//tmpVector[pIDX].angle[1] = ph;
		//tmpVector[pIDX].angle[2] = t2;
		//tmpVector[pIDX].coord[0] = px+lxs+1;
		//tmpVector[pIDX].coord[1] = py+1;
		//tmpVector[pIDX].coord[2] = pz+1;
	//	tmpVector[pIDX].jgr = grain_f[pIDX];
		//tmpVector[pIDX].jph = phase_f[pIDX];
                tmpVector[pIDX].SS1 = displacement_fluct[pIDX][0];
                tmpVector[pIDX].SS2 = displacement_fluct[pIDX][1];
                tmpVector[pIDX].SS3 = displacement_fluct[pIDX][2];
               
               
	}

	sprintf(fname, "%s_S%04d.iout", s, step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, mpirank*lsize*sizeof(EntryType), TexEntry,
			TexEntry, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, TexEntry, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}/*end WriteSLIPMPI()*/

void WriteDDMPI(char *s, int step)
{
	typedef struct{
		//real angle[3];	// Euler angles
		//int coord[3];		// coordinates
		int jgr;			// grain type
		//int jph;			//phase type
		//real shear;
        real SS1,SS2,SS3,SS4,SS5,SS6,SS7,SS8,SS9,SS10,SS11,SS12;
                // buffer zone
	} EntryType;
	EntryType value = {0};
	MPI_Datatype TexEntry;

	EntryType *tmpVector;
	//real t1,ph,t2;
	//ten2nd sa2xt;
	char fname[200] = {0};
	MPI_File fp;
	int i;

	int count = 13;	
	int block_lens[13] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
	MPI_Aint indices[13];
	MPI_Datatype old_types[13] = {MPI_INT,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real,MPI_real};
	
	//MPI_Address(&value.coord, &indices[0]);
	//MPI_Address(value.coord, &indices[1]);
	MPI_Address(&value.jgr, &indices[0]);
	//MPI_Address(&value.jph, &indices[3]);
        MPI_Address(&value.SS1, &indices[1]);
        MPI_Address(&value.SS2, &indices[2]);
        MPI_Address(&value.SS3, &indices[3]);
        MPI_Address(&value.SS4, &indices[4]);
        MPI_Address(&value.SS5, &indices[5]);
        MPI_Address(&value.SS6, &indices[6]);
        MPI_Address(&value.SS7, &indices[7]);
        MPI_Address(&value.SS8, &indices[8]);
        MPI_Address(&value.SS9, &indices[9]);
        MPI_Address(&value.SS10, &indices[10]);
        MPI_Address(&value.SS11, &indices[11]);
        MPI_Address(&value.SS12, &indices[12]);
       
	// make relative
	for(i=count-1;i>=0;i--)
		indices[i] -= indices[0];
	MPI_Type_struct(count,block_lens,indices,old_types,&TexEntry);
	MPI_Type_commit(&TexEntry);

	AllocMem(tmpVector, lsize, EntryType);
	local_loop{
		//T2_loop{
		//	sa2xt[mi][mj] = TranMat_xt2sa[pIDX][mj][mi];	// transpose
		//}
		//EulerToTransMatrix(&t1,&ph,&t2,sa2xt,1);
		//tmpVector[pIDX].angle[0] = t1;
		//tmpVector[pIDX].angle[1] = ph;
		//tmpVector[pIDX].angle[2] = t2;
		//tmpVector[pIDX].coord[0] = px+lxs+1;
		//tmpVector[pIDX].coord[1] = py+1;
		//tmpVector[pIDX].coord[2] = pz+1;
		tmpVector[pIDX].jgr = grain_f[pIDX];
		//tmpVector[pIDX].jph = phase_f[pIDX];
                tmpVector[pIDX].SS1 = rho_s[pIDX][0];
                tmpVector[pIDX].SS2 = rho_s[pIDX][1];
                tmpVector[pIDX].SS3 = rho_s[pIDX][2];
                tmpVector[pIDX].SS4 = rho_s[pIDX][3];
                tmpVector[pIDX].SS5 = rho_s[pIDX][4];
                tmpVector[pIDX].SS6 = rho_s[pIDX][5];
                tmpVector[pIDX].SS7 = rho_s[pIDX][6];
                tmpVector[pIDX].SS8 = rho_s[pIDX][7];
                tmpVector[pIDX].SS9 = rho_s[pIDX][8];
                tmpVector[pIDX].SS10 = rho_s[pIDX][9];
                tmpVector[pIDX].SS11 = rho_s[pIDX][10];
                tmpVector[pIDX].SS12 = rho_s[pIDX][11];
               
	}

	sprintf(fname, "%s_S%04d.iout", s, step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, mpirank*lsize*sizeof(EntryType), TexEntry,
			TexEntry, "native", MPI_INFO_NULL);
	MPI_File_write(fp, tmpVector, lsize, TexEntry, &status);
	MPI_File_close(&fp);

	free(tmpVector);

	return;
}


#ifdef DD_BASED_FLAG
void WriteRhoMPI(char *s, char *type, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	int i;
	int jph;

	local_loop{
		jph = phase_f[pIDX];
		rho_tot[pIDX] = 0.0;
		/* SSD */
		if(type[0]=='s'||type[0]=='S'){
			for(i=0;i<nSYS[jph-1];i++){
				rho_tot[pIDX] += rho_s[pIDX][i];
			}
		}
		else if(type[0]=='m'||type[0]=='M'){
			for(i=0;i<nSYS[jph-1];i++){
				rho_tot[pIDX] += rho_m[pIDX][i];
			}
		}
		else if(type[0]=='d'||type[0]=='D'){
			for(i=0;i<nSYS[jph-1];i++){
				rho_tot[pIDX] += rho_dipole[pIDX][i];
			}
		}
		else if(type[0]=='c'||type[0]=='C'){
			for(i=0;i<nSYS[jph-1];i++){
				rho_tot[pIDX] += rho_forest[pIDX][i];
			}
		}
		else if(type[0]=='g'||type[0]=='G'){
			for(i=0;i<nSYS[jph-1];i++){
				rho_tot[pIDX] += sqrt(rho_g1[pIDX][i]*rho_g1[pIDX][i]+rho_g2[pIDX][i]*rho_g2[pIDX][i]+rho_g3[pIDX][i]*rho_g3[pIDX][i]);
			}
		}
		else{
			PError("Wrong type of dislcoations to write!",720);
		}
	}

	sprintf(fname, "%s_%s_S%06d.iout", s, type, step);
	MPI_File_open(MPI_COMM_WORLD, fname,
			MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
	MPI_File_set_view(fp, mpirank*lsize*sizeof(real), MPI_DOUBLE,
			MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_File_write(fp, rho_tot, lsize, MPI_DOUBLE, &status);
	MPI_File_close(&fp);

	return;
}/*end WriteRhoMPI()*/

void WriteRhoDotMPI(char *s, int step)
{
	char fname[100] = {0};
	MPI_File fp;
	real *tmpVector;
	int i;

	AllocMem(tmpVector, lsize, real);
	
	for(i=0;i<NSYSMX;i++){
    // SSD
		local_loop{
			tmpVector[pIDX] = rho_dot_s[pIDX][i];
		}
	
		sprintf(fname, "%s_SSD_slip%02d_S%06d.iout", s,i+1,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, mpirank*lsize*sizeof(real), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
    
    // GND-I
		local_loop{
			tmpVector[pIDX] = rho_dot_g1[pIDX][i];
		}
	
		sprintf(fname, "%s_GND1_slip%02d_S%06d.iout", s,i+1,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, mpirank*lsize*sizeof(real), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);
    
    // GND-II
		local_loop{
			tmpVector[pIDX] = rho_dot_g2[pIDX][i];
		}
	
		sprintf(fname, "%s_GND2_slip%02d_S%06d.iout", s,i+1,step);
		MPI_File_open(MPI_COMM_WORLD, fname,
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
		MPI_File_set_view(fp, mpirank*lsize*sizeof(real), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(fp, tmpVector, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&fp);

  }
	free(tmpVector);

	return;
}/*end WriteRhoDotMPI()*/



#endif
