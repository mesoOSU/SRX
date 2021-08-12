#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

int mpirank = 0, NumPE = 0;
int lsize = 0;
FILE *fp;
MPI_Status status;
MPI_File FP;
MPI_Offset filesize;

int main(int argc, char **argv){
	int i, fsize;
	double buffer = 0.0;
	double *array;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NumPE);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

	assert(argc==4);
	if(argv[3][0]=='f'){	// from out to iout
		fp = fopen(argv[1], "r");
		if(fp==NULL){
			perror("Error opening file\n");
			exit(3);
		}
		fsize = 0;
		while(1){
			fsize++;
			if(fscanf(fp, "%lf", &buffer) ==  EOF){
				fsize--;
				break;
			}
		}
		rewind(fp);
	
		lsize = fsize/NumPE;
		array = (double*)malloc(lsize*sizeof(double));
		for(i=0;i<mpirank*lsize;i++){
			fscanf(fp, "%lf", &buffer);
		}
	
		for(i=0;i<lsize;i++){
			if(fscanf(fp, "%lf\t", array+i)!=1){
				perror("Error of reading data file\n");
				exit(7);
			}
		}
		fclose(fp);
	
		MPI_File_open(MPI_COMM_WORLD, argv[2],
				MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &FP);
		MPI_File_set_view(FP, mpirank*lsize*sizeof(double), MPI_DOUBLE,
				MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_write(FP, array, lsize, MPI_DOUBLE, &status);
		MPI_File_close(&FP);
	
		free(array);

	}
	else if(argv[3][0]=='b'){	// from iout to out

		assert(NumPE==1);	// multiple PEs requires data communication
		if(mpirank==0){
			MPI_File_open(MPI_COMM_WORLD, argv[2],
					MPI_MODE_RDONLY, MPI_INFO_NULL, &FP);
			MPI_File_get_size(FP, &filesize);
			filesize = filesize/sizeof(double);
			array = (double*)malloc(filesize*sizeof(double));
			MPI_File_set_view(FP, 0, MPI_DOUBLE,
					MPI_DOUBLE, "native", MPI_INFO_NULL);
			MPI_File_read(FP, array, filesize, MPI_DOUBLE, &status);
			MPI_File_close(&FP);
	
			fp = fopen(argv[1], "w");
			for(i=0;i<filesize;i++)
				fprintf(fp, "%f\n", array[i]);
			fclose(fp);

			free(array);
		}
	}

	MPI_Finalize();

	return 0;
}
