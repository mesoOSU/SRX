#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

typedef double real;
#define MPI_real MPI_DOUBLE
int rank = 0, NumPE = 0;
int lsize = 0;
FILE *fp;
MPI_Status status;
MPI_File FP;
MPI_Offset filesize;

int main(int argc, char **argv)
{
	typedef struct{
		real angle[3];	// Euler angles
		int coord[3];		// coordinates
		int jgr;			// grain type
		int jph;			//phase type
		int buff;			// buffer zone
	} EntryType;
	EntryType value;
	EntryType *array;
	MPI_Datatype TexEntry;
	int i, fsize;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NumPE);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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


	assert(argc==3);
	assert(NumPE==1);	// multiple PEs requires data communication
	if(rank==0){
		MPI_File_open(MPI_COMM_WORLD, argv[1],
				MPI_MODE_RDONLY, MPI_INFO_NULL, &FP);
		MPI_File_get_size(FP, &filesize);
		filesize /= sizeof(EntryType);
		array = (EntryType*)malloc(filesize*sizeof(EntryType));
		MPI_File_set_view(FP, 0, TexEntry,
				TexEntry, "native", MPI_INFO_NULL);
		MPI_File_read(FP, array, filesize, TexEntry, &status);
		MPI_File_close(&FP);

		fp = fopen(argv[2], "w");
		for(i=0;i<filesize;i++)
			fprintf(fp, "%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\n",
					array[i].angle[0],array[i].angle[1],array[i].angle[2],
					array[i].coord[0],array[i].coord[1],array[i].coord[2],
					array[i].jgr, array[i].jph);
		fclose(fp);

		free(array);
	}

	return 0;
}

