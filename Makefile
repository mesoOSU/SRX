CC = mpicxx
HEADER = evp.h evpFlags.h V3math.h
CFLAGS = -Wall -O3
SEEKFLAG = -DMPICH_IGNORE_CXX_SEEK
#INCLUDES = -I/public/home/pzhao/Soft_mu01/include 
#LFLAGS = -L/public/home/pzhao/Soft_mu01/lib
MPI_CFLAGS = -I/users/PAS1064/osu10047/local/mpich/3.2.1/include
MPI_LIBS = -L/users/PAS1064/osu10047/local/mpich/3.2.1/lib
FFTW_CFLAGS = -I/users/PAS1064/osu10047/local/fftw/2.1.5/include
FFTW_LIBS = -L/users/PAS1064/osu10047/local/fftw/2.1.5/lib
LIBS = -lfftw_mpi -lfftw -lgsl -lgslcblas -lm
#LIBS = -lfftw3_mpi -lfftw3 -lgsl -lgslcblas -lm
SRCS = evp.c io.c init.c kinematics.c evolution.c constitutive.c
OBJS = $(SRCS:.c=.o)
MAIN = trial
.PHONY: depend clean
#all: $(MAIN)

$(MAIN): $(OBJS) Makefile
	$(CC) $(SEEKFLAG) $(FFTW_CFLAGS) $(MPI_CFLAGS)  $(CFLAGS)  -o $(MAIN) $(OBJS) $(FFTW_LIBS) $(MPI_LIBS) $(LIBS)

%.o: %.c $(HEADER) Makefile
	$(CC) $(SEEKFLAG) $(FFTW_CFLAGS) $(MPI_CFLAGS) $(CFLAGS)  -c $< -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $^
