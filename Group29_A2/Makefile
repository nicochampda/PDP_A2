CC=mpicc
CFLAGS=-O3 
LIBS= -lmpi -lm

MPI_wave: MPI_wave.c
	$(CC) $(CFLAGS) -o MPI_wave MPI_wave.c $(LIBS)


clean:
	$(RM) MPI_wave
