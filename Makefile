CC=mpicc
CFLAGS=-Wall -O3 -g
LIBS= -lmpi -lm

MPI_wave: MPI_wave.c
	$(CC) $(CFLAGS) -o MPI_wave MPI_wave.c $(LIBS)

wave: wave.c
	gcc $(CFLAGS) $< -o $@ -lm

clean:
	$(RM) MPI_wave
	$(RM) wave

