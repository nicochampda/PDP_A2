#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <unistd.h>

//#define VERIFY  

void Prvalues(int length, int heigth,  double matrix[length * heigth]);
double timer();
double initialize(double x, double y, double t);
void save_solution(double *u, int Ny, int Nx, int n);

int main(int argc, char *argv[]){
    int Nx,Ny,Nt;
    double dt, dx, lambda_sq;
    double x, y;
    double *u;
    double *u_old;
    double *u_new;
    double *u_old_blocks;
    double *u_blocks;
    double *u_new_blocks;
    double begin,end;
    int nprocs, rank;

  
    int p1 = 0, p2 = 0; // number of x-tiles and y-tiles meaning there are p1*p2 processors

    if(argc < 4){
        printf("usage : mpirun -np p1*p2 ./MPI_wave N p1 p2\n");
        return -1;
    }

    //reading inputs
    Nx=atoi(argv[1]);
    p1=atoi(argv[2]);
    p2=atoi(argv[3]);

    Ny=Nx;
    Nt=Nx;
    dx=1.0/(Nx-1);
    dt=0.50*dx;
    lambda_sq = (dt/dx)*(dt/dx);


    //calcul of the length and heigth of a block even if N non divisible
    int block_heigth = (Ny - 2)/p2 + 2; 
    int block_length = (Nx- 2)/p1 + 2;

    int mod_heigth = (Ny - 2)%p2; 
    int mod_length = (Nx - 2)%p1; 
    if(mod_heigth != 0){
        block_heigth += 1;
        mod_heigth = p2 - mod_heigth;
    }

    if(mod_length != 0){
        block_length += 1;
        mod_length = p1 - mod_length;
    }


    MPI_Request req_u_send;
    MPI_Request req_u_recv;
    MPI_Request req_u_new_send;
    MPI_Request req_u_new_recv;

    MPI_Request req_shift[8];



    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);


    if(rank == 0){
        if(p1*p2 != nprocs){
            printf("the number of procs use must be equal to p1 * p2 \n");
            return -1;
        }

        printf("Approximation of the solution of the wave equation with MPI\n");
        printf("Using %d processors. \n",nprocs);
        printf("Using a total of %d points and time steps \n", Nx);
        printf("\n");

    }


    MPI_Comm old_comm = MPI_COMM_WORLD;
	int ndims = 2;
	int dim_size[2];
	dim_size[0] = p1;
	dim_size[1] = p2;
	int periods[2];
	periods[0] = 0;
	periods[1] = 0;
	int reorder = 1;

	int row_rank, col_rank, temp_rank;

	int coords[2];
    int temp_coords[2];

	MPI_Comm grid_comm, row_comm, col_comm;
        

	// Create cartesian communicator
	MPI_Cart_create (old_comm , ndims , dim_size, periods, reorder, &grid_comm);

	// Set new rank and coordinates
	MPI_Comm_rank(grid_comm, &rank);
	MPI_Cart_coords(grid_comm, rank, ndims, coords);

	MPI_Comm_split(grid_comm, coords[0], coords[1], &row_comm);
	MPI_Comm_rank(row_comm, &col_rank);

	MPI_Comm_split(grid_comm, coords[1], coords[0], &col_comm);
	MPI_Comm_rank(col_comm, &row_rank);

    /* Define new type to send datas */
	MPI_Datatype blocktype, blockselect, rowtype, coltype;

	MPI_Type_contiguous(block_heigth * block_length, MPI_DOUBLE, &blocktype);
    MPI_Type_commit(&blocktype);

	MPI_Type_vector(block_heigth, block_length, Nx, MPI_DOUBLE, &blockselect);
	MPI_Type_commit(&blockselect);

    MPI_Type_contiguous(block_length, MPI_DOUBLE, &rowtype);
    MPI_Type_commit(&rowtype);

	MPI_Type_vector(block_heigth, 1, block_length, MPI_DOUBLE, &coltype);
	MPI_Type_commit(&coltype);

    //Set mod to 0 in all procs except those on the last column/row
    if(row_rank != p1-1 ){
        mod_length = 0;
    }

    if(col_rank != p2-1 ){
        mod_heigth = 0;
    }

    if(rank == 0){
        u = malloc(Nx*Ny*sizeof(double));
        u_old = malloc(Nx*Ny*sizeof(double));
        u_new = malloc(Nx*Ny*sizeof(double));

        /* Setup IC */

        memset(u,0,Nx*Ny*sizeof(double));
        memset(u_old,0,Nx*Ny*sizeof(double));
        memset(u_new,0,Nx*Ny*sizeof(double));


        /* Initialization p */
        for(int i = 1; i < (Ny-1); ++i) {
            for(int j = 1; j < (Nx-1); ++j) {
                x = j*dx;
                y = i*dx;

                /* u0 */
                u[i*Nx+j] = initialize(x,y,0);

                /* u1 */
                u_new[i*Nx+j] = initialize(x,y,dt);
           }
        }


        /* Send all parts of u_old, u and u_new to the right processors */
        for(int i = 0; i < p1; i++) {
            for(int j = 0; j < p2; j++) {
           
                temp_coords[0] = i ;
                temp_coords[1] = j ;

                MPI_Cart_rank(grid_comm,temp_coords,&temp_rank);

                MPI_Isend(&u[i*(block_length-2) + j*(block_heigth-2)*Nx], 1, blockselect, temp_rank, 0, grid_comm, &req_u_send);
                MPI_Isend(&u_new[i*(block_length-2) + j*(block_heigth-2)*Nx], 1, blockselect, temp_rank, 1, grid_comm, &req_u_new_send);

            }
        } 
    }

#ifdef VERIFY
        double max_error=0.0;
        double local_max_error = 0.0;
#endif



    /* Each processor receive block of u0 and u1 as initialization */
    u_blocks = (double *)malloc(block_heigth * block_length * sizeof(double));
    u_new_blocks = (double *)malloc(block_heigth * block_length * sizeof(double));
    u_old_blocks = (double *)malloc(block_heigth * block_length * sizeof(double));

    memset(u_old_blocks,1,block_heigth*block_length*sizeof(double));

    MPI_Irecv(u_blocks, 1, blocktype, 0, 0, grid_comm, &req_u_recv);
    MPI_Irecv(u_new_blocks, 1, blocktype, 0, 1, grid_comm, &req_u_new_recv);

    MPI_Wait(&req_u_recv, MPI_STATUS_IGNORE);
    MPI_Wait(&req_u_new_recv, MPI_STATUS_IGNORE);


    /* integration of the solution on each processors */

    begin=MPI_Wtime();

    for(int n=2; n<Nt; ++n) {
        // Swap ptrs 
        double *tmp_blocks = u_old_blocks;
        u_old_blocks = u_blocks;
        u_blocks = u_new_blocks;
        u_new_blocks = tmp_blocks;

        // Apply stencil
        for(int i = 1; i < (block_heigth-1-mod_heigth); i++) {
            for(int j = 1; j < (block_length-1-mod_length); j++) {
                u_new_blocks[i*block_length+j] = 2*u_blocks[i*block_length+j] - u_old_blocks[i*block_length+j] + lambda_sq*(u_blocks[(i+1)*block_length+j] + u_blocks[(i-1)*block_length+j] + u_blocks[i*block_length+j+1] + u_blocks[i*block_length+j-1] - 4*u_blocks[i*block_length+j]);
        

#ifdef VERIFY
                double e = fabs(u_new_blocks[i*block_length + j] - initialize((j+row_rank*(block_length-2))*dx, (i+col_rank*(block_heigth-2))*dx, n*dt));
                if(e>local_max_error){
                      local_max_error = e;
                }
#endif

            }
        }
        
        /* exchange ghost values */
        int source, dest;
        //sending left column
        MPI_Cart_shift(grid_comm, 0, 1, &source, &dest);
        MPI_Isend(&u_new_blocks[1], 1, coltype, source, 11, grid_comm, &req_shift[0]);
        MPI_Irecv(&u_new_blocks[block_length-1], 1, coltype, dest, 11, grid_comm, &req_shift[1]);

        //sending right column
        MPI_Cart_shift(grid_comm, 0, -1, &source, &dest);
        MPI_Isend(&u_new_blocks[block_length-2], 1, coltype, source, 22, grid_comm, &req_shift[2]);
        MPI_Irecv(&u_new_blocks[0], 1, coltype, dest, 22, grid_comm, &req_shift[3]);

        //sending bottom row
        MPI_Cart_shift(grid_comm, 1, -1, &source, &dest);
        MPI_Isend(&u_new_blocks[(block_heigth - 2)*block_length], 1, rowtype, source, 33, grid_comm, &req_shift[4]);
        MPI_Irecv(&u_new_blocks[0], 1, rowtype, dest, 33, grid_comm, &req_shift[5]);

        //sending top row
        MPI_Cart_shift(grid_comm, 1, 1, &source, &dest);
        MPI_Isend(&u_new_blocks[block_length], 1, rowtype, source, 44, grid_comm, &req_shift[6]);
        MPI_Irecv(&u_new_blocks[(block_heigth-1)*block_length], 1, rowtype, dest, 44, grid_comm, &req_shift[7]);

        MPI_Waitall(8, req_shift, MPI_STATUS_IGNORE);

    }
    end = MPI_Wtime();
    if(rank == 0){
        printf("Time : %f\n", end - begin);
    }

#ifdef VERIFY
    MPI_Reduce(&local_max_error, &max_error, 1, MPI_DOUBLE, MPI_MAX, 0, grid_comm);
    if(rank == 0){
        printf("Max error : %f\n", max_error);
    }
#endif
  
//        sleep(1+rank);
//        printf("coords %i %i\n", row_rank, col_rank);
//        Prvalues(block_length, block_heigth, u_new_blocks);
 

    MPI_Type_free(&blockselect);
    MPI_Type_free(&blocktype);
    MPI_Type_free(&rowtype);
    MPI_Type_free(&coltype);

    MPI_Comm_free(&row_comm);   
    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&grid_comm);
  

    MPI_Finalize();

    return 0;
}


double timer()
{ 
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

double initialize(double x, double y, double t)
{
  double value = 0;
#ifdef VERIFY
  /* standing wave */
  value=sin(3*M_PI*x)*sin(4*M_PI*y)*cos(5*M_PI*t);
#else
  /* squared-cosine hump */
  const double width=0.1;

  double centerx = 0.25;
  double centery = 0.5;

  double dist = sqrt((x-centerx)*(x-centerx) +
                     (y-centery)*(y-centery));
  if(dist < width) {
    double cs = cos(M_PI_2*dist/width);
    value = cs*cs;
  }
#endif
  return value;
}

void save_solution(double *u, int Ny, int Nx, int n)
{
  char fname[50];
  sprintf(fname,"solution-%d.dat",n);
  FILE *fp = fopen(fname,"w");

  fprintf(fp,"%d %d\n",Nx,Ny);

  for(int j = 0; j < Ny; ++j) {
    for(int k = 0; k < Nx; ++k) {
      fprintf(fp,"%e\n",u[j*Nx+k]);
    }
  }

  fclose(fp);
}


void Prvalues(int length, int heigth,  double matrix[length * heigth]){   
    int i, j;
    printf("\n");
    for (i = 0; i < heigth; i++){
        for (j = 0; j < length; j++){
            printf("%.3f\t", matrix[i*length + j]);
        }
        printf("\n");
    }
    printf("\n");
}
