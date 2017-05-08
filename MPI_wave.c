#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <unistd.h>

#define WRITE_TO_FILE 
#define VERIFY 

void Prvalues(int length, int heigth,  double matrix[length * heigth]);
double timer();
double initialize(double x, double y, double t);
void save_solution(double *u, int Ny, int Nx, int n);

int main(int argc, char *argv[]){
    int Nx,Ny,Nt;
    double dt, dx, lambda_sq;
    double x, y ;
    double *u;
    double *u_old;
    double *u_new;
    double *u_old_blocks;
    double *u_blocks;
    double *u_new_blocks;
    double begin,end;
    int nprocs, rank, i_global, j_global, i_local, j_local, i_local_min, i_local_max, j_local_min, j_local_max ;

  
    int p1,p2 ; // number of x-tiles and y-tiles meaning there are p1*p2 processors


    p1=atoi(argv[2]);
    p2=atoi(argv[3]);

    nprocs=p1*p2;
 
    Nx=128;
    if(argc>1)
        Nx=atoi(argv[1]);
    Ny=Nx;
    Nt=Nx;
    dx=1.0/(Nx-1);
    dt=0.50*dx;
    lambda_sq = (dt/dx)*(dt/dx);

    int block_length = (Ny - 2)/p2 + 2; 
    int block_heigth = (Nx- 2)/p1 + 2;

    printf("len %i hei %i\n", block_length, block_heigth);

    MPI_Request req_u_send;
    MPI_Request req_u_recv;
    MPI_Request req_u_new_send;
    MPI_Request req_u_new_recv;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    printf("Approximation of the solution of the wave equation with MPI\n");
    printf("Using %d processors. \n",nprocs);
    printf("Using a total of %d points and time steps \n", Nx);
    printf("\n");




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
	MPI_Datatype blocktype, blockselect;

	MPI_Type_contiguous(block_length * block_heigth, MPI_DOUBLE, &blocktype);
    MPI_Type_commit(&blocktype);

	MPI_Type_vector(block_length, block_heigth, Nx, MPI_DOUBLE, &blockselect);
	MPI_Type_commit(&blockselect);

    printf("coords %i %i rank %i row %i\n ", coords[0], coords[1], rank, row_rank);
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
                u[i*Nx+j] = i+ 0.1*j; //initialize(x,y,0);

                /* u1 */
                u_new[i*Nx+j] = i+ 0.1*j; //initialize(x,y,dt);
           }
        }


        Prvalues(Nx, Ny, u);
        /* Send all parts of u_old, u and u_new to the right processors */
        for(int i = 0; i < p1; i++) {
            for(int j = 0; j < p2; j++) {
           
                temp_coords[0] = i ;
                temp_coords[1] = j ;

                MPI_Cart_rank(grid_comm,temp_coords,&temp_rank);
                printf("send to %i %i %i\n", i, j, temp_rank);
                MPI_Isend(&u[i*(block_heigth-2) + j*(block_length-2)*Nx], 1, blockselect, temp_rank, 0, grid_comm, &req_u_send);
                MPI_Isend(&u_new[i*(block_heigth-2) + j*(block_length-2)*Nx], 1, blockselect, temp_rank, 1, grid_comm, &req_u_new_send);
            }
        } 
    }


    /* Each processor receive block of u0 and u1 as initialization */
    u_blocks = (double *)malloc(block_length * block_heigth * sizeof(double));
    u_new_blocks = (double *)malloc(block_length * block_heigth * sizeof(double));
    u_old_blocks = (double *)malloc(block_length * block_heigth * sizeof(double));

    MPI_Irecv(u_blocks, 1, blocktype, 0, 0, grid_comm, &req_u_recv);
    MPI_Irecv(u_new_blocks, 1, blocktype, 0, 1, grid_comm, &req_u_new_recv);

    MPI_Wait(&req_u_recv, MPI_STATUS_IGNORE);
    MPI_Wait(&req_u_new_recv, MPI_STATUS_IGNORE);


    /* Printing part */
    sleep(rank);
    printf("coords %i %i\n", row_rank, col_rank);
    Prvalues(block_heigth, block_length, u_blocks);
//    Prvalues(block_heigth, block_length, u_new_blocks);
  
    /*
    each processors receive its right part of u_old, u, u_new from the root processor (plus the halo points already?) and store them in variables called u_old_local, u_local, u_new_local 
    */

    /* integration of the solution on each processors */

    for(int n=2; n<Nt; ++n) {
        // Swap ptrs 
        double *tmp_blocks = u_old_blocks;
        u_old_blocks = u_blocks;
        u_blocks = u_new_blocks;
        u_new_blocks = tmp_blocks;

        // Apply stencil 
        for(int i = 1; i < (block_heigth-1); ++i) {
            for(int j = 1; j < (block_length-1); ++j) {
                u_new_blocks[i*Nx+j] = 2*u_blocks[i*Nx+j] - u_old_blocks[i*Nx+j] + lambda_sq*(u_blocks[(i+1)*Nx+j] + u_blocks[(i-1)*Nx+j] + u_blocks[i*Nx+j+1] + u_blocks[i*Nx+j-1] - 4*u_blocks[i*Nx+j]);
        
            }
        }
        // Exchange ghost values
        // if processor not on the rigth/left --> send column
        if (col_rank =! 0){
            //send rigth column
        }
        if (col_rank =! p1){
            //send left column
        }
        // if processor not on the top/bottom --> send row
        if (row_rank =! 0){
            //send rigth column
        }
        if (row_rank =! p2){
            //send left column
        }
    }


  
 

    MPI_Type_free(&blockselect);
    MPI_Type_free(&blocktype);

    MPI_Comm_free(&row_comm);   
    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&grid_comm);
  

    MPI_Finalize();

    return 0;
}

/********************************************************/
/*This part below is taken from the serial part wave.c */
/******************************************************/

/*  for(int i = 1; i < (Ny-1); ++i) {
    for(int j = 1; j < (Nx-1); ++j) {
      double x = j*dx;
      double y = i*dx;

      // u0 
      u[i*Nx+j] = initialize(x,y,0);

      // u1 
      u_new[i*Nx+j] = initialize(x,y,dt);
    }
  }

#ifdef WRITE_TO_FILE
  save_solution(u_new,Ny,Nx,1);
#endif
#ifdef VERIFY
  double max_error=0.0;
#endif

  // Integrate 

  begin=timer();
  for(int n=2; n<Nt; ++n) {
    // Swap ptrs 
    double *tmp = u_old;
    u_old = u;
    u = u_new;
    u_new = tmp;

    // Apply stencil 
    for(int i = 1; i < (Ny-1); ++i) {
      for(int j = 1; j < (Nx-1); ++j) {

        u_new[i*Nx+j] = 2*u[i*Nx+j]-u_old[i*Nx+j]+lambda_sq*(u[(i+1)*Nx+j] + u[(i-1)*Nx+j] + u[i*Nx+j+1] + u[i*Nx+j-1] -4*u[i*Nx+j]);
      }
    }

#ifdef VERIFY
    double error=0.0;
    for(int i = 0; i < Ny; ++i) {
      for(int j = 0; j < Nx; ++j) {
        double e = fabs(u_new[i*Nx+j]-initialize(j*dx,i*dx,n*dt));
        if(e>error)
          error = e;
      }
    }
    if(error > max_error)
      max_error=error;
#endif

#ifdef WRITE_TO_FILE
    save_solution(u_new,Ny,Nx,n);
#endif

  }
  end=timer();

  printf("Time elapsed: %g s\n",(end-begin));

#ifdef VERIFY
  printf("Maximum error: %g\n",max_error);
#endif

  free(u);
  free(u_old);
  free(u_new);

  return 0;
}
*/

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
            printf("%.1f\t", matrix[i*length + j]);
        }
        printf("\n");
    }
    printf("\n");
}
