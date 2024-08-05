#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#define N 4

double** Init_double(int r, int c)
{
    double** M;
    M = (double**)malloc(r * sizeof(double*));
    for (int i = 0; i < r; i++)
        M[i] = (double*)malloc(c * sizeof(double));
    return M;
}

int** Init_int( int r, int c)
{
    int**M;
    M = (int**)malloc(r * sizeof(int*));
    for (int i = 0; i < r; i++)
        M[i] = (int*)malloc(c * sizeof(int));
    return M;
}

int U_1(int x)
{
    return 1;
}

int U_2(int x)
{
    return 1;
}

typedef struct GRID_INFO_T
{
    int p; // # of processes
    MPI_Comm comm; // grid communicator
    MPI_Comm row_comm; // communicator for my row
    MPI_Comm col_comm; // communicator for my column
    int q; // rows in grid
    int c; //columns in grid
    int my_row; // my row’s coordinate
    int my_col; // my column’s coordinate
    int my_rank; // my rank in the grid
} GRID_INFO_T;


/* FUNCTION TO SETUP THE GRID */
void test_2D_grid(GRID_INFO_T* grid)
{
  int dimensions[2] ;   
  int wrap_around[2];
  int coordinates[2];
  int free_coords[2];
 
  /* set up global grid information */
  grid->q = 2;
  grid->c = (grid->p)/grid->q;
  dimensions[0] = grid->q;
  dimensions[1] = grid->c;

  /* circular shift in second dimension, also in first just because */
  wrap_around[0] = 0;
  wrap_around[1] = 0;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, 1, &grid->comm);
  MPI_Comm_rank(grid->comm, &grid->my_rank);

  /* get process coordinates in grid communicator */
  MPI_Cart_coords(grid->comm, grid->my_rank, 2, coordinates);
  grid->my_row = coordinates[0];
  grid->my_col = coordinates[1];
}

void Neighbor (int **A, int **M, int local_m, int local_n, int m, int n, int my_row, int my_col, int my_rank, int *b1, int *b2, int *b3, int *b4)
{
    for(int i=0;i<local_m;i++)
    {
        for(int j=0;j<local_n;j++)
        {
            int k=0;
            if(j==local_n-1)
            {
                A[k][i*local_n+j]=b1[i];
                k++;

            }
            else
            {
                A[k][i*local_n+j]=M[i][j+1];
                k++;
                
            }
            if(i==0)
            {
                A[k][i*local_n+j]=b2[j];
                
                k++;
            }
            else
            {
                A[k][i*local_n+j]=M[i-1][j];
                k++;
            }
            if(j==0)
            {
                A[k][i*local_n+j]=b3[i];
                k++;
            }
            else
            {
                A[k][i*local_n+j]=M[i][j-1];
                k++;
            }
            if(i==local_m-1)
            {
                A[k][i*local_n+j]=b4[j];
                k++;
            }
            else
            {
                A[k][i*local_n+j]=M[i+1][j];
                k++;
            }
        
        }
        
    }

    for(int i=0;i<4;i++)
        {
            for(int j=0;j<local_n*local_m;j++)
            {
                printf("%d     ",A[i][j]);
            }
            printf("\n");
        }

}



void Initialize(int **A, int local_m, int local_n, int m, int n, int my_row, int my_col, int my_rank)
{
  int *b1,*b2,*b3,*b4;
  b1 = (int*)malloc(local_m * sizeof(int));
  b3 = (int*)malloc(local_m * sizeof(int));

  b2 =(int*)malloc(local_n * sizeof(int));
  b4 = (int*)malloc(local_n * sizeof(int));

  int **M;
  M=Init_int(local_m,local_n);

  int C,D=1;
  int index[local_m*local_n];
  MPI_Status status;
  C = 1+(my_col*(local_n))+(my_row*local_m*N);

  printf("\np = %d\n",my_rank);
  for(int i=0; i<local_m; i++)
  {
    for(int j=0; j<local_n;j++)
    {
        M[i][j]=C;
        C++;
    }
    C=C+(n-1)*local_n;
  }

  if(my_col==n-1)
  {
    for(int i=0;i<local_m; i++)
    {
       
        MPI_Send(&M[i][0],1,MPI_INT, my_rank-1,1,MPI_COMM_WORLD);
        MPI_Recv(&b3[i], 1, MPI_INT, my_rank-1,3,MPI_COMM_WORLD, &status);
    }
    for(int i=0;i<local_m; i++)
    {
        MPI_Send(&M[i][local_n-1], 1, MPI_INT,my_rank-n+1, 2, MPI_COMM_WORLD );
        MPI_Recv(&b1[i],1,MPI_INT,my_rank-n+1,1,MPI_COMM_WORLD,&status);
    }
  }

  else if(my_col==0)
  {
    for(int i=0;i<local_m; i++)
    {
        MPI_Recv(&b1[i], 1, MPI_INT, my_rank+1,1,MPI_COMM_WORLD, &status);
        MPI_Send(&M[i][local_n-1],1,MPI_INT, my_rank+1,3,MPI_COMM_WORLD);
    
    }
    for(int i=0;i<local_m; i++)
    {
        MPI_Recv(&b3[i], 1, MPI_INT, my_rank+n-1,2,MPI_COMM_WORLD, &status);
        MPI_Send(&M[i][0],1,MPI_INT, my_rank+n-1,1,MPI_COMM_WORLD);
    
    }

  }
  else
  {
    for(int i=0;i<local_m; i++)
    {
        MPI_Recv(&b1[i], 1, MPI_INT, my_rank+1,1,MPI_COMM_WORLD, &status);
        MPI_Send(&M[i][local_n-1],1,MPI_INT, my_rank+1,3,MPI_COMM_WORLD);   
    }
    for(int i=0;i<local_m; i++)
    {    
        MPI_Send(&M[i][0],1,MPI_INT, my_rank-1,1,MPI_COMM_WORLD);
        MPI_Recv(&b3[i], 1, MPI_INT, my_rank-1,3,MPI_COMM_WORLD, &status);
    }
  }

  if(my_row == 0)
  {
    MPI_Send(&M[0][0], local_n, MPI_INT, n*(m-1)+my_rank, 4, MPI_COMM_WORLD );
    MPI_Recv(&b2[0],local_n, MPI_INT, n*(m-1)+my_rank, 2, MPI_COMM_WORLD,&status);

    MPI_Send(&M[local_m-1][0], local_n, MPI_INT, n+my_rank, 2, MPI_COMM_WORLD );
    MPI_Recv(&b4[0], local_n, MPI_INT, n+my_rank, 4, MPI_COMM_WORLD,&status );
  }
  else if(my_row == m-1 )
  {
    MPI_Send(&M[local_m-1][0], local_n, MPI_INT, my_rank-n*(m-1), 2, MPI_COMM_WORLD );
    MPI_Recv(&b4[0], local_n, MPI_INT, my_rank-n*(m-1),4,MPI_COMM_WORLD, &status);
    
    MPI_Send(&M[0][0], local_n, MPI_INT, my_rank-n, 4, MPI_COMM_WORLD );
    MPI_Recv(&b2[0], local_n, MPI_INT, my_rank-n,2,MPI_COMM_WORLD, &status);
  }
  else
  {
    MPI_Send(&M[0][0], local_n, MPI_INT, my_rank-n, 4, MPI_COMM_WORLD );
    MPI_Recv(&b2[0], local_n, MPI_INT, my_rank-n,2,MPI_COMM_WORLD, &status);

    MPI_Send(&M[local_m-1][0], local_n, MPI_INT, n+my_rank, 2, MPI_COMM_WORLD );
    MPI_Recv(&b4[0], local_n, MPI_INT, n+my_rank, 4, MPI_COMM_WORLD,&status );
  }
    Neighbor(A,M,local_m,local_n,m,n,my_row,my_col,my_rank,b1,b2,b3,b4);

}


/*void print_dirac(int **A, int local_m, int local_n)
{
    int i, j, k,l;
    l=0;
    double **Dirac,**v, **Dirac_matrix;
    Dirac=Init_double(local_m-2,local_n-2);
    v= Init_double((local_m-2)*(local_n-2),1);
    Dirac_matrix= Init_double((local_m-2)*(local_n-2),(local_m-2)*(local_n-2));
    
    for (i = 0; i < (local_m-2)*(local_n-2); i++) 
    {
        for (j = 0; j < (local_m-2)*(local_n-2); j++)
        {

            if(i==j)
            {
                v[j][0]=1;
            }
            else
            {
                v[j][0]=0;
            }
            
        }
        Dirac= dirac(v, A, local_m, local_n);
        
        for(k=0;k<(local_m-2)*(local_n-2);k++)
        {
             Dirac_matrix[k][l]=Dirac[k][0];
        //printf("V = %f", sol.M[][0]);
        }
        l++;              
    }

    for(k=0;k<(local_m-2)*(local_n-2); k++)
    {
        
        for (l = 0; l < (local_m-2)*(local_n-2); l++)
        {  
            printf("%2.2f  ",Dirac_matrix[k][l]);
        }
        printf("\n");
        
        
        }
}*/



 /****************   MAIN   *****************/
int main(int argc, char* argv[])
{
  int local_m, local_n,C=0;
  int my_rank;

  MPI_Init(&argc, &argv);
  struct GRID_INFO_T *g;
  g =(struct GRID_INFO_T*) malloc(sizeof(struct GRID_INFO_T));
  MPI_Comm_size (MPI_COMM_WORLD, &g->p );
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  test_2D_grid(g);

  local_n=2;
  local_m = (N*N)/(g->p*local_n);



  int **A;
  A= Init_int(4,(local_m)*(local_n));
  //double **Dirac;
  //Dirac=Init_double((local_m-2)*(local_n-2),1);
  

  Initialize(A,local_m,local_n,g->q,g->c,g->my_row,g->my_col,g->my_rank);



  /*printf("myArray from processor %d, %d, %d = \n",my_rank,g->my_row,g->my_col);
  for (int k = 0; k < 4; k++) 
    {
        
        for (int l = 0; l < local_n*local_m; l++) 
        {
            printf("%d ",A[k][l]);
        }
        printf("\n");
    }*/
    

  MPI_Finalize();
  return 0;
}
