#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#define N 50

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


double **Dirac(double **v, int **A, int *index, int local_m, int local_n, int my_row, int my_col, int m, int n)
{
    int j=0,x1;
    double **Sol;
    Sol = Init_double(local_m*local_n,1);
    
        int x,y,flag=0,count=1,in;
        float del,del1y,del_1y,del2y,del_2y;
        float mu =0.4;
        double S=0, *D=(double*)malloc(N*N * sizeof(double));
        if(local_m%2!=0)
        {
            if(my_row%2==0)
            x1=1;
            else
            x1=2;
        }
        else
            x1 = 1;
       
        for(int i=1; i<=local_m*local_n;i++)
        {
            in= index[i-1];
            
            for(int j=1; j<=N*N;j++)
            {
                flag=0;
                del=0; 
                del1y=0;
                del_1y=0;
                del2y=0;
                del_2y=0;

                if(in==j)
                {
                    del=1;
                    flag=1;
                }
 
                else if(j==abs((A[0][i-1])))
                {
                    if((in)%N==0)
                    {
                        del2y=-1;

                    }
                    else
                        del2y=1;
                    flag=1;
                }
                else if(j==abs(A[1][i-1]))
                {
                    if(in<=N)
                    {
                        del1y=-1;

                    }
                    else
                        del1y=1;
                    flag=1;
                }
                else if(j==abs(A[2][i-1]))
                {
                    if((in-1)%N==0)
                        del_2y=-1;
                    else
                        del_2y=1;
                    flag=1;
                }
                else if(j==abs(A[3][i-1]))
                {
                    if(in>N*(N-1))
                        del_1y=-1;
                    else
                        del_1y=1;
                    flag=1;
                }

                if(flag==1)
                {
                    
                    D[j-1]=mu*del+(U_1(index[i-1])*del1y-U_1(A[3][i-1])*del_1y)/2+pow((-1),x1)*(U_2(index[i-1])*del2y-U_2(A[2][i-1])*del_2y)/2;

                    S=S+D[j-1]*v[j-1][0]; 

                }

            }
            Sol[i-1][0]=S;

            S=0;
            count++;
            if((count-1)%local_n==0)
            {
                x1++;
            }
        }
        return Sol;    
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
}



void Initialize(int **A, int *index, int local_m, int local_n, int m, int n, int my_row, int my_col, int my_rank)
{
  int *b1,*b2,*b3,*b4;
  b1 = (int*)malloc(local_m * sizeof(int));
  b3 = (int*)malloc(local_m * sizeof(int));

  b2 =(int*)malloc(local_n * sizeof(int));
  b4 = (int*)malloc(local_n * sizeof(int));

  int **M;
  M=Init_int(local_m,local_n);

  int C,D=1;
  
  MPI_Status status;
  C = 1+(my_col*(local_n))+(my_row*local_m*N);

  //printf("\np = %d\n",my_rank);
  for(int i=0; i<local_m; i++)
  {
    for(int j=0; j<local_n;j++)
    {
        M[i][j]=C;
        index[i*local_n+j]=M[i][j];
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






 /****************   MAIN   *****************/
int main(int argc, char* argv[])
{
  int local_m, local_n,C=0;
  int my_rank;
  double starttime, endtime;

  MPI_Init(&argc, &argv);
  struct GRID_INFO_T *g;
  g =(struct GRID_INFO_T*) malloc(sizeof(struct GRID_INFO_T));
  MPI_Comm_size (MPI_COMM_WORLD, &g->p );
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  test_2D_grid(g);

  local_n=25;
  local_m = (N*N)/(g->p*local_n);

  int **A;
  A= Init_int(4,(local_m)*(local_n));
  int *index;
  index = (int*)malloc(local_m*local_n* sizeof(int));
  MPI_Barrier(MPI_COMM_WORLD);
  starttime = MPI_Wtime();
  
  Initialize(A,index,local_m,local_n,g->q,g->c,g->my_row,g->my_col,g->my_rank);
  MPI_Barrier(MPI_COMM_WORLD);


  double** Dirac_mat, **v, **total;
  Dirac_mat=Init_double(local_m*local_n,1);
  total=Init_double(local_m*local_n,N*N);
  v = Init_double(N*N,1);

  int l=0;
  for(int i=0;i<N*N;i++)
  {
    for(int j=0;j<N*N;j++)
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
    Dirac_mat= Dirac(v,A,index,local_m,local_n,g->my_row,g->my_col,g->q,g->c);
    for(int k=0;k<local_m*local_n;k++)
    {
        
        total[k][l]=Dirac_mat[k][0];
        //printf("V = %f", sol.M[][0]);
    }
    l++; 

  }
    endtime = MPI_Wtime();
    /*printf("\nMatrix Dv for p = %d\n", g->my_rank);
    MPI_Barrier(MPI_COMM_WORLD);

    for(int k=0;k<local_m*local_n;k++)
    {
        for (int l = 0; l < N*N; l++)
                printf("%2.2f   ",total[k][l]);
            printf("\n");
    }*/

    if (my_rank == 0){
    printf("\nThe time took for n=100 with %d process is %f seconds\n",g->p,endtime-starttime);

  }
    

  MPI_Finalize();
  return 0;
}
