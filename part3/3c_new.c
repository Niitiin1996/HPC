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

double **Dirac_transpose(double **v, int **A, int *index, int local_m, int local_n, int my_row, int my_col, int m, int n)
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
 
                else if(j==(A[2][i-1]))
                {
                    if((A[2][i-1])%N==0)
                    {
                        del2y=-1;
                    }
                    else
                        del2y=1;
                    flag=1;


                }
                else if(j==A[3][i-1])
                {
                    if(A[3][i-1]>N*(N-1))
                    {
                        del1y=-1;

                    }
                    else
                        del1y=1;
                    flag=1;
                }
                else if(j==A[0][i-1])
                {
                    if((A[0][i-1]-1)%N==0)
                        del_2y=-1;
                    else
                        del_2y=1;
                    flag=1;
                }
                else if(j==A[1][i-1])
                {
                    if(A[1][i-1]>N*(N-1))
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
                count++;
                if((count-1)%N==0)
                {
                    x1++;
                }
            }
            Sol[i-1][0]=S;
            S=0;

            
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

double **orient(double array[N*N][1], int *index)
{
    double **new;
    new= Init_double(N*N,1);

    for(int i=0; i<N*N; i++)
    {
        for(int j=0;j<N*N;j++)
        {
            if(index[j]==i+1)
                new[i][0]=array[j][0];
                
        }
    }

    return new;

}

double dotProduct(double A[N*N][1], double B[N*N][1]) 
{
    double dotProduct = 0;
    for(int i = 0; i < N*N; i++) 
    {
        dotProduct += A[i][0] * B[i][0];
    }
    return dotProduct;
}

double **CG(double **v, int **A, int *index, int local_m, int local_n, int my_row, int my_col, int m, int n, int my_rank)
{
    double **b,**x,**Dx,Dx_obs[N*N][1],**Dx_obs_trans, Dc[local_m*local_n][1],**r_new;
    double **DDx,r[N*N][1];
    int *index_g,*in;
    in =(int*)malloc(local_m*local_n *sizeof(int));

    
    //INDEX GLOBAL

    if(my_rank==0)
    {   
        index_g=(int*)malloc(N*N* sizeof(int));
    }
    MPI_Gather(index,local_m*local_n,MPI_INT,index_g,local_m*local_n,MPI_INT,0,MPI_COMM_WORLD);
    if(my_rank==0)
    {
        for(int i=0;i<N*N;i++)
        printf("Index: %d ",index_g[i]);}

    b=Init_double(local_m*local_n,1);
    
    Dx=Init_double(local_m*local_n,1);
    r_new=Init_double(local_m*local_n,1);
    DDx=Init_double(local_m*local_n,1);
    x= Init_double(N*N,1);


    b = Dirac_transpose(v,A,index,local_m,local_n,my_row,my_col,m,n);
    for(int i=0;i<N*N;i++)
    {
        x[i][0] = 1;
    }
    Dx=Dirac(x,A,index,local_m,local_n,my_row,my_col,m,n);

    for(int i=0;i<local_m*local_n;i++)
    {
        Dc[i][0]=Dx[i][0];
    }
    Dx_obs_trans= Init_double(N*N,1);
    
    MPI_Gather(Dc,local_m*local_n*1,MPI_DOUBLE,Dx_obs,local_m*local_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(my_rank==0)
    {
        Dx_obs_trans=orient(Dx_obs,index_g);
        for(int i=0;i<N*N;i++)
            {
                //printf("Index: %f ",Dx_obs_trans[i][0]);
            }

    }
    MPI_Bcast(*Dx_obs_trans,N*N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    DDx=Dirac_transpose(Dx_obs_trans,A,index,local_m,local_n,my_row,my_col,m,n);
    double norm1=0, norm_f1;
    int flag;
    for(int i=0; i<local_m*local_n;i++)
    {
        r[i][0]= b[i][0] - DDx[i][0];
        norm1+=r[i][0]*r[i][0];
        
    }
    MPI_Reduce(&norm1, &norm_f1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(my_rank==0)
    {
        norm_f1=sqrt(norm_f1);
        if(norm_f1<0.001)
            flag=1;
        else
            flag=0;
            

    }
    MPI_Bcast(&flag,1,MPI_INT,0,MPI_COMM_WORLD);


    double p[local_m*local_n][1],rt[N*N][1],**pt,r_new_t[N*N][1], rs[N*N][1],ps[N*N][1];

    pt=Init_double(N*N,1);
    double alpha, norm2, beta;

    int count =0;
    if(flag==0)
    {
        for(int i=0;i<local_m*local_n;i++)
            {
                p[i][0]=r[i][0];
            }

        do{
            MPI_Gather(p,local_m*local_n*1,MPI_DOUBLE,ps,local_m*local_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Gather(r,local_m*local_n*1,MPI_DOUBLE,rs,local_m*local_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
            if(my_rank==0)
            for(int i=0;i<N*N;i++)
            {
                pt[i][0]=ps[i][0];
            }

            Dx=Dirac(pt,A,index,local_m,local_n,my_row,my_col,m,n);
            for(int j=0;j<local_m*local_n;j++)
            {
                Dc[j][0]=Dx[j][0];

            }
            MPI_Gather(Dc,local_m*local_n*1,MPI_DOUBLE,Dx_obs,local_m*local_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
            if(my_rank==0)
            {
                Dx_obs_trans=orient(Dx_obs,index_g);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(*Dx_obs_trans,N*N,MPI_DOUBLE,0,MPI_COMM_WORLD);
            DDx=Dirac_transpose(Dx_obs_trans,A,index,local_m,local_n,my_row,my_col,m,n);
            for(int i=0;i<local_m*local_n;i++)
            {
                Dc[i][0]=DDx[i][0];

            } 
            MPI_Gather(Dc,local_m*local_n*1,MPI_DOUBLE,Dx_obs,local_m*local_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
            if(my_rank==0)
            {
                alpha = dotProduct(rs,rs)/dotProduct(ps,Dx_obs);
                //printf("%f",alpha);
            }
            MPI_Bcast(&alpha,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            for(int i=0;i<local_m*local_n;i++)
            {
                p[i][0] = alpha*p[i][0];
            }
            for(int i=0;i<local_m*local_n;i++)
            {
            x[i][0] = x[i][0] + p[i][0];
            }
            for(int i=0;i<local_m*local_n;i++)
            {
            r_new[i][0] = r[i][0] - alpha*DDx[i][0];
                //printf("%f   ",r_new[i][0]);
            }
            
            double norm=0;
            
            for(int i=0;i<local_m*local_n;i++)
            {
            norm = norm + (r_new[i][0]*r_new[i][0]);
            }
            printf("%f   ",norm);
            MPI_Reduce(&norm, &norm2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            
            if(my_rank==0)
            {
                //printf(  "norm  %f", norm);
                norm2=sqrt(norm2);
                if(norm2 < 0.001)
                flag = 1;
                else
                flag = 0;
            }
            MPI_Bcast(&flag,1,MPI_INT,0,MPI_COMM_WORLD);
            if(flag ==1)
            {
                break;
            }
            for(int j=0;j<local_m*local_n;j++)
            {
                Dc[j][0]=r_new[j][0];

            }
            
            MPI_Gather(Dc,local_m*local_n*1,MPI_DOUBLE,r_new_t,local_m*local_n,MPI_DOUBLE,0,MPI_COMM_WORLD);
            if(my_rank==0)
            {
                beta = dotProduct(r_new_t,r_new_t)/dotProduct(rt,rt);

            }
            MPI_Bcast(&beta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            for(int i=0;i<local_m*local_n;i++)
            {
                p[i][0] = r_new[i][0] + beta*p[i][0];
            }
            for(int i=0;i<local_m*local_n;i++)
            {
                r[i][0] = r_new[i][0];
            }
        count++;
        }while(count<2);
        //printf("COUNTTTTTTT =%d", count);

    }











    return x;

}




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
  int *index;
  index = (int*)malloc(local_m*local_n* sizeof(int));

  
  Initialize(A,index,local_m,local_n,g->q,g->c,g->my_row,g->my_col,g->my_rank);
  MPI_Barrier(MPI_COMM_WORLD);


  double** Dirac_mat, **v, **total;
  Dirac_mat=Init_double(local_m*local_n,1);
  total=Init_double(local_m*local_n,N*N);
  v = Init_double(N*N,1);

  int l=0;
  for(int i=0;i<1;i++)
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
    
    Dirac_mat= Dirac_transpose(v,A,index,local_m,local_n,g->my_row,g->my_col,g->q,g->c);
    MPI_Barrier(MPI_COMM_WORLD);
    for(int k=0;k<local_m*local_n;k++)
    {
        
        total[k][l]=Dirac_mat[k][0];
        //printf("V = %f", sol.M[][0]);
    }
    l++; 

  }

    /*printf("\nMatrix Dv for p = %d\n", g->my_rank);
    MPI_Barrier(MPI_COMM_WORLD);

    for(int k=0;k<local_m*local_n;k++)
    {
        for (int l = 0; l < N*N; l++)
                printf("%2.2f   ",total[k][l]);
            printf("\n");
    }*/

      
    double **x;
    x = CG(v,A,index,local_m,local_n,g->my_row,g->my_col,g->q,g->c,my_rank);
    /*int *x;
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank==0)
    {   
        x=(int*)malloc(local_m*local_n*g->p* sizeof(int));
    }
    MPI_Gather(index,4,MPI_INT,x,4,MPI_INT,0,MPI_COMM_WORLD);
    if(g->my_rank==0)
    {
        for(int i=0;i<N*N;i++)
        printf("Index: %d ",x[i]);
    }*/

  MPI_Finalize();
  return 0;
}
