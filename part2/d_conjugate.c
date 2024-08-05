#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct vector
{
    double M[64][1];

};
int U1(int x)
{
    float i=1.0;
    return i;
}

int U2(int x)
{
    float i=1.0;
    return i;
}



struct vector dirac_transpose(int A[4][64],int N,double Sol[64][1] )
{
    int i, j, x, y, x1=1, flag=0, ini=0;
    float del,del1y,del_1y,del2y,del_2y;
    float m =0.4;
    double S=0, D[64];
    struct vector D_t;
    int count = 1;
    for(i=0;i<N*N;i++)
    {
        D_t.M[i][0]=0;
    }

    for(i=1 ; i<=N*N; i++)
    {
        for(j=1; j<=N*N; j++)
        {
            flag=0;
            del=0; 
            del1y=0;
            del_1y=0;
            del2y=0;
            del_2y=0;
            if(i==j)
            {
                del=1;
                flag=1;    
            }
            else if(j==abs(A[0][i-1]))
            {
                if(A[0][i-1]<0)
                    del2y=-1;
                else
                    del2y=1;
                flag=1;
            }
            else if(j==abs(A[1][i-1]))
            {
                if(A[1][i-1]<0)
                    del1y=-1;
                else
                    del1y=1;
                flag=1;
            }
            else if(j==abs(A[2][i-1]))
            {
                if(A[2][i-1]<0)
                    del_2y=-1;
                else
                    del_2y=1;
                flag=1;
            }
            else if(j==abs(A[3][i-1]))
            {
                if(A[3][i-1]<0)
                    del_1y=-1;
                else
                    del_1y=1;
                flag=1;
            }
            if(flag==1)
            {
                D[j-1]=m*del+(U1(i-1)*del1y-U1(abs(A[3][i-1]))*del_1y)/2+pow((-1),x1)*(U2(i-1)*del2y-U2(abs(A[2][i-1]))*del_2y)/2;
                

              //  D_t.M[j-1][0]=D_t.M[j-1][0]+D[j-1]*Sol[j-1][0];
              D_t.M[j-1][0]=D_t.M[j-1][0]+D[j-1]*Sol[i-1][0];

              
            }           
        }
        count++;
        if((count-1)%N==0)
        {
            x1++;
        }

    }
    //printf("\n");
    return D_t;
}

struct vector dirac (double v[64][1], int A[4][64],int N )
{   int i, j, x, y, x1=1, flag=0, ini=0;
    float del,del1y,del_1y,del2y,del_2y;
    float m =0.4;
    double S=0, D[64], D_trans[64][1]; 
    struct vector Sol;
    int count = 1;
    for(i=1 ; i<=N*N; i++)
    {
        for(j=1; j<=N*N; j++)
        {
            flag=0;
            del=0; 
            del1y=0;
            del_1y=0;
            del2y=0;
            del_2y=0;

            if(i==j)
            {
                del=1;
                flag=1;    
            }

            else if(j==abs(A[0][i-1]))
            {
                if(A[0][i-1]<0)
                    del2y=-1;
                else
                    del2y=1;
                flag=1;
            }

            else if(j==abs(A[1][i-1]))
            {
                if(A[1][i-1]<0)
                    del1y=-1;
                else
                    del1y=1;
                flag=1;
            }

            else if(j==abs(A[2][i-1]))
            {
                if(A[2][i-1]<0)
                    del_2y=-1;
                else
                    del_2y=1;
                flag=1;
            }

            else if(j==abs(A[3][i-1]))
            {
                if(A[3][i-1]<0)
                    del_1y=-1;
                else
                    del_1y=1;
                flag=1;
            }

            if(flag==1)
            {
                D[j-1]=m*del+(U1(i-1)*del1y-U1(abs(A[3][i-1]))*del_1y)/2+pow((-1),x1)*(U2(i-1)*del2y-U2(abs(A[2][i-1]))*del_2y)/2;
                
                S=S+D[j-1]*v[j-1][0];               
            }
            
        }
        Sol.M[i-1][0]=S;
        S=0;
        count++;
        if((count-1)%N==0)
        {
            x1++;
        }  
    }
    return Sol;
}

double dotProduct(double A[64][1], double B[64][1], int N) 
{
    double dotProduct = 0;
    for(int i = 0; i < N*N; i++) 
    {
        dotProduct += A[i][0] * B[i][0];
    }
    return dotProduct;
}

struct vector conjugate_gradient(int A[4][64], int N, double e[64][1])
{
    double alpha, r[N*N][1], p[N*N][1], AAp[N*N][1], r_new[N*N][1];
    struct vector b , x, Dx, DDx, Dp, DDp, x_new;
    double norm1=0,norm2=0, beta ;
    int k=0;

    b = dirac_transpose(A,N,e);
    
    for(int i=0;i<N*N;i++)
    {
        x.M[i][0] = 0;
    }

    Dx = dirac(x.M, A, N);
    
    DDx = dirac_transpose(A,N,Dx.M);

    for(int i=0;i<N*N;i++)
    {
        r[i][0] = b.M[i][0] - DDx.M[i][0];
        //printf("r=  %f ",r[i][0]);
    }

    for(int i=0;i<N*N;i++)
    {
        norm1 = norm1 + (r[i][0]*r[i][0]); 
            
    }
    norm1=sqrt(norm1); 
    //printf(" Norm1 =%f", norm1);

    if(norm1>0.000001)
    {

        for(int i=0;i<N*N;i++)
        {
         p[i][0] = r[i][0];
        }
        do
        {
            norm2=0;
            Dp = dirac(p,A,N);
            DDp = dirac_transpose(A,N,Dp.M);

            alpha = dotProduct(r,r,N)/dotProduct(p,DDp.M,N);
           // printf("\nalpha = %f",alpha);

            for(int i=0;i<N*N;i++)
            {
                p[i][0] = alpha*p[i][0];
            }
            /*for(int i=0;i<N*N;i++)
            {
                printf(" P = %f",p[i][0]);
            }*/

            for(int i=0;i<N*N;i++)
            {
            x.M[i][0] = x.M[i][0] + p[i][0];
            }
            /*for(int i=0;i<N*N;i++)
            {
                printf(" X = %f",x.M[i][0]);
            }*/

            for(int i=0;i<N*N;i++)
            {
            r_new[i][0] = r[i][0] - alpha*DDp.M[i][0];
            }
            /*for(int i=0;i<N*N;i++)
            {
                printf(" R_new = %f",r_new[i][0]);
            }*/

            for(int i=0;i<N*N;i++)
            {
            norm2 = norm2 + (r_new[i][0]*r_new[i][0]);
            
            }
            norm2=sqrt(norm2);
            //printf(" Norm2 = %f",norm2);

            if(norm2 < 0.001)
            {
                break;
            }
            beta = dotProduct(r_new,r_new,N)/dotProduct(r,r,N);
            //printf(" Beta = %f",beta);

            for(int i=0;i<N*N;i++)
            {
                p[i][0] = r_new[i][0] + beta*p[i][0];
            }
            for(int i=0;i<N*N;i++)
            {
                r[i][0] = r_new[i][0];
            }
            /* code */
            k++;
        } while (k<100); 
        //printf("K= %d\n",k);
        return x;    
    }

}

double Corr(int A[4][64], int N, int t)
{
    int i,j, x, x_t;
    double e[64][1];
    struct vector D_inv;
    double correlation = 0.0;


    for( i=0; i<N*N; i++)
    {
        x = i+1;
        x_t = x; 

        for(int z=0;z<t;z++)
        {
            x_t = abs(A[1][x_t-1]);
        }
        //printf("\nx_t = %d",x_t);
        

        for(j=0; j<N*N; j++)
        {
            if(j==x_t-1)
            {
                e[j][0]=1;
            }
            else
            {
                e[j][0]=0;
            }
        }

        D_inv = conjugate_gradient(A,N,e);

        correlation += pow(D_inv.M[x-1][0],2)/(N*N);
    }
    return correlation;

}

int main() 
{
   int i,j=0,N=8;
   int C=1,D=1,t=1;
   int A[4][64];
   int k,l;
   double Dirac[64][64], Dirac_inv[64][64],  v[64][1];
   struct vector sol, sol2, sol_con;
   double e[64][1];
   double correlation;
   FILE *ptr;
   
   // Initialize U as unit matrix
   for(i = 0; i < N*N; i++)
   {
    
    //1st pont
    if(C<D*N)
    {

        A[j][i]=C+1;
        j++;
    }
    else
    {
        A[j][i]=-(C-N+1);
        j++;
    }
   
    //2nd point
    if(D-1==0)
    {
        A[j][i]=-((N-1)*N+C);
        j++;
    }
    else
    {
        A[j][i]=C-N;
        j++;
    }

    //3rd point
    if((C-1)%N==0)
    {
        A[j][i]=-(C-1+N);
        j++;
        
    }
    else
    {
        A[j][i]=C-1;
        j++;
    }
    //4th point
    if(D==N)
    {
        A[j][i]=-(C-(N-1)*N);
        j++;
    }
    else
    {
        A[j][i]=C+N;
        j++;
    }
    C++;
    j=0;
    if((i+1)%N==0)
        D++;  
   }

    printf("myArray =");
    for (k = 0; k < 4; k++) 
    {
        printf("\n");
        for (l = 0; l < N*N; l++) 
        {
            printf("%d      ", A[k][l]);
        }
    }

    /*l=0;
    for (i = 0; i < N*N; i++) 
    {
        for (j = 0; j < N*N; j++)
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
    
        sol= dirac_transpose(A, N, v);
        
        for(k=0;k<N*N;k++)
        {        
            Dirac[k][l]=sol.M[k][0];
        }
        l++;              
    }    
    //sol2 = dirac_transpose(v,A,N,sol.M);
    printf("Matrix Dv = \n");
     for (int k = 0; k < N*N; k++) 
    
    {
        for(int l=0;l<64;l++)
        {            
            printf("%2.2f  ",Dirac[k][l]);
        }
        
        printf("\n");
    }*/
    l=0;
    for (i = 0; i < N*N; i++) 
    {
        for (j = 0; j < N*N; j++)
        {

            if(i==j)
            {
                e[j][0]=1;
            }
            else
            {
                e[j][0]=0;
            }
        }
    
        sol_con= conjugate_gradient(A, N, e);
        
        for(k=0;k<N*N;k++)
        {
        
            Dirac_inv[k][l]=sol_con.M[k][0];
        }
        l++;           
    }

    for(int t=1;t<N;t++)
    {
        correlation = Corr(A, N, t);
        printf("\n%f",correlation);

    }


    ptr= fopen("Matrix_Inverse.txt","w");
     for (int k = 0; k < N*N; k++) 
    
    {
        for(int l=0;l<N*N;l++)
        {
            fprintf(ptr, "%2.2f ",Dirac_inv[k][l]);
        }
        fprintf(ptr,"\n");
    }
    return 0;
}