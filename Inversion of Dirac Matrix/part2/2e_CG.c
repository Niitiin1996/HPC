#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

double *a1, *a2;
int Coun=0;

/// @brief 
double** Initialize(double**M)
{
    M = (double**)malloc(576 * sizeof(double*));
    for (int i = 0; i < 576; i++)
        M[i] = (double*)malloc(1 * sizeof(double));
    return M;

}

double complex** Initialize_v(double complex**M)
{
    M = (double complex**)malloc(576 * sizeof(double complex*));
    for (int i = 0; i < 576; i++)
        M[i] = (double complex*)malloc(1 * sizeof(double complex));
    return M;

};

void Read(int N)
{
    a1 = (double*)malloc(N*N * sizeof(double));
    a2 = (double*)malloc(N*N * sizeof(double));    
    FILE *f;
    f=fopen("U_1","r");
    if(f== NULL )
    {
        printf("File not found" );
        exit(1);

    }
    int x1,x2,mu;
    for(int i=0; i<N*N; i++)
    {
        fscanf(f, "%d%d%d%lf",&x1,&x2,&mu,&a1[i]);
        fscanf(f, "%d%d%d%lf",&x1,&x2,&mu,&a2[i]);
    }

    fclose(f);
}

double complex U_1(int x)
{
    double complex z;
    z=0 + cexp(a1[x]*I); 
    return z;
}

double complex U_2(int x)
{
    double complex z;
    z=0 + cexp(a2[x]*I); 
    return z;
}

double complex** dirac_transpose(int** A,int N,double complex** Sol)
{
    int i, j, x, y, x1=1, flag=0, ini=0;
    float del,del1y,del_1y,del2y,del_2y;
    float m =0.24;
    double S=0;
    double complex *D=(double complex*)malloc(N*N * sizeof(double complex));
    double complex **D_t;
    D_t=Initialize_v(D_t);
    int count = 1;
    
    for(i=0;i<N*N;i++)
    {
        D_t[i][0]=0;
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
                {
                    del2y=-1;
                    A[0][i-1]=A[0][i-1]*-1;
                }
                else
                    del2y=1;
                flag=1;
            }
            else if(j==abs(A[1][i-1]))
            {
                if(A[1][i-1]<0)
                {
                    del1y=-1;
                    A[1][i-1]=A[1][i-1]*-1;
                }
                    
                else
                    del1y=1;
                flag=1;
            }
            else if(j==abs(A[2][i-1]))
            {
                if(A[2][i-1]<0)
                {
                    del_2y=-1;
                    A[2][i-1]=A[2][i-1]*-1;
                }
                else
                    del_2y=1;
                flag=1;
            }
            else if(j==abs(A[3][i-1]))
            {
                //printf("\n%d",A[3][i-1]);
                if(A[3][i-1]<0)
                {
                    del_1y=-1;
                    A[3][i-1]=A[3][i-1]*-1;
                }
                else
                    del_1y=1;
                flag=1;
            }
            if(flag==1)
            {
                //printf("\n%d",A[3][i-1]);

                D[j-1]=m*del+(U_1(i-1)*del1y-conjf(U_1(A[3][i-1]))*del_1y)/2+pow((-1),x1)*(U_2(i-1)*del2y-conjf(U_2(A[2][i-1]))*del_2y)/2;


              //  D_t.M[j-1][0]=D_t.M[j-1][0]+D[j-1]*Sol[j-1][0];
              D_t[j-1][0]=D_t[j-1][0]+D[j-1]*Sol[i-1][0];

              
            }           
        }
        count++;

        if((count-1)%N==0)
        {
            x1++;
        }

    }
    
    return D_t;
}

double complex** dirac (double complex **v, int **A,int N )
{   int i, j, x, y, x1=1, flag=0, ini=0,C=0;
    float del,del1y,del_1y,del2y,del_2y;
    float m =0.24;
    double complex S=0, *D = (double complex*)malloc(N*N * sizeof(double complex)) ; 
    double complex **Sol;
    Sol = Initialize_v(Sol);
    
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
                {
                    del2y=-1;
                    A[0][i-1]=abs(A[0][i-1]);
                }
                else
                    del2y=1;
                flag=1;
            }

            else if(j==abs(A[1][i-1]))
            {
                if(A[1][i-1]<0)
                {
                    del1y=-1;
                    A[1][i-1]=abs(A[1][i-1]);
                }
                else
                    del1y=1;
                flag=1;
            }

            else if(j==abs(A[2][i-1]))
            {
                if(A[2][i-1]<0)
                {
                    del_2y=-1;
                    A[2][i-1]=abs(A[2][i-1]);
                }
                else
                    del_2y=1;
                flag=1;
            }

            else if(j==abs(A[3][i-1]))
            {
                if(A[3][i-1]<0)
                {
                    del_1y=-1;
                    A[3][i-1]=abs(A[3][i-1]);
                }
                else
                    del_1y=1;
                flag=1;
            }

            if(flag==1)
            {
                D[j-1]=m*del+(U_1(i-1)*del1y-conjf(U_1(A[3][i-1]))*del_1y)/2+pow((-1),x1)*(U_2(i-1)*del2y-conjf(U_2(A[2][i-1]))*del_2y)/2;
                S=S+D[j-1]*v[j-1][0];               
            }
            
        }
        Sol[i-1][0]=S;
        S=0;
        count++;
        if((count-1)%N==0)
        {
            x1++;
        }  
    }
    
    free(D);
    return Sol;
}

double complex dotProduct(double complex **A, double complex **B, int N) 
{
    double complex dotProduct = 0;
    for(int i = 0; i < N*N; i++) 
    {
        dotProduct += A[i][0] * B[i][0];
    }
    return dotProduct;
}

double complex ** conjugate_gradient(int **A, int N, double complex **e)
{
    int i,C=0;
    
    double complex alpha;
    double complex **r, **p, **AAp, **r_new;
    r = Initialize_v(r);

    p = Initialize_v(p);
    AAp = Initialize_v(AAp);
    r_new = Initialize_v(r_new);
    
    double complex**b , **x, **Dx, **DDx, **Dp, **DDp, **x_new;
    b=Initialize_v(b);
    x=Initialize_v(x);
    Dx=Initialize_v(Dx);
    DDx=Initialize_v(DDx);
    Dp=Initialize_v(Dp);
    DDp=Initialize_v(DDp);
    x_new=Initialize_v(x_new);
    
    double norm1=0,norm2=0, beta ;
    int k=0;

    b = dirac_transpose(A,N,e);

    
    for(int i=0;i<N*N;i++)
    {
        x[i][0] = 0;
    }


    Dx = dirac(x, A, N);
  
    DDx = dirac_transpose(A,N,Dx);
    
    for(int i=0;i<N*N;i++)
    {
        r[i][0] = b[i][0] - DDx[i][0];
        //printf("r=  %f ",r[i][0]);
    }

    for(int i=0;i<N*N;i++)
    {

        norm1 = norm1 + conjf(r[i][0]*r[i][0]);
             
    }
    
    
        

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
            DDp = dirac_transpose(A,N,Dp);
            
           

            alpha = dotProduct(r,r,N)/dotProduct(p,DDp,N);
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
            x[i][0] = x[i][0] + p[i][0];
            }
            /*for(int i=0;i<N*N;i++)
            {
                printf(" X = %f",x.M[i][0]);
            }*/

            for(int i=0;i<N*N;i++)
            {
            r_new[i][0] = r[i][0] - alpha*DDp[i][0];
            }
            /*for(int i=0;i<N*N;i++)
            {
                printf(" R_new = %f",r_new[i][0]);
            }*/

            for(int i=0;i<N*N;i++)
            {
            norm2 = norm2 + sqrt((r_new[i][0]*r_new[i][0]));
            }
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
            //code 
            k++;
        } while (k<100);
        
        
        //printf("K= %d\n",k);
        free(r); free(p); free(AAp); free(r_new);
        return x; 
           
    }
    else
    {
        return x; 

    }

}

double Corr(int **A, int N, int t)
{
    
    int i,j, x, x_t;
    double complex**e;
    e=Initialize_v(e);
    double complex **D_inv;
    D_inv=Initialize_v(D_inv);
    double correlation = 0.0;

    for( i=0; i<N*N; i++)
    {
        printf("i=%d",i);
        x = i+1;
        x_t = x; 

        for(int z=0;z<t;z++)
        {
            x_t = abs(A[1][x_t-1]);
        }
        //printf("x_t =%d\n",x_t);

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
        
        correlation += pow(cabsf(D_inv[x-1][0]),2)/(N*N);
        //printf("%f",correlation);

    }
    return correlation;

}

int main() 
{

   int i,j=0,N=24;
   int C=1,D=1,t=1;
   int** A = (int**)malloc(24 * sizeof(int*));
   double complex**v;
   v=Initialize_v(v);
   double correlation;
   int k,l;
   double complex **Dirac,**Sol1;
   Dirac = Initialize_v(Dirac);
   Sol1= (double complex**)malloc(N*N * sizeof(double complex*));
   for (int i = 0; i < N*N; i++)
        Sol1[i] = (double complex*)malloc(N*N * sizeof(double complex));
        
   //struct vector sol, sol2, sol_con;
   //double e[N*N][1];
   //double correlation;
   FILE *ptr;
   ptr= fopen("Correlation.txt","w");
   for (i = 0; i < 24; i++)
        A[i] = (int*)malloc(N*N * sizeof(int));

   Read(N);

   // Initialize Border matrix
   for(i = 0; i < N*N; i++)
   {
    
    //1st point
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
        A[j][i]=-(C+N-1);
        j++;       
    }
    else
    {
        A[j][i]=C-1;
        j++;
    }

    //24th point
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

   

    /*printf("myArray =");
    for (k = 0; k < 24; k++) 
    {
        printf("\n");
        for (l = 0; l < 1; l++) 
        {
            printf("%d      ", A[k][l]);
        }
    }*/

    //fclose(ptr);

    for (j = 0; j < N*N; j++)
    {

            if(j==0)
            {
                v[j][0]=1;
            }
            else
            {
                v[j][0]=0;
            }
            
    }

     
    //Dirac= dirac_transpose(A, N, v);

    for(int t=0;t<N;t+50)
    {
        printf("%d:-",(t+1));
        correlation = Corr(A, N, (t+1));
        printf("%2.5f\n",correlation);
    }
    /*ptr= fopen("Matrix_Inverse.txt","w");*/
    /*for (int k = 0; k < N*N; k++) 
    {
        
        printf("%2.2f%+.2fi ",Dirac[k][0]);
        printf("\n");
    }
    free(A);
    free(Dirac);
    return 0;*/

}