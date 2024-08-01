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

struct vector dirac (int v[64][1], int A[4][64],int N )
{   int i,j,x,y,x1=1,flag=0,count=1;
    float del,del1y,del_1y,del2y,del_2y;
    float m =0.4;
    double S=0, D; 
    struct vector Sol;
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
                
                D=m*del+(U1(i-1)*del1y-U1(abs(A[3][i-1]))*del_1y)/2+pow((-1),x1)*(U2(i-1)*del2y-U2(abs(A[2][i-1]))*del_2y)/2;
                S=S+D*v[j-1][0];                
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

int main() 
{
   int i,j=0, N=8, v[64][1];
   int C=1,D=1;
   int A[4][64];
   int k,l;
   struct vector sol;
   FILE *ptr;

   double Dirac[64][64];
// BORDER POINTS
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

   //Border Arrays
   printf("myArray =");
   for (k = 0; k < 4; k++) 
    {
        printf("\n");
        for (l = 0; l < N*N; l++) 
        {
            printf("%d ",A[k][l]);
        }
    }
    l=0;
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
        sol= dirac(v, A, N);
        
        for(k=0;k<N*N;k++)
        {
        
        Dirac[k][l]=sol.M[k][0];
        //printf("V = %f", sol.M[][0]);
        }
        l++;              
    }

    ptr= fopen("Matrix.txt","w");


    printf("Matrix Dv =");
    for (int k = 0; k < N*N; k++) 
    {
        

        for (int l = 0; l < N*N; l++) 
            fprintf(ptr,"%2.2f  ",Dirac[k][l]);
        fprintf(ptr,"\n");

    }
    fclose(ptr);
    return 0;
}