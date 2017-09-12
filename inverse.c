/*=================================================
Inverse matrix calculation by LU decomposition
Author:Bernie Huang
Date:9/11/2017
Version:1.0
===================================================*/
#include<stdio.h>
#include<stdlib.h>

void LU_comp(double* A, double* L, double* U,int n);
void inverse(double* L, double* U, double* AI,int n);

int main(){
    int n,i,j,k,d;
    double *A,*L,*U,*AI;
    printf("please enter the dimension of matrix:");
    scanf("%d",&n);

    d=n*n;
    A=(double*)malloc(d*sizeof(double));
    L=(double*)malloc(d*sizeof(double));
    U=(double*)malloc(d*sizeof(double));
    AI=(double*)malloc(d*sizeof(double));

    int ct=0;
    for(i=0;i<d;i++){
        A[i]=0.0;
        L[i]=0.0;
        U[i]=0.0;
        if(i==ct*(n+1)){
           AI[i]=1.0;
           ct++;
        }
        else
        {
           AI[i]=0.0;
        }
    }

    printf("please enter the A elements:\n");

    for(i=0;i<d;i++){
        scanf("%lf",&A[i]);
    }

    LU_comp(A,L,U,n);
    inverse(L,U,AI,n);

    free(A);
    free(L);
    free(U);
    free(AI);
    return 0;
}

void LU_comp(double* A, double* L, double* U,int n){
    int i,j,k,d;
    double coeff;
    d=n*n;
    //Lower matrix
    for(i=0;i<n;i++){
        L[i*(n+1)]=1.0;
    }

    for(i=0;i<n;i++){
        for(j=i+1;j<n;j++){
            coeff=A[j*n+i]/A[i*n+i];
            L[j*n+i]=coeff;
            for(k=i+1;k<n;k++)
                A[j*n+k]=A[j*n+k]-coeff*A[i*n+k];
        }
    }

    //Upper matrix
    for(i=0;i<n;i++){
        for(j=i;j<n;j++)
            U[i*n+j]=A[i*n+j];
    }
}

void inverse(double* L, double* U, double* AI,int n){
    int i,j,k;
    double *X,*Y;
    X=(double*)malloc(n*sizeof(double));
    Y=(double*)malloc(n*sizeof(double));
    for(i=0;i<n;i++) X[i]=0.0;
    for(i=0;i<n;i++) Y[i]=0.0;

    for(k=0;k<n;k++){
        //forward substitution
        Y[0]=AI[k*n+0];
        for(i=1;i<n;i++){
            Y[i]=AI[k*n+i];
            for(j=0;j<i;j++)
                Y[i]=Y[i]-L[i*n+j]*Y[j];
        }

        //Backward substitution
        X[n-1]=Y[n-1]/U[n*n-1];
        for(i=n-2;i>-1;i--){
            X[i]=Y[i];
            for(j=n-1;j>i;j--)
                X[i]=X[i]-U[i*n+j]*X[j];
            X[i]=X[i]/U[i*(n+1)];
        }
        printf("%d row:\n",k+1);
        for(i=0;i<n;i++)
            printf("%f \n",X[i]);
        printf("\n");
    }

    free(X);
    free(Y);
}
