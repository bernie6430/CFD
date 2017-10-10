#include<stdio.h>
#include<stdlib.h>
#include<math.h>

const double l = 1.0;
const double gam= 0.1;
const double rho  = 1.0;
const double finA = 1.0;
const double finB = 0.0;

typedef struct tmp_info{
    double ae;
    double aw;
    double ap;
    double aww;
}info;

void TridiagonalSolve(const double *, const double *, double *, double *, double *, unsigned int);
void Exact_Sol(double *,double,int);

int main(){
    int N;
    double u;

    FILE *output = fopen("Output_Hw5.txt","w");
    FILE *output2 = fopen("Output_Hw5_Exact.txt","w");
    FILE *err = fopen("Output_err.txt","w");
    if(output==0 || output2==0){
        printf("error message!!\n");
        exit(1);
    }

    //Define the node number
    printf("enter node number:\n");
    scanf("%d",&N);
    printf("enter velocity:\n");
    scanf("%lf",&u);

    double F=rho*u;
    double D=(gam*N)/l;
    double Pe=F/D;
    //construct array for storage
    double *aw =(double*)malloc(sizeof(double)*N); //superdiagonal elements
    double *ap =(double*)malloc(sizeof(double)*N); //diagonal elements
    double *ae =(double*)malloc(sizeof(double)*N); //subdiagonal elements
    double *su =(double*)malloc(sizeof(double)*N); //rhs elements
    double *fin=(double*)malloc(sizeof(double)*N); //solution elements
    double *exact = (double*)malloc(sizeof(double)*N); //exact solution

    int i,j;
    for(i=0;i<N;i++)
        fin[i]=1.0;

    //iterative TDMA
    for(j=0;j<500;j++){
        for(i=0;i<N;i++){
            //First node
            if(i==0){
                aw[i]=0.0;
                ap[i]=4.0*D+F;
                ae[i]=-(4.0/3.0)*D;
                su[i]=((8.0/3.0)*D+1.25*F)*finA+0.125*F*(fin[i]-3.0*fin[i+1]);
            }
            //second node
            else if(i==1){
                aw[i]=-(D+F);
                ap[i]=2.0*D+F;
                ae[i]=-D;
                su[i]=-0.25*F*finA+0.125*F*(5.0*fin[i]-3.0*fin[i+1]);
            }
            //Last Node
            else if(i==(N-1)){
                aw[i]=-((4.0/3.0)*D+F);
                ap[i]=4.0*D;
                ae[i]=0.0;
                //printf("%f\n",tem[i-2]);
                su[i]=((8.0/3.0)*D-F)*finB+0.125*F*(3.0*fin[i]-2.0*fin[i-1]-fin[i-2]);
            }
            //Other Nodes
            else{
                aw[i]=-(D+F);
                ap[i]=2.0*D+F;
                ae[i]=-D;
                su[i]=0.125*F*(5.0*fin[i]-fin[i-1]-fin[i-2]-3.0*fin[i+1]);
            }
        }

        TridiagonalSolve(aw,ap,ae,su,fin,N);

    }



    Exact_Sol(exact,u,N);
    printf("\n");

    printf("Peclet Number: %f\n",Pe);

    printf("numerical solution :            analytical solution:\n");
    //Save data to txt file
    for(i=0;i<N;i++){
        printf("node %d: %f %25f\n",i+1,fin[i],exact[i]);
        fprintf(output,"%f\n",fin[i]);
        fprintf(output2,"%f\n",exact[i]);
        fprintf(err,"%f\n",fabs((fin[i]-exact[i])/exact[i]));
    }

    //release memory
    free(ae);
    free(ap);
    free(aw);
    free(su);
    free(fin);
    free(exact);

    fclose(output);
    fclose(output2);
    fclose(err);
    return 0;
}

//exact solution
void Exact_Sol(double* exact,double u,int N){
    int i;
    double x;
    for(i=0;i<N;i++){
        if(i==0)
            x=l/(2*N);
        else
            x+=l/N;
        exact[i]=((exp(rho*u*x/gam)-1)/(exp(rho*u*l/gam)-1))*(finB-finA)+finA;
    }
}
//Thomos method
void TridiagonalSolve(const double *a, const double *b, double *c, double *d, double *x, unsigned int n){
    int i;

    c[0] = c[0]/b[0];
    d[0] = d[0]/b[0];
    double id;
    for(i=1;i!=n;i++){
        id=1.0/(b[i]-c[i-1]*a[i]);
        c[i]=c[i]*id;
        d[i]=(d[i] - a[i]*d[i - 1])*id;
    }

    x[n-1]=d[n-1];
    for(i=(n-2);i!=-1;i--)
        x[i]=d[i]-c[i]*x[i+1];
}
