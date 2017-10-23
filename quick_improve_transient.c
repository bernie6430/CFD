#include<stdio.h>
#include<stdlib.h>
#include<math.h>

const double l = 1.5;
const double gam= 0.03;
const double rho  = 1.0;
const double finA = 0.0;
const double finB = 0.0;

void TridiagonalSolve(const double *, const double *, double *, double *, double *, unsigned int);
void Exact_Sol(double *,double,int);

int main(){
    int N;
    double u;

    FILE *output = fopen("Output_Hwt.txt","w");
    FILE *output2 = fopen("Output_Hwt_Exact.txt","w");
    FILE *err = fopen("Output_err.txt","w");
    if(output==0 || output2==0){
        printf("error message!!\n");
        exit(1);
    }

    //Define the node number
    N=45;
    u=2.0;

    double delt=0.01;
    double F=rho*u;
    double D=(gam*N)/l;
    double ap0=(rho*l)/(N*delt);
    double Pe=F/D;
    //construct array for storage
    double *aw =(double*)malloc(sizeof(double)*N); //superdiagonal elements
    double *ap =(double*)malloc(sizeof(double)*N); //diagonal elements
    double *ae =(double*)malloc(sizeof(double)*N); //subdiagonal elements
    double *su =(double*)malloc(sizeof(double)*N); //rhs elements
    double *fin=(double*)malloc(sizeof(double)*N); //solution elements
    double *exact = (double*)malloc(sizeof(double)*N); //exact solution

    double *fin0 =(double*)malloc(sizeof(double)*N);
    double *source =(double*)malloc(sizeof(double)*N);

    int i,j,k,s;
    double delx=l/(2.0*N);

    for(i=0;i<N;i++){
        fin0[i]=0.0;
        fin[i]=0.0;
        if(delx<0.6){
            source[i]=-200.0*(delx)+100.0;
            delx+=l/N;
        }
        else if(delx>0.6 && delx<0.8){
            source[i]=100.0*(delx)-80.0;
            delx+=l/N;
        }
        else{
            source[i]=0.0;
            delx+=l/N;
        }
    }


    //iterative TDMA
    for(k=0;k<180;k++){
        for(s=0;s<N;s++)
            fin[s]=0.0;
    for(j=0;j<80;j++){
        for(i=0;i<N;i++){
            //First node
            if(i==0){
                aw[i]=0.0;
                ap[i]=4.0*D+F+ap0;
                ae[i]=-(4.0/3.0)*D;
                su[i]=0.125*F*(fin[i]-3.0*fin[i+1])+(l/N)*source[i]+ap0*fin0[i];
            }
            //second node
            else if(i==1){
                aw[i]=-(D+F);
                ap[i]=2.0*D+F+ap0;
                ae[i]=-D;
                su[i]=0.125*F*(5.0*fin[i]-3.0*fin[i+1])+(l/N)*source[i]+ap0*fin0[i];
            }
            //Last Node
            else if(i==(N-1)){
                aw[i]=-(D+F);
                ap[i]=D+F+ap0;
                ae[i]=0.0;
                su[i]=0.125*F*(3.0*fin[i]-2.0*fin[i-1]-fin[i-2])+(l/N)*source[i]+ap0*fin0[i];
            }
            //Other Nodes
            else{
                aw[i]=-(D+F);
                ap[i]=2.0*D+F+ap0;
                ae[i]=-D;
                su[i]=0.125*F*(5.0*fin[i]-fin[i-1]-fin[i-2]-3.0*fin[i+1])+(l/N)*source[i]+ap0*fin0[i];
            }
        }
        TridiagonalSolve(aw,ap,ae,su,fin,N);
    }
    for(s=0;s<N;s++){
        fin0[s]=fin[s];
    }
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
    free(fin0);
    free(source);
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
