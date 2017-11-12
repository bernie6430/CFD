#include<stdio.h>
#include<stdlib.h>
#include<math.h>

const double l = 1.0;
const double gam= 1.0;
const double rho  = 1.0;

void TridiagonalSolve(const double *, const double *, double *, double *, double *, unsigned int);


int main(){
    int N;
    int ca;
    int step;
    double delt;

    FILE *output = fopen("Output_project3.txt","w");
    FILE *err = fopen("Output_err.txt","w");
    if(output==0){
        printf("error message!!\n");
        exit(1);
    }

    //Define the node number
    N=20;
    printf("enter the case:");
    scanf("%d",&ca);
    printf("enter the step:");
    scanf("%d",&step);
    printf("enter the delt:");
    scanf("%lf",&delt);

    double D=(gam*N)/l;
    double ap0=(rho*l)/(N*delt);
    //construct array for storage
    double *aw =(double*)malloc(sizeof(double)*N); //superdiagonal elements
    double *ap =(double*)malloc(sizeof(double)*N); //diagonal elements
    double *ae =(double*)malloc(sizeof(double)*N); //subdiagonal elements
    double *su =(double*)malloc(sizeof(double)*N); //rhs elements
    double *u=(double*)malloc(sizeof(double)*N); //solution elements

    double *u0 =(double*)malloc(sizeof(double)*N);

    int i,j,k,s;
    double delx=l/(2.0*N);

    for(i=0;i<N;i++){
        u[i]=0.0;
        if(delx<0.5)
            u0[i]=2*delx;
        else
            u0[i]=2-2*delx;
        delx+=l/N;
    }


    //explicit
    if(ca==1)
        for(k=0;k<step;k++){
            for(i=0;i<N;i++){
                //First node
                if(i==0){
                    aw[i]=0.0;
                    ap[i]=ap0-3.0*D;
                    ae[i]=D;
                    su[i]=0.0;
                    u[i]=(ap[i]*u0[i]+ae[i]*u0[i+1])/ap0;
                }
                //Last Node
                else if(i==(N-1)){
                    aw[i]=D;
                    ap[i]=ap0-3.0*D;
                    ae[i]=0.0;
                    su[i]=0.0;
                    u[i]=(ap[i]*u0[i]+aw[i]*u0[i-1])/ap0;
                }
                //Other Nodes
                else{
                    aw[i]=D;
                    ap[i]=ap0-2.0*D;
                    ae[i]=D;
                    su[i]=0.0;
                    u[i]=(ap[i]*u0[i]+aw[i]*u0[i-1]+ae[i]*u0[i+1])/ap0;
                }
            }
            for(s=0;s<N;s++)
                u0[s]=u[s];
        }
    //implicit
    else if(ca==2){
        for(k=0;k<step;k++){
            for(i=0;i<N;i++){
                //First node
                if(i==0){
                    aw[i]=0.0;
                    ap[i]=(ap0+3.0*D);
                    ae[i]=-D;
                    su[i]=ap0*u0[i];
                }
                //Last Node
                else if(i==(N-1)){
                    aw[i]=-D;
                    ap[i]=(ap0+3.0*D);
                    ae[i]=0.0;
                    su[i]=ap0*u0[i];
                }
                //Other Nodes
                else{
                    aw[i]=-D;
                    ap[i]=(ap0+2.0*D);
                    ae[i]=-D;
                    su[i]=ap0*u0[i];
                }
            }
            TridiagonalSolve(aw,ap,ae,su,u,N);
            for(s=0;s<N;s++)
                u0[s]=u[s];
        }
    }
    //CN
    else if(ca==3){
        for(k=0;k<step;k++){
            for(i=0;i<N;i++){
                //First node
                if(i==0){
                    aw[i]=0.0;
                    ap[i]=(ap0+1.5*D);
                    ae[i]=-0.5*D;
                    su[i]=(ap0-1.5*D)*u0[i]+0.5*D*u0[i+1];
                }
                //Last Node
                else if(i==(N-1)){
                    aw[i]=-0.5*D;
                    ap[i]=(ap0+1.5*D);
                    ae[i]=0.0;
                    su[i]=(ap0-1.5*D)*u0[i]+0.5*D*u0[i-1];
                }
                //Other Nodes
                else{
                    aw[i]=-0.5*D;
                    ap[i]=(ap0+D);
                    ae[i]=-0.5*D;
                    su[i]=(ap0-D)*u0[i]+0.5*D*u0[i-1]+0.5*D*u0[i+1];
                }
            }
            TridiagonalSolve(aw,ap,ae,su,u,N);
            for(s=0;s<N;s++)
                u0[s]=u[s];
        }
    }
    else{
        printf("none");
    }

    printf("\n");

    printf("numerical solution :");

    printf("\n");
    //Save data to txt file
    for(i=0;i<N;i++){
        printf("node %d: %10.9f \n",i+1,u0[i]);
        fprintf(output,"%10.9f\n",u0[i]);
    }

    //release memory
    free(ae);
    free(ap);
    free(aw);
    free(su);
    free(u);
    free(u0);
    fclose(output);
    fclose(err);
    return 0;
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
