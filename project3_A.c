#include<stdio.h>
#include<stdlib.h>
#include<math.h>

const double l = 1.0;
const double gam= 0.1;
const double rho  = 1.0;
const double finA = 0.0;
const double finB = 0.0;

void TridiagonalSolve(const double *, const double *, double *, double *, double *, unsigned int);
void Exact_Sol(double *,double,int);

int main(){
    int N;
    int ca;
    int step;
    double delt;
    double it;  //instant time

    FILE *output = fopen("Output_project3A.txt","w");
    if(output==0){
        printf("error message!!\n");
        exit(1);
    }

    //Define the node number
    N=50;

    //define scheme and time step
    printf("select scheme 1.explicit 2.implicit 3.CN:");
    scanf("%d",&ca);
    printf("enter the time:");
    scanf("%lf",&it);
    printf("enter the delt:");
    scanf("%lf",&delt);


    double D=(gam*N)/l;
    double ap0=(rho*l)/(N*delt);
    double Fe;
    double Fw;
    double Fe0;
    double Fw0;
    //construct array for storage
    double *aw =(double*)malloc(sizeof(double)*N); //superdiagonal elements
    double *ap =(double*)malloc(sizeof(double)*N); //diagonal elements
    double *ae =(double*)malloc(sizeof(double)*N); //subdiagonal elements
    double *su =(double*)malloc(sizeof(double)*N); //rhs elements
    double *u=(double*)malloc(sizeof(double)*N); //solution elements

    double *u0 =(double*)malloc(sizeof(double)*N);

    int i,j,k,s;

    //initialize
    for(i=0;i<N;i++){
        u[i]=0.0;
        u0[i]=0.0;
    }

    step=it/delt;
    double checku;

    //explicit
    if(ca==1)
        for(k=0;k<step;k++){
            for(i=0;i<N;i++){
                //First node
                if(i==0){
                    Fw=0.5*(u0[i]+0.1*(1+sin(6.0*k*delt)));
                    Fe=0.5*(u0[i]+u0[i+1]);
                    aw[i]=0.0;
                    ap[i]=ap0-3.0*D-0.5*Fe;
                    ae[i]=D-0.5*Fe;
                    su[i]=0.1*(1.0+sin(6.0*k*delt))*(2.0*D+0.1*(1.0+sin(6.0*k*delt)));
                    u[i]=(ap[i]*u0[i]+ae[i]*u0[i+1]+su[i])/ap0;
                    //printf("%f\n",Fe);
                }
                //Last Node
                else if(i==(N-1)){
                    Fw=0.5*(u0[i]+u0[i-1]);
                    Fe=0.5*(u0[i]+1.0);
                    aw[i]=D+0.5*Fw;
                    ap[i]=ap0-3.0*D+0.5*Fw;
                    ae[i]=0.0;
                    su[i]=(2.0*D-1.0);
                    u[i]=(ap[i]*u0[i]+aw[i]*u0[i-1]+su[i])/ap0;
                    //printf("%f\n",Fw);
                }
                //Other Nodes
                else{
                    Fw=0.5*(u0[i]+u0[i-1]);
                    Fe=0.5*(u0[i]+u0[i+1]);
                    aw[i]=D+0.5*Fw;
                    ap[i]=ap0-2.0*D-0.5*(Fe-Fw);
                    ae[i]=D-0.5*Fe;
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
            while(1){
                checku=u[N-1];
                for(i=0;i<N;i++){
                    //First node
                    if(i==0){
                        Fw=0.5*(u[i]+0.1*(1+sin(6.0*(k+1)*delt)));
                        Fe=0.5*(u[i]+u[i+1]);
                        aw[i]=0.0;
                        ap[i]=ap0+3.0*D+0.5*Fe;
                        ae[i]=-(D-0.5*Fe);
                        su[i]=ap0*u0[i]+0.1*(1.0+sin(6.0*(k+1)*delt))*(2.0*D+0.1*(1.0+sin(6.0*(k+1)*delt)));
                        //printf("%f\n",Fe);
                    }
                    //Last Node
                    else if(i==(N-1)){
                        Fw=0.5*(u[i]+u[i-1]);
                        Fe=0.5*(u[i]+1.0);
                        aw[i]=-(D+0.5*Fw);
                        ap[i]=ap0+3.0*D-0.5*Fw;
                        ae[i]=0.0;
                        su[i]=ap0*u0[i]+(2.0*D-1.0);
                        //printf("%f\n",Fw);
                    }
                    //Other Nodes
                    else{
                        Fw=0.5*(u[i]+u[i-1]);
                        Fe=0.5*(u[i]+u[i+1]);
                        aw[i]=-(D+0.5*Fw);
                        ap[i]=ap0+2.0*D+0.5*(Fe-Fw);
                        ae[i]=-(D-0.5*Fe);
                        su[i]=ap0*u0[i];
                    }
                }
                TridiagonalSolve(aw,ap,ae,su,u,N);
                if((fabs(checku-u[N-1]))<0.00000001)
                   break;
            }
            for(s=0;s<N;s++)
                u0[s]=u[s];
        }
    }
    //CN
    else if(ca==3){
        for(k=0;k<step;k++){
            while(1){
                checku=u[N-1];
                for(i=0;i<N;i++){
                    //First node
                    if(i==0){
                        Fw0=0.5*(u0[i]+0.1*(1+sin(6.0*k*delt)));
                        Fe0=0.5*(u0[i]+u0[i+1]);
                        Fw=0.5*(u[i]+0.1*(1+sin(6.0*(k+1)*delt)));
                        Fe=0.5*(u[i]+u[i+1]);
                        aw[i]=0.0;
                        ap[i]=ap0+0.5*(3.0*D+0.5*Fe);
                        ae[i]=-0.5*(D-0.5*Fe);
                        su[i]=0.5*(D-0.5*Fe0)*u0[i+1]+(ap0-0.5*(3.0*D+0.5*Fe0))*u0[i]+0.5*(0.1*(1.0+sin(6.0*(k+1)*delt))*(2.0*D+0.1*(1.0+sin(6.0*(k+1)*delt))))+0.5*(0.1*(1.0+sin(6.0*k*delt))*(2.0*D+0.1*(1.0+sin(6.0*k*delt))));
                        //printf("%f\n",Fe);
                    }
                    //Last Node
                    else if(i==(N-1)){
                        Fw0=0.5*(u0[i]+u0[i-1]);
                        Fe0=0.5*(u0[i]+1.0);
                        Fw=0.5*(u[i]+u[i-1]);
                        Fe=0.5*(u[i]+1.0);
                        aw[i]=-0.5*(D+0.5*Fw);
                        ap[i]=ap0+0.5*(3.0*D-0.5*Fw);
                        ae[i]=0.0;
                        su[i]=0.5*(D+0.5*Fw0)*u0[i-1]+(ap0-0.5*(3.0*D-0.5*Fw0))*u0[i]+(2.0*D-1.0);
                        //printf("%f\n",Fw);
                    }
                    //Other Nodes
                    else{
                        Fw0=0.5*(u0[i]+u0[i-1]);
                        Fe0=0.5*(u0[i]+u0[i+1]);
                        Fw=0.5*(u[i]+u[i-1]);
                        Fe=0.5*(u[i]+u[i+1]);
                        aw[i]=-0.5*(D+0.5*Fw);
                        ap[i]=ap0+0.5*(2.0*D+0.5*(Fe-Fw));
                        ae[i]=-0.5*(D-0.5*Fe);
                        su[i]=0.5*(D-0.5*Fe0)*u0[i+1]+0.5*(D+0.5*Fw0)*u0[i-1]+(ap0-0.5*(2.0*D+0.5*Fe0-0.5*Fw0))*u0[i];
                    }
                }
                TridiagonalSolve(aw,ap,ae,su,u,N);
                if((fabs(checku-u[N-1]))<0.00000001)
                   break;
            }
            for(s=0;s<N;s++)
                u0[s]=u[s];
        }
    }
    else{
        exit(1);
    }

    printf("\n");

    printf("numerical solution :");

    printf("\n");
    //Save data to txt file

    for(i=0;i<N;i++){
        printf("node %d: %10.9f \n",i+1,u[i]);
        fprintf(output,"%10.9f\n",u[i]);
    }

    //release memory
    free(ae);
    free(ap);
    free(aw);
    free(su);
    free(u);
    free(u0);
    fclose(output);
    return 0;
}
//void Exact_Sol(double,double,int N){

//}

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
