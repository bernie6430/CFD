#include<stdio.h>
#include<stdlib.h>
#include<math.h>

const double l = 1.0;
const double gam= 0.1;
const double rho  = 1.0;
const double finA = 1.0;
const double finB = 0.0;

void TridiagonalSolve(const double *, const double *, double *, double *, double *, unsigned int);
void Exact_Sol(double *,double,int);

int main(){
    int N;
    int scheme;
    int cas;
    double u;

    FILE *output = fopen("Output_Hw2.txt","w");
    FILE *output2 = fopen("Output_Hw2_Exact.txt","w");
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
    printf("enter scheme:(1.CD 2.FOU 3.SOU)\n");
    scanf("%d",&scheme);


    double F=rho*u;
    double D=(gam*N)/l;
    double Pe=F/D;
    //construct array for storage
    double *ae =(double*)malloc(sizeof(double)*N); //suberdiagonal elements
    double *ap =(double*)malloc(sizeof(double)*N); //diagonal elements
    double *aw =(double*)malloc(sizeof(double)*N); //superdiagonal elements
    double *su =(double*)malloc(sizeof(double)*N); //rhs elements
    double *fin=(double*)malloc(sizeof(double)*N); //solution elements
    double *exact = (double*)malloc(sizeof(double)*N); //exact solution

    double *aww =(double*)malloc(sizeof(double)*(N-2));
    double *tem =(double*)malloc(sizeof(double)*(N-2));

    int i,j,k;
    for(i=0;i<(N-2);i++)
        tem[i]=1.0;
    //scheme 1;central difference
    if(scheme==1){
        for(i=0;i<N;i++){
            //First node
            if(i==0){
                ae[i]=0.0;
                ap[i]=3.0*D+F/2.0;
                aw[i]=-(D-F/2.0);
                su[i]=(2.0*D+F)*finA;
                fin[i]=0.0;
            }
            //Last Node
            else if(i==(N-1)){
                ae[i]=-(D+F/2.0);
                ap[i]=3.0*D-F/2.0;
                aw[i]=0.0;
                su[i]=(2.0*D-F)*finB;
                fin[i]=0.0;
            }
            //Other Nodes
            else{
                ae[i]=-(D+F/2.0);
                ap[i]=2.0*D;
                aw[i]=-(D-F/2.0);
                su[i]=0.0;
                fin[i]=0.0;
            }
        }
        TridiagonalSolve(ae,ap,aw,su,fin,N);
    }
    //scheme 2: first order upwind
    else if(scheme==2){
        for(i=0;i<N;i++){
            //First node
            if(i==0){
                ae[i]=0.0;
                ap[i]=3.0*D+F;
                aw[i]=-D;
                su[i]=(2.0*D+F)*finA;
                fin[i]=0.0;
            }
            //Last Node
            else if(i==(N-1)){
                ae[i]=-(D+F);
                ap[i]=3.0*D+F;
                aw[i]=0.0;
                su[i]=2.0*D*finB;
                fin[i]=0.0;
            }
            //Other Nodes
            else{
                ae[i]=-(D+F);
                ap[i]=2.0*D+F;
                aw[i]=-D;
                su[i]=0.0;
                fin[i]=0.0;
            }
        }
        TridiagonalSolve(ae,ap,aw,su,fin,N);
    }
    //scheme 3:second order upwind
    else if(scheme==3){
        printf("choose boundary condition:1.fixed 2.FOU\n");
        scanf("%d",&cas);
        //iterative TDMA
        for(j=0;j<1;j++){
            for(i=0;i<N;i++){
                //First node
                if(i==0){
                    ae[i]=0.0;
                    ap[i]=3.0*D+2.0*F;
                    aw[i]=-D;
                    su[i]=2.0*F+2.0*D;
                    fin[i]=0.0;
                }
                //second node
                else if(i==1){
                    ae[i]=-(D+2.0*F);
                    ap[i]=2.0*D+1.5*F;
                    aw[i]=-D;
                    su[i]=-0.5*F;
                    fin[i]=0.0;
                }
                //Last Node
                else if(i==(N-1)){
                    if(cas==1){
                        ae[i]=-(D+1.5*F);
                        ap[i]=3.0*D;
                        aw[i]=0.0;
                        aww[i-2]=0.5*F;
                        fin[i]=0.0;
                        su[i]=-0.5*F*tem[i-2];
                    }
                    else if(cas==2){
                        ae[i]=-(D+F);
                        ap[i]=F+3.0*D;
                        aw[i]=0.0;
                        aww[i-2]=0.0;
                        fin[i]=0.0;
                        su[i]=-aww[i-2]*tem[i-2];
                    }
                    else{
                        printf("error message!!\n");
                        exit(1);
                    }
                }
                //Other Nodes
                else{
                    ae[i]=-(D+2.0*F);
                    ap[i]=2.0*D+1.5*F;
                    aw[i]=-D;
                    aww[i-2]=0.5*F;
                    fin[i]=0.0;
                    su[i]=-0.5*F*tem[i-2];
                }
            }
            TridiagonalSolve(ae,ap,aw,su,fin,N);
            for(k=0;k<(N-2);k++)
                tem[k]=fin[k];
        }
    }
    else{
        printf("error message!!\n");
        exit(1);
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
    free(aww);
    free(tem);

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



