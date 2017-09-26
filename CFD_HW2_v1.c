#include<stdio.h>
#include<stdlib.h>
#include<math.h>

const double l = 1.0;
const double gam= 0.1;
const double rho  = 1.0;
const double finA = 1.0;
const double finB = 0.0;

void TridiagonalSolve(const double *, const double *, double *, double *, double *, unsigned int, FILE *);
void Exact_Sol(double *,double,int,FILE *);

int main(){
    int N;
    int scheme;
    double u;

    FILE *output = fopen("Output_Hw2.txt","w");
    FILE *output2 = fopen("Output_Hw2_Exact.txt","w");
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
    //construct array for storage
    double *ae =(double*)malloc(sizeof(double)*N); //suberdiagonal elements
    double *ap =(double*)malloc(sizeof(double)*N); //diagonal elements
    double *aw =(double*)malloc(sizeof(double)*N); //superdiagonal elements
    double *su =(double*)malloc(sizeof(double)*N); //rhs elements
    double *fin=(double*)malloc(sizeof(double)*N); //solution elements
    double *exact = (double*)malloc(sizeof(double)*N); //exact solution

    int i;
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
    }
    //scheme 3:second order upwind
    else if(scheme==3){
        printf("under construction\n");
        exit(1);
        /*for(i=0;i<N;i++){
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
        }*/
    }
    else{
        printf("error message!!\n");
        exit(1);
    }

    Exact_Sol(exact,u,N,output2);
    printf("\n");
    TridiagonalSolve(ae,ap,aw,su,fin,N,output);

    //release memory
    free(ae);
    free(ap);
    free(aw);
    free(su);
    free(fin);
    free(exact);

    fclose(output);
    fclose(output2);
    return 0;
}

//exact solution
void Exact_Sol(double* exact,double u,int N,FILE *output2){
    int i;
    double x;
    printf("the exact solution is:\n");
    for(i=0;i<N;i++){
        if(i==0)
            x=l/(2*N);
        else
            x+=l/N;
        exact[i]=((exp(rho*u*x/gam)-1)/(exp(rho*u*l/gam)-1))*(finB-finA)+finA;
        printf("node %d: %f\n",i+1,exact[i]);
        fprintf(output2,"%f\n",exact[i]);
    }
}
//Thomos method
void TridiagonalSolve(const double *a, const double *b, double *c, double *d, double *x, unsigned int n,FILE *output){
    int i;

    c[0] = c[0]/b[0];
    d[0] = d[0]/b[0];
    double id;
    for(i = 1; i != n; i++){
        id = 1.0/(b[i] - c[i - 1]*a[i]);
        c[i] = c[i]*id;
        d[i] = (d[i] - a[i]*d[i - 1])*id;
    }

    x[n - 1] = d[n - 1];
    for(i = n - 2; i != -1; i--)
        x[i] = d[i] - c[i]*x[i + 1];

    printf("the numerical solution is:\n");
    //Save data to txt file
    for(i=0;i<n;i++){
        printf("node %d: %f\n",i+1,x[i]);
        fprintf(output,"%f\n",x[i]);
    }
}


