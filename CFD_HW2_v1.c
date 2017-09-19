#include<stdio.h>
#include<stdlib.h>

void TridiagonalSolve(const double *, const double *, double *, double *, double *, unsigned int, FILE *output);

int main(){
    int N;
    const double l = 1.0;
    const double gamma= 0.1;
    const double rho  = 1.0;
    const double finA = 1.0;
    const double finB = 0.0;
    double u;

    FILE *output = fopen("Output_Hw2.txt","w");
    if(output==0){
        printf("error message!!\n");
        exit(1);
    }

    //Define the node number
    printf("enter node number:\n");
    scanf("%d",&N);
    printf("enter velocity:\n");
    scanf("%lf",&u);

    double F=rho*u;
    double D=(gamma*N)/l;
    //construct array for storage
    double *aw =(double*)malloc(sizeof(double)*N); //superdiagonal elements
    double *ap =(double*)malloc(sizeof(double)*N); //diagonal elements
    double *ae =(double*)malloc(sizeof(double)*N); //subdiagonal elements
    double *su =(double*)malloc(sizeof(double)*N); //rhs elements
    double *fin=(double*)malloc(sizeof(double)*N); //solution elements

    int i;
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

    //for(i=0;i<N;i++){
    //    printf("%f\n",ae[i]);
    //}
    TridiagonalSolve(ae,ap,aw,su,fin,N,output);

    //release memory
    free(ae);
    free(ap);
    free(aw);
    free(su);
    free(fin);

    fclose(output);
    return 0;
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

    //Save data to txt file
    for(i=0;i<n;i++){
        printf("%f\n",x[i]);
        fprintf(output,"%f\n",x[i]);
    }
}
