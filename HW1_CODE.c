#include<stdio.h>
#include<stdlib.h>


void TridiagonalSolve(const double *, const double *, double *, double *, double *, unsigned int, FILE *output);

int main(){
    int N;
    int Ca;
    const double k = 0.5;
    const double S = 1000000.0;
    const double TA = 100.0;
    const double TB2 = 100.0;
    const double TB = 200.0;
    const double Tin = 20.0;
    FILE *output = fopen("Output.txt","w");
    if(output==0){
        printf("error message!!\n");
        exit(1);
    }

    //Define the case
    printf("enter the case number:");
    scanf("%d",&Ca);
    //Define the node number
    printf("enter node number:");
    scanf("%d",&N);

    //construct array for storage
    double *a=(double*)malloc(sizeof(double)*N); //subdiagonal elements
    double *b=(double*)malloc(sizeof(double)*N); //diagonal elements
    double *c=(double*)malloc(sizeof(double)*N); //superdiagonal elements
    double *d=(double*)malloc(sizeof(double)*N); //rhs elements
    double *x=(double*)malloc(sizeof(double)*N); //solution elements

    //Case 1
    if(Ca==1){
        int i;
        double l = 0.02;
        double delx = l/N;
        for(i=0;i<N;i++){
            //First node
            if(i==0){
                a[i]=0;
                b[i]=(3*k)/delx;
                c[i]=-k/delx;
                d[i]=S*delx+(2*k*TA)/delx;
                x[i]=0;
            }
        //Last Node
            else if(i==(N-1)){
                a[i]=-k/delx;
                b[i]=(3*k)/delx;
                c[i]=0;
                d[i]=S*delx+(2*k*TB)/delx;
                x[i]=0;
            }
            //Other Nodes
            else{
                a[i]=-k/delx;
                b[i]=(2*k)/delx;
                c[i]=-k/delx;
                d[i]=S*delx;
                x[i]=0;
            }
        }
    }
    //Case 2
    else if(Ca==2){
        int i;
        double l = 1;
        double delx = l/N;
        double n = 5.0;                  //n^2=(hP/kA)=25
        for(i=0;i<N;i++){
            //First node
            if(i==0){
                a[i]=0;
                b[i]=3/delx+n*n*delx;
                c[i]=-1/delx;
                d[i]=n*n*delx*Tin+(2*TB2)/delx;
                x[i]=0;
            }
            //Last Node
            else if(i==(N-1)){
                a[i]=-1/delx;
                b[i]=1/delx+n*n*delx;
                c[i]=0;
                d[i]=n*n*delx*Tin;
                x[i]=0;
            }
            //Other nodes
            else{
                a[i]=-1/delx;
                b[i]=2/delx+n*n*delx;
                c[i]=-1/delx;
                d[i]=n*n*delx*Tin;
                x[i]=0;
            }
        }
    }
    else{
        printf("error message!!\n");
    }

    TridiagonalSolve(a,b,c,d,x,N,output);

    //release memory
    free(a);
    free(b);
    free(c);
    free(d);
    free(x);

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
