#include<stdio.h>
#include<stdlib.h>
#include<math.h>

const double k = 1000;
const double delx = 0.1;
const double dely = 0.1;
const double A=0.001;
const double T = 100.0;
const double Q = 500000.0;

void TridiagonalSolve(const double *, const double *, double *, double *, double *, unsigned int);
void Exact_Sol(double *,double,int);

int main(){
    int N;
    int col;
    int row;

    FILE *output = fopen("Output_2D.txt","w");
    if(output==0){
        printf("error message!!\n");
        exit(1);
    }

    //Define the node number
    printf("enter row node number:\n");
    scanf("%d",&row);
    printf("enter column node number:\n");
    scanf("%d",&col);

    N=row*col;
    double C=(k/delx)*A;
    double D=((2*k)/dely)*A;

    //construct array for storage
    double *aw =(double*)malloc(sizeof(double)*row); //superdiagonal elements
    double *ap =(double*)malloc(sizeof(double)*row); //diagonal elements
    double *ae =(double*)malloc(sizeof(double)*row); //subdiagonal elements
    double *su =(double*)malloc(sizeof(double)*row); //rhs elements
    double *fin_col=(double*)malloc(sizeof(double)*col);
    double *fin=(double*)malloc(sizeof(double)*N); //solution elements

    int i,j,k,s;
    for(i=0;i<N;i++)
        fin[i]=0.0;


    //iterative TDMA
    for(k=0;k<100;k++){
        for(i=0;i<row;i++){
            if(i==0){
                for(j=0;j<col;j++){
                //First node
                    if(j==0){
                        aw[j]=0.0;
                        ap[j]=2.0*C;
                        ae[j]=-C;
                        su[j]=C*fin[(i+1)*col+j]+Q*A;
                    }
                    //Last Node
                    else if(j==(col-1)){
                        aw[j]=-C;
                        ap[j]=2.0*C;
                        ae[j]=0.0;
                        su[j]=C*fin[(i+1)*col+j];
                    }
                    //Other Nodes
                    else{
                        aw[j]=-C;
                        ap[j]=3.0*C;
                        ae[j]=-C;
                        su[j]=C*fin[(i+1)*col+j];
                    }
                }
                for(s=0;s<col;s++)
                    fin_col[s]=0.0;
                TridiagonalSolve(aw,ap,ae,su,fin_col,col);
                for(s=0;s<col;s++)
                    fin[i*col+s]=fin_col[s];
            }
            else if(i==(row-1)){
                for(j=0;j<col;j++){
                //First node
                    if(j==0){
                        aw[j]=0.0;
                        ap[j]=4.0*C;
                        ae[j]=-C;
                        su[j]=C*fin[(i-1)*col+j]+Q*A+D*T;
                    }
                    //Last Node
                    else if(j==(col-1)){
                        aw[j]=-C;
                        ap[j]=4.0*C;
                        ae[j]=0.0;
                        su[j]=C*fin[(i-1)*col+j]+D*T;
                    }
                    //Other Nodes
                    else{
                        aw[j]=-C;
                        ap[j]=5.0*C;
                        ae[j]=-C;
                        su[j]=C*fin[(i-1)*col+j]+D*T;
                    }
                }
                for(s=0;s<col;s++)
                    fin_col[s]=0.0;
                TridiagonalSolve(aw,ap,ae,su,fin_col,col);
                for(s=0;s<col;s++)
                    fin[i*col+s]=fin_col[s];
            }
            else{
                for(j=0;j<col;j++){
                //First node
                    if(j==0){
                        aw[j]=0.0;
                        ap[j]=3.0*C;
                        ae[j]=-C;
                        su[j]=C*fin[(i+1)*col+j]+C*fin[(i-1)*col+j]+Q*A;
                    }
                    //Last Node
                    else if(j==(col-1)){
                        aw[j]=-C;
                        ap[j]=3.0*C;
                        ae[j]=0.0;
                        su[j]=C*fin[(i+1)*col+j]+C*fin[(i-1)*col+j];
                    }
                    //Other Nodes
                    else{
                        aw[j]=-C;
                        ap[j]=4.0*C;
                        ae[j]=-C;
                        su[j]=C*fin[(i+1)*col+j]+C*fin[(i-1)*col+j];
                    }
                }
                for(s=0;s<col;s++)
                    fin_col[s]=0.0;
                TridiagonalSolve(aw,ap,ae,su,fin_col,col);
                for(s=0;s<col;s++)
                    fin[i*col+s]=fin_col[s];
            }
        }
    }


    printf("\n");


    printf("numerical solution :            analytical solution:\n");
    //Save data to txt file
    for(i=0;i<N;i++){
        printf("node %d: %f %25f\n",i+1,fin[i]);
        fprintf(output,"%f\n",fin[i]);
    }

    //release memory
    free(ae);
    free(ap);
    free(aw);
    free(su);
    free(fin_col);
    free(fin);

    fclose(output);
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
