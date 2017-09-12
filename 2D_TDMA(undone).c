#include<stdio.h>
#include<stdlib.h>

typedef struct node_info{
    double an;
    double as;
    double ae;
    double aw;
    double ap;
    double du;
}info;

void define_node(info* node_con,int node);


int main(){
    int i,j;
    int node;
    double *beta,*diago,*alpha,*A,*C;
    info *node_con;
    printf("enter the node number:");
    scanf("%d",&node);

    beta=(double*)malloc(node*sizeof(double));
    diago=(double*)malloc(node*sizeof(double));
    alpha=(double*)malloc(node*sizeof(double));
    A=(double*)malloc(node*sizeof(double));
    C=(double*)malloc(node*sizeof(double));
    node_con=(info*)malloc(node*sizeof(info));


    for(i=0;i<node;i++){
        beta[i]=0.0;
        diago[i]=0.0;
        alpha[i]=0.0;
        A[i]=0.0;
        C[i]=0.0;
    }


    free(beta);
    free(diago);
    free(alpha);
    free(A);
    free(C);
    free(node_con);
    return 0;
}

void define_node(info* node_con,int node){
    const double k=1000.0;
    const double q=500000.0;
    const double temp=100.0;
    const double tck=0.01;
    const double delx=0.1;
    const double dely=0.1;
    int n,s,w,e,p,a;
    int i,j;

    for(i=0;i<node;i++){
        if((i+1)%node==0){
            node_con[i].an=0.0;
            node_con[i].su=(2*k*0.01*0.1*temp)/dely;
        }
        else if((i-1)%node==0)
            node_con[i].as=0.0;
        else if(i<4){
            node_com[i].aw=0.0;
            node_com[i].su=q*0.01*0.1;
        }
        else if(i>7)
            node_com[i].ae=0.0;
        else{
            node_com[i].an=10.0;
            node_com[i].as=10.0;
            node_com[i].aw=10.0;
            node_com[i].ae=10.0;
        }
    }
}
