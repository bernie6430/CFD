/*===========================================================
Simple example for Newton Raphson method for solving equation
Author:Bo Han Huang
Date:9/12/2017
Version:1.0
=============================================================*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>


int main(){
    double xc;
    double xt;
    double y;
    double y_dif;
    double record[2];
    int i;

    printf("enter initial value:");
    scanf("%lf",&xt);

   while(1){
        record[0]=xt;
        y=xt+2-exp(xt);
        y_dif=1-exp(xt);
        xc=xt-y/y_dif;
        record[1]=xc;
        if(record[0]-record[1]<0.0001)
            break;
        xt=xc;
    }
    printf("%f",xc);
    return 0;
}
