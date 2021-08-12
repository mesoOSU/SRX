#include"stdio.h"
#include"math.h"
#include"stdlib.h"

#define MESHX 128
#define MESHY 128
#define MESHZ 128
#define MESHXYZ MESHX*MESHY*MESHZ
#define t1 (1/sqrt(2))
#define saveT 2000

void readfrmfile (double *x, int s, long t);
//double c_mu(double phi, double mu);
void main(int argc, char *argv[]) {

 int i, j,z, index,s;
 long t;
 double y[MESHXYZ];
 double x[MESHXYZ];
 double EVM[MESHXYZ];
 double SVM[MESHXYZ];
 double eps[MESHXYZ][6];
 double sig[MESHXYZ][6];


 for (t=7000; t <= atol(argv[1]); t+=saveT) {
        for (s=0;s<6;s++){
    readfrmfile (x,s,t);
    for (i=0;i<MESHX;i++){
        for(j =0;j<MESHY;j++) {
            for(z =0;z<MESHZ;z++){
                index = z + (j)*128 + (i)*128*128;
  //  eps[index][s]=y[index];
    sig[index][s]=x[index];
    }
        }
    }
        }

    //Calculate VM
    FILE *fp;
 long count=0;
 char filename[1000];

// sprintf(filename,"EVM_hex_S%04d.txt",t);
// fp=fopen(filename,"w");

FILE *fp1;
 long count1=0;
 char filename1[1000];

 sprintf(filename1,"SVM_491_%04d.txt",t);
 fp1=fopen(filename1,"w");

    for (i=0;i<MESHX;i++){
        for(j =0;j<MESHY;j++) {
            for(z =0;z<MESHZ;z++){
                index = z + (j)*128 + (i)*128*128;
           // if(eps[index][2]<0.0){
    SVM[index]= (t1)*sqrt(pow((sig[index][0]-sig[index][1]),2)+pow((sig[index][0]-sig[index][2]),2)+pow((sig[index][1]-sig[index][2]),2)+6*pow((sig[index][3]),2)+6*pow(sig[index][4],2)+6*pow(sig[index][5],2));
           // else
            //{
    // SVM[i]=(t1)*sqrt(pow((sig[i][0]-sig[i][1]),2)+pow((sig[i][0]-sig[i][2]),2)+pow((sig[i][1]-sig[i][2]),2)+6*pow((sig[i][3]),2)+6*pow(sig[i][4],2)+6*pow(sig[i][5],2));
           // }
           // if(eps[i][2]<0){
  //  EVM[index]= (t1*2/3)*sqrt(pow((eps[index][0]-eps[index][1]),2)+pow((eps[index][0]-eps[index][2]),2)+pow((eps[index][1]-eps[index][2]),2)+6*pow(eps[index][3],2)+6*pow(eps[index][4],2)+6*pow(eps[index][5],2));
           // }
           // else{
   // EVM[i]= (t1*2/3)*sqrt(pow((eps[i][0]-eps[i][1]),2)+pow((eps[i][0]-eps[i][2]),2)+pow((eps[i][1]-eps[i][2]),2)+6*pow(eps[i][3],2)+6*pow(eps[i][4],2)+6*pow(eps[i][5],2));
       /*   if(sig[index][2]<0.0){
              SVM[index] = 0.0 - SVM[index];
          }
          if(eps[index][2]<0.0){
              EVM[index] = 0.0 - EVM[index];
          } */
          
              
              
    // fprintf(fp,"%le\n",EVM[index]);
     fprintf(fp1,"%le\n",SVM[index]);
 }
        }
    }
//  fclose(fp);
  fclose(fp1);
}
}
void readfrmfile(double *x, int s, long t) {
  int i,j,z,index;
  FILE *fp;
  FILE *fp1;
  char filename[1000];
  char filename1[1000];
 /* sprintf(filename,"eps_%d_S%04d.txt",s+1,t);
  fp = fopen(filename,"r");

    for(i=0;i<MESHX;i++) {
        for(j =0;j<MESHY;j++) {
            for(z =0;z<MESHZ;z++){
                index = z + (j)*128 + (i)*128*128;

      fscanf(fp,"%le\n",&y[index]);
    }
        }
    } */
    sprintf(filename1,"sig_%d_491_%04d.txt",s+1,t);
  fp1 = fopen(filename1,"r");

    for(i=0;i<MESHX;i++) {
        for(j =0;j<MESHY;j++) {
            for(z =0;z<MESHZ;z++){
                index = z + (j)*128 + (i)*128*128;

      fscanf(fp1,"%le\n",&x[index]);
    }
        }
    }

//     fscanf(fp,"\n");
  //fclose(fp);
  fclose(fp1);
}

