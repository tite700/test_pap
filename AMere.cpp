#include "AMere.h"


Aclass::Aclass(double T, double S, int M, int N, double mu){
    T_= T;
    S_=S;
    M_=M;
    N_=N;
    mu_=mu;
}


double Aclass::get_T() const{return T_;}
double Aclass::get_S() const{return S_;}
int Aclass::get_M() const{return M_;}
int Aclass::get_N() const{return N_;}
double Aclass::get_mu() const{return mu_;}

void Aclass::vecmult(double **A1, double*C,double *k,double *vec1,int N){
    for(int i=0;i<(N-1);i++){
        double tmp=0;
        for(int j=0;j<(N-1);j++){
            tmp=tmp+A1[i][j]*C[j];
            vec1[i]= tmp+k[i];
        }

    }
}


double Aclass::max1(double a, double b){return ((a>b)?a:b);}


void Aclass::Resolution(double **A,double*b,double *x, int n){
    n=n-1;// on va donner N en argument alors qu'on veut N-1

    //Allocation mémoire:
    double *z=new double[n];
    double **L=new double*[n];
    for(int i=0;i<n;i++){L[i]=new double[n];}
    double **U=new double*[n];
    for(int i=0;i<n;i++){U[i]=new double[n];}

    // Facorisation LU:
    L[0][0]=A[0][0];
    U[0][1]=A[0][1]/L[0][0];
    for(int i=1;i<=(n-2);i++){
        L[i][i-1]=A[i][i-1];
        L[i][i]=A[i][i]-L[i][i-1]*U[i-1][i];
        U[i][i+1]=A[i][i+1]/L[i][i];
    }
    L[n-1][n-2]=A[n-1][n-2];
    L[n-1][n-1]=A[n-1][n-1]-L[n-1][n-2]*U[n-2][n-1];
    // On resout Lz=b avec z=Ux
    z[0]=b[0]/L[0][0];
    for(int i=1;i<=(n-1);i++){
        z[i]=(b[i]-L[i][i-1]*z[i-1])/L[i][i];
    }
    // On resout Ux=z
    x[n-1]=z[n-1];
    for(int i=(n-2);i>=0;i--){
        x[i]=z[i]-U[i][i+1]*x[i+1];
    }
    // on désalloue la mémoire
    delete[]z;
    for(int j=0;j<n;j++){
        delete[] L[j];
        delete[] U[j];
    }
    delete[] L;
    delete[] U;

 }


