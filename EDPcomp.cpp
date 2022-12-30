#include "EDPcomp.h"
#include <math.h>



//////////////////////////////////////////////////////////////////////////////////////////////
///////////////// PARTIE POUR LE PUT///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

// Constructeur                                                                                                                               
CranckPut::CranckPut(double T, double S,int M, int N,double mu): Aclass(T,  S, M,  N,  mu) {

    data_= new double*[M+1];
    for (int i = 0; i < M+1; i++){
        data_[i] = new double[N+1];
    }


    if(!data_){ throw ("Erreur lors de l'allocation mémoire");}
}


//Destructezur                                                                                                                                
CranckPut::~CranckPut(){
    for (int i = 0; i < (get_M()+1); i++){delete[] data_[i];}
    delete[] data_;

}

void CranckPut::Methode (double K,double r,double sigma){

    int N=get_N();
    int M=get_M();
    //double mu=get_mu();
    double T=get_T();
    double S=get_S();
    double dt=(double) T/M;
    double ds=(float) S/N;

    /// calcul des points aux limites ///                                                                                                         

    // 2ème et 3 ème CL                                                                                                                           
    for(int i = 0; i < (M+1); i++){// on prends pas la derniere ligne d'où le M                                                                       
        data_[i][0]=K*exp((-1)*r*(T-(dt*i)));
        data_[i][N]=0;
    }
    //1ère CL:                                                                                                                                    
    for(int j = 0; j < (N+1); j++){data_[M][j]=max1(0.0,K-(ds*j));}



    //Allocations mémoires:
    double *a=new double[N-1];
    double *b=new double[N-1];
    double *c=new double[N-1];
    double *d=new double[N-1];

    double **A1=new double*[N-1];
    double **A2=new double*[N-1];
    for (int i = 0; i <(N-1); i++){
        A1[i] = new double[N-1];
        A2[i] = new double[N-1];
    }

    double *C=new double[N-1];
    double *k=new double[N-1];
    double *vec1=new double[N-1];

    // calcul des coefficients a,b,c,d puis construction des matrices A1 et B2 et K

    // a,b,c,d
    for (int j = 0; j < (N-1); j++){
        a[j]=0.25*(j+1)*dt*(pow(sigma,2)*(j+1)-r);
        b[j]=(1-0.5*pow(sigma*(j+1),2)*dt);
        c[j]=0.25*(j+1)*dt*(pow(sigma,2)*(j+1)-r);
        d[j]=(1+(r+0.5*pow(sigma*(j+1),2))*dt);
    }

    // A1 et A2:
    for (int i = 0; i < (N-1); i++){
        for (int j = 0; j < (N-1); j++){
            A1[i][j]=0.0;
            A2[i][j]=0.0;
        }
    }
    for (int i = 0; i < (N-1); i++){A1[i][i]=b[i];A2[i][i]=d[i];}
    for (int i = 0; i < (N-2); i++){
        A1[i][i+1]=c[i];
        A1[i+1][i]=a[i+1];
        A2[i][i+1]=(-1.0)*c[i];
        A2[i+1][i]=(-1.0)*a[i+1];
    }


    // initianilisation de C avce la première CL:
    for(int j=0; j<(N-1); j++){C[j]=data_[j+1][M];}

    // la matrice K:
    //k[0]=a[0]*(2.0*K);
    for (int j = 0; j < (N-1); j++){k[j]=0.0;}


    //// Boucle calculant les Pm au fur et à mesure///
    // et les stockans dans data_

    for(int j=M; j>0;j--){
        k[0]=a[0]*(data_[j-1][0]+data_[j][0]);
        k[N-2]=c[N-2]*(data_[j-1][N]+data_[j][N]);

        vecmult(A1,C,k,vec1,M);
        Resolution(A2,vec1,C,M);

        for (int i = 0; i < (N-1); i++){data_[j-1][i+1]=C[i];}
        
    }
    
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] C;
    delete[] k;
    delete[] vec1;
    for(int j=0; j<(N-1);j++){delete[] A1[j]; delete[] A2[j];}
    delete[] A1;
    delete[] A2;


}


//////////////////////////////////////////////////////////////////////////////////////////////
///////////////// PARTIE POUR LE CALLL///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////


// Constructeur                                                                                                                               
CranckCall::CranckCall(double T, double S,int M, int N,double mu): Aclass(T,  S, M,  N,  mu) {

    data_= new double*[M+1];
    for (int i = 0; i < M+1; i++){
        data_[i] = new double[N+1];
    }


    if(!data_){ throw ("Erreur lors de l'allocation mémoire");}
}


//Destructezur                                                                                                                                
CranckCall::~CranckCall(){
    for (int i = 0; i < (get_M()+1); i++){delete[] data_[i];}
    delete []data_;

}

void CranckCall::Methode (double K,double r,double sigma){

    int N=get_N();
    int M=get_M();
    //double mu=get_mu();
    double T=get_T();
    double S=get_S();
    double dt=(double) T/M;
    double ds=(float) S/N;

    /// calcul des points aux limites ///                                                                                                         

    // 2ème et 3 ème CL                                                                                                                           
    for(int i = 0; i < (M+1); i++){// on prends pas la derniere ligne d'où le M                                                                       
        data_[i][0]=0;
        data_[i][N]=K*exp((-1)*r*(T-(dt*i)));
    }
    //1ère CL:                                                                                                                                    
    for(int j = 0; j < (N+1); j++){data_[M][j]=max1(0.0,(ds*j)-K);}



    //Allocations mémoires:
    double *a=new double[N-1];
    double *b=new double[N-1];
    double *c=new double[N-1];
    double *d=new double[N-1];

    double **A1=new double*[N-1];
    double **A2=new double*[N-1];
    for (int i = 0; i <(N-1); i++){
        A1[i] = new double[N-1];
        A2[i] = new double[N-1];
    }

    double *C=new double[N-1];
    double *k=new double[N-1];
    double *vec1=new double[N-1];

    // calcul des coefficients a,b,c,d puis construction des matrices A1 et B2 et K

    // a,b,c,d
    for (int j = 0; j < (N-1); j++){
        a[j]=0.25*(j+1)*dt*(pow(sigma,2)*(j+1)-r);
        b[j]=(1-0.5*pow(sigma*(j+1),2)*dt);
        c[j]=0.25*(j+1)*dt*(pow(sigma,2)*(j+1)-r);
        d[j]=(1+(r+0.5*pow(sigma*(j+1),2))*dt);
    }

    // A1 et A2:
    for (int i = 0; i < (N-1); i++){
        for (int j = 0; j < (N-1); j++){
            A1[i][j]=0.0;
            A2[i][j]=0.0;
        }
    }
    for (int i = 0; i < (N-1); i++){A1[i][i]=b[i];A2[i][i]=d[i];}
    for (int i = 0; i < (N-2); i++){
        A1[i][i+1]=c[i];
        A1[i+1][i]=a[i+1];
        A2[i][i+1]=(-1.0)*c[i];
        A2[i+1][i]=(-1.0)*a[i+1];
    }


    // initianilisation de C avce la première CL:
    for(int j=0; j<(N-1); j++){C[j]=data_[j+1][M];}

    // la matrice K:
    //k[0]=a[0]*(2.0*K);
    for (int j = 0; j < (N-1); j++){k[j]=0.0;}


    //// Boucle calculant les Pm au fur et à mesure///
    // et les stockans dans data_

    for(int j=M; j>0;j--){
        k[0]=a[0]*(data_[j-1][0]+data_[j][0]);
        k[N-2]=c[N-2]*(data_[j-1][N]+data_[j][N]);

        vecmult(A1,C,k,vec1,M);
        Resolution(A2,vec1,C,M);

        for (int i = 0; i < (N-1); i++){data_[j-1][i+1]=C[i];}
        
    }

    // On désalloue la mémoire

    delete []a;
    delete []b;
    delete []c;
    delete []d;
    delete []C;
    delete []k;
    delete []vec1;
    for(int j=0; j<(N-1);j++){delete A1[j]; delete A2[j];}
    delete []A1;
    delete []A2;


}