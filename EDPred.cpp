#include "EDPred.h"
#include <math.h>



//https://h-deb.clg.qc.ca/Sujets/Divers--cplusplus/CPP--Tableaux-2D.html                                                                      



// Constructeur                                                                                                                               
DiffImpPut::DiffImpPut(double T, double S,int M, int N,double mu): Aclass(T,  S, M,  N,  mu) {

    data_= new double*[M+1];
    for (int i = 0; i < M+1; i++){
        data_[i] = new double[N+1];
    }


    if(!data_){ throw ("Erreur lors de l'allocation mémoire");}
}


//Destructezur                                                                                                                                
DiffImpPut::~DiffImpPut(){
    for (int i = 0; i < (get_M()+1); i++){delete[] data_[i];}
    delete[] data_;

}



// méthode qui va bien définr data_ avec les vraies valeurs (data_ est notre maytrice Z)                                                      
// on va faire de l'allocation dynamique et retourner un (double?) pointeur                                                                   
// il faudra ensuite afficher la courbe avec une autre fonction
void DiffImpPut::Methode(double K,double r,double sigma){



/// Definition de la matrice A / des listes coreepondant aux coeffs ab,c (b=c) ///                                                            
int N=get_N();
int M=get_M();
double mu=get_mu();
double T=get_T();
double S=get_S();
double dt=(double) T/M;
double ds=(float) S/N;

double** A=new double*[N+1];
for (int i = 0; i < N+1; i++){
    A[i] = new double[N+1];
}
double *vec1=new double[N+1];
double *C=new double[N+1];


// Initialisations des variables

// initialisation de A:
for(int i = 0; i < (N+1); i++){
    for(int j = 0; j < (M+1); j++){
        if(i==j){A[i][j]= 1 + (2*mu*dt*dt)/(ds*ds);}
	    else if(i==j+1 || j==i+1){A[i][j]=(-1)*(dt*dt)/(ds*ds);}
	    else {A[i][j]= 0;}
    }
}





/// calcul des points aux limites ///                                                                                                         

// CL                                                                                                                           
for(int i = 0; i < (M+1); i++){data_[0][i]=exp( ((2*r/sigma*sigma)-1)*i*ds/2*max1(0.0,(exp(i*ds)-1) ) );}


///  Calcul des éléments possibles grâce à la méthode ///                                                                                     

 // M=N dans nos essais                                                                                                   
for(int i = 1 ; i < M+1; i++){
    for(int j = 0; j <(M+1); j++){vec1[j]=data_[i-1][j];}
    Resolution(A,vec1,C,M+2);
    for (int j = 0; j < (M+1); j++){data_[i][j]=C[j];}

}




/// suppresion allocations mémoires ///                                                                                                       


for (int i = 0; i < N-1; i++){delete[] A[i];}
delete[] A;
delete[] C;
delete[] vec1;

}


