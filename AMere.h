#ifndef _MERE_H_
#define _MERE_H_


/*                                                                                                                                            
    @requires:rien                                                                                                                            
    @assigns:rien                                                                                                                             
    @ensures:                                                                                                                                 
*/


class Aclass {

    double T_; // temps max                                                                                                                   
    double S_; // veleur spatial maximum                                                                                                      
    int M_; // nombre de points temporels (du maillage)                                                                                       
    int N_; // nombre de points spatiales(du maillage)                                                                                        
    double mu_;


    public:
        Aclass(double T, double S, int M, int N, double mu);
        // ~Aclass(){};                                                                                                                        
        double get_T() const;
        double get_S() const;
        int get_M() const;
        int get_N() const;
        double get_mu() const;

        void vecmult(double **A1, double*C,double *k,double *vec1,int N); // pour multiplication entre
        //une matrice carré et un vecteur et le k ?
        double max1(double a, double b);
        void Resolution(double **A,double*b,double *x, int n);// resout Ax=b pour A tridiagonale


        virtual void Methode (double K,double r,double sigma)=0; //=0va etre implementer dans les filles pr compléter et ontenir Z à la fin                                                                                        

};



#endif