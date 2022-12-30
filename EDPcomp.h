
#ifndef _EDPcomp_H_
#define _EDPcomp_H_

#include "AMere.h"


class CranckPut : public Aclass {
    double **data_;

    public:
        CranckPut(double T, double S,int M, int N,double mu);
        ~CranckPut();

        double operator()(int l, int c) const { return data_[l][c]; }
        

        void Methode (double K,double r,double sigma);

};



class CranckCall : public Aclass {
    double **data_;

    public:
        CranckCall(double T, double S,int M, int N,double mu);
        ~CranckCall();
        double operator()(int l, int c) const { return data_[l][c]; }

        void Methode (double K,double r,double sigma);

};


#endif