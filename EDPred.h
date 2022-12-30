#ifndef _EDPred_H_
#define _EDPred_H_

#include "AMere.h"


class DiffImpPut : public Aclass {
    double **data_;

    public:
        DiffImpPut(double T, double S,int M, int N,double mu);
        ~DiffImpPut();
        double operator()(int l, int c) const { return data_[l][c]; }


        void Methode (double K,double r,double sigma);

};

#endif