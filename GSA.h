#ifndef GSA_H
#define GSA_H

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "omprng.h"

using namespace std;

class GSA
{
private:
	typedef struct {
        double value;
        int num;
    } Fset;
    bool terminaton(Fset*);
	void init_parents(Fset*,int);
	void qs(Fset*,int,int);
	int psearch(int,double);
    int N;
    long long tmax;
    double To;
    double pmo;
    double alpha;
    double beta;
    double K;
    bool paramsset;
    int Niter;
    int Nitermax;
    long long t;
    omprng MyRng;
public:
    GSA(void);
    ~GSA(void);
    void optimize(double f(double*,int),int,int,double*, double*,double*);
    void set_parameters(double*);
	void optimizeSA(double f(double*,int), double*, int, double, double, double*, double*, double*, double&);
};

#endif
