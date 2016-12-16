#include "GSA.h"
#include <iostream>
#include <cmath>

using namespace std;

double f(double *x,int N){ //Rastrigin function
    double A = 10;
    double s = 0;
    double pi = acos(-1);
    for(int i = 0;i < N; ++i)
        s += x[i]*x[i]-A*cos(2*pi*x[i]);
    return A*N + s;
}

double f2(double *x,int N){ //sum of squares
    double s = 0;
    for(int i = 0; i < N; ++i)
        s += x[i]*x[i];
    return s + 1;
}

double f3(double *x,int N){ //Goldstein-Price
    double s1 = (x[0] + x[1] + 1)*(x[0] + x[1] + 1);
    double s2 = 19-14*x[0] + 3*x[0]*x[0] - 14*x[1] + 6*x[0]*x[1] + 3*x[1]*x[1];
    double s3 = (2*x[0]-3*x[1])*(2*x[0]-3*x[1]);
    double s4 = 18-32*x[0]+12*x[0]*x[0]+48*x[1]-36*x[0]*x[1]+27*x[1]*x[1];
    return (1+s1*s2)*(30+s3*s4);
}

double f4(double *x,int N){ //Rosenbrock
    double s = 0;
    for(int i = 0; i < N -1;++i){
        s += 100*(x[i+1]-x[i]*x[i])*(x[i+1]-x[i]*x[i])+(x[i]-1)*(x[i]-1);
    }
    return s;
}

double f5(double *x,int N){ //Bukin N6
    return 100*sqrt(abs(x[1] - 0.01*x[0]*x[0]))+0.01*abs(x[0]+10);
}

double f6(double *x,int N){ //Cross-in-tray function
    double c = abs(sin(x[0])*sin(x[1])*exp(abs(100-sqrt(x[0]*x[0]+x[1]*x[1])/acos(-1)))+1);
    return -0.0001*pow(c,0.1);
}

double f7(double *x,int N){ //Styblinski–Tang function
    double s = 0;
    for(int i = 0; i < N; ++i)
        s += pow(x[i],4) - 16*x[i]*x[i] + 5*x[i];
    return s/2;
}

double f8(double *x,int N){ //Three-hump camel function
    return 2*x[0]*x[0] - 1.05*pow(x[0],4) + pow(x[0],6)/6+x[0]*x[1]+x[1]*x[1];
}

double f9(double *x,int N){ //Easom function
    return -cos(x[0])*cos(x[1])*exp(-((x[0]-acos(-1))*(x[0]-acos(-1))+(x[1]-acos(-1))*(x[1]-acos(-1))));
}

int main(){
    double *params = new double[8];
    params[0] = 200; //N
    params[1] = 10000; //tmax
    params[2] = 100; //To
    params[3] = 0.9; //pmo
    params[4] = 0.99; //alpha
    params[5] = 0.6; //beta
    params[6] = 5; //K
    params[7] = 100;   //Nitermax
    int ndim = 2;
    double *min = new double[ndim];
    double *max = new double[ndim];
    double *out = new double[ndim];
    for(int i = 0; i < ndim; ++i){
        min[i] = -5;
        max[i] = 5;
    }
    min[0] = -5;
    min[1] = -5;
    max[0] = 5;
    max[1] = 5;
    GSA gsa;
    gsa.set_parameters(params);
    double sum = 0;
    double outt[2];
    outt[0] = 0;
    outt[1] = 0;
    int NN = 200;
    double t =omp_get_wtime();
    for(int i =0; i< NN; ++i){
        gsa.optimize(f,1,ndim,min,max,out);
        sum += f(out,ndim);
        cout << i << endl;
        outt[0] += out[0];
        outt[1] += out[1];
    }
    t = omp_get_wtime() - t;
    cout << t << endl;
    cout << sum/NN << endl;
    for(int i = 0; i < 2; ++i)
        cout << outt[i]/NN << " ";
    delete[] params;
    delete[] min;
    delete[] max;
    delete[] out;
    return 0;
}
