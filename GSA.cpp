#include "GSA.h"

using namespace std;
GSA::GSA(void){
    paramsset = false;
}


GSA::~GSA(void)
{
}

void GSA::set_parameters(double *params){
    paramsset = true;
    N = int(params[0]);
    if (N%2 == 1) ++N;
    tmax = params[1];
    To = params[2];
    pmo = params[3];
    alpha = params[4];
    beta = params[5];
    K = params[6];
    Nitermax = int(params[7]);
}

bool GSA::terminaton(Fset *fs){
    bool ans = true;
    for(int i = 1; i < N; ++i){
        if(abs(fs[i].value-fs[i-1].value) > 1e-10){
            Niter = 0;
            ans = false;
            break;
        }
    }
    ++Niter;
    if(Niter < Nitermax && t < tmax)
        return false;
    return true;
}

void GSA::qs(Fset* s_arr, int first, int last)
{
    int i = first, j = last;
    double x = s_arr[(first + last) / 2].value;

    do {
        while (s_arr[i].value > x) i++;
        while (s_arr[j].value < x) j--;

        if(i <= j) {
            if (s_arr[i].value < s_arr[j].value){
                s_arr[i].value += s_arr[j].value;
                s_arr[j].value = s_arr[i].value - s_arr[j].value;
                s_arr[i].value -= s_arr[j].value;
                s_arr[i].num += s_arr[j].num;
                s_arr[j].num = s_arr[i].num - s_arr[j].num;
                s_arr[i].num -= s_arr[j].num;
            }
            i++;
            j--;
        }
    } while (i <= j);

    if (i < last)
        qs(s_arr, i, last);
    if (first < j)
        qs(s_arr, first, j);
}

int GSA::psearch(int N,double x){
    double c = 2.0/(N*(N+1));
    for(int i = 0; i < N; ++i){
        x -= (i+1) * c;
        if(x < 0)
            return i;
    }
}

//TODO: переложить параметры в параметры класса +
//Распараллеливание
//Критерий завершения цикла переписать
//Проверить на традиционных функциях
void GSA::optimize(double f(double*,int),int Nthreads,int dim,double *min, double *max,double *out){
    if(!paramsset){
        cout <<"GSA parameters not set. Optimization failed" << endl;
        return;
    }
    Niter = 0;
    t = 0;
    omp_set_num_threads(Nthreads);
    //init population
    double ** population = new double*[N];
    for(int i = 0; i < N; ++i){
        population[i] = new double[dim];
		for(int j = 0; j < dim; ++j)
            population[i][j] = MyRng.runif(min[j],max[j]);
	}
    //function values
	Fset * parents = new Fset[N];
	double ** children = new double*[N];
	double * fvalues = new double[N];
	for(int i = 0; i < N; ++i){
		children[i] = new double[dim];
		fvalues[i] = f(population[i],dim);
	}
    //sobsna, cikl
    while(!terminaton(parents)){
        ++t;
        //Ranking selection method
        for(int i = 0; i < N; ++i){
            parents[i].value = fvalues[i];
            parents[i].num = i;
        }
        qs(parents,0,N-1);
#pragma omp parallel for
        for(int i = 0; i < N/2; ++i){
            //choose two random parents
            double r3 = MyRng.runif();
            double r4 = MyRng.runif();
            int p1 = psearch(N,r3);
            int p2 = psearch(N,r4);
            //arithmetic crossover
            double r1 = MyRng.runif();
            for(int j = 0; j < dim; ++j){
                children[2*i][j] = r1 * population[p1][j] + (1-r1) * population[p2][j];
                children[2*i + 1][j] = (1-r1) * population[p1][j] + r1 * population[p2][j];
            }
            //nonuniform mutation
            for(int j = 0; j < 2; ++j){
                double r2 = MyRng.runif();
                if(r2 < pmo){
                    int gene = MyRng.runifint(0,dim-1);
                    int sign = MyRng.runifint(0,1);
                    double r3 = MyRng.runif();
                    double delta = 1-pow(r3,pow((1-t/tmax),beta));
                    delta *= sign == 0 ? (min[gene] - children[2*i + j][gene]) :  (max[gene] - children[2*i + j][gene]);
                    children[2*i + j][gene] += delta;
                }
            }
            //Boltzmann Trial
            double f1 = f(children[2*i],dim);
            double f2 = f(children[2*i+1],dim);
            double r = MyRng.runif();
            if(exp((parents[p1].value-f1)/To) < r){
                for(int j = 0; j < dim; ++j)
                    children[2*i][j] = population[parents[p1].num][j];
				fvalues[2*i] = parents[p1].value;
			}
			else
				fvalues[2*i] = f1;
            r = MyRng.runif();
            if(exp((parents[p2].value-f2)/To) < r){
                for(int j = 0; j < dim; ++j)
                    children[2*i+1][j] = population[parents[p2].num][j];
				fvalues[2*i+1] = parents[p2].value;
			}
			else
				fvalues[2*i+1] = f2;
        }
        for(int i = 0; i < N; ++i)
            for(int j = 0; j < dim; ++j)
                population[i][j] = children[i][j];
        //decrease Temperature and mutation probability
        To *= alpha;
        if( t % int(K) == 0 && pmo > 1.0/N)
            pmo *= alpha;
    }
    int nbest = 0;
    double fvalmin = fvalues[0];
    for(int i = 1; i < N; ++i)
        if(fvalues[i] < fvalmin)
            nbest = i;
    for(int i = 0; i < dim; ++i)
        out[i] = population[nbest][i];
    for(int i = 0; i < N; ++i){
        delete[] population[i];
        delete[] children[i];
    }
    delete[] population;
    delete[] children;
    delete[] parents;
    delete[] fvalues;
}

/*
void GSA::optimizeSA(double f(double*,int), double *init, int dim, int numch, int bins, double dt, double T, double *params, double *vector, double *out, double &minf){
    //define params
    int N = 1;
    int Ns = params[0];
    double To = params[1];
    int Nt = params[2];
    double alpha = params[3];
    double c = params[4];
    int Ne = params[5];
    omprng MyRng;
    //init population
    double ** population = new double*[1];
    population[0] = new double[dim];
    for(int i = 0; i < dim; ++i)
        population[0][i] = init[i];
    //function values
    double* fs = new double[1];
    double* fsch = new double[1];
    fit.GSAiteration(population,N,numch,bins,dt,T,fs,1,&MyRng);
    cout << "Init : " << endl;
    for(int i = 0; i < dim; ++i)
        cout << population[0][i] << endl;
    cout << fs[0] << " ";
    cout << endl;
    double** children = new double*[1];
    children[0] = new double[dim];
    //sobsna, cikl
    int Ncur = 0;
    double ff = fs[0];
    while(Ncur < Ne){
        int Nacc = 0;
        for(int i = 1; i <= Ns; ++i){
            for(int j = 1; j <= Nt; ++j){
                //generate new vector
                for(int z = 0; z < dim; ++z){
                    if(z == 2){
                        double k = children[0][0];
                        double a = children[0][1];
                        children[0][z] = MyRng.runif((6.6+0.5*a/k)*a,8.65*a);
                        continue;
                    }
                    if(z == 9){
                        double k = children[0][0];
                        children[0][z] = MyRng.runif(0,2*sqrt(k));
                        continue;
                    }
                    double r = MyRng.runif(-1,1);
                    children[0][z] = population[0][z] + r*vector[z];
                    while(children[0][z] <= 0){
                        r = MyRng.runif(-1,1);
                        children[0][z] = population[0][z] + r*vector[z];
                    }
                }
                //Boltzmann
                fit.GSAiteration(children,N,numch,bins,dt,T,fsch,1,&MyRng);
                double r = MyRng.runif();
                if(exp((fs[0]-fsch[0])/To) > r){
                    ++Nacc;
                    for(int z = 0; z < dim; ++z)
                        population[0][z] = children[0][z];
                    fs[0] = fsch[0];
                    }
                for(int z = 0; z < dim; ++z)
                    cout << population[0][z] << endl;
                cout << "function:" << fs[0] << endl;
                cout << "--------------------------" << endl;
            }
            To*=alpha;
        }
        //change vector
        double a = 2*(1-2*c + c*c)/c, b = 4 - 3/c - c, d = 1/c;
        double n = double(Nacc)/(Ns*Nt);
        double g = a*n*n + b*n + d;
        for(int z = 0; z < dim; ++z)
            vector[z] *= g;
        if(abs(ff - fs[0]) < 1e-3)
            ++Ncur;
        ff = fs[0];
    }
    minf = ff;
    for(int i = 0; i < dim; ++i)
        out[i] = population[0][i];
    delete[] population[0];
    delete[] children[0];
    delete[] population;
    delete[] children;
    delete[] fs;
    delete[] fsch;
}
*/
