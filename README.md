# GSA
A C++ implementation of GSA(Genetic simulated annealing) global optimization algorithm. Implementation consists of GSA class.
Usage:
1. Set GSA optimizer parameters:N,tmax,To,pmo,alpha,beta,K,Nitermax (see paper for more details)
2. Run optimize method. Parameters: double function(double*,int), number of threads, dimension of  argument vector, minimal values, maximal values, output array.
For random number class uses omprng: http://homepage.divms.uiowa.edu/~mbognar/omprng/
See example(GSA.cpp) for more details.

Reference:
Qiaoling Xu, Gongwang Zhang, Chao Zhao and Aimin An, A Robust Adaptive Hybrid  Genetic  Simulated  Annealing Algorithm forthe  global optimization  of multimodal functions, 2011 Chinese Control and Decision Conference (CCDC), doi:10.1109/CCDC.2011.5968132
