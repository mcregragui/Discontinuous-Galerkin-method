#ifndef PARAMETERS_H
#define PARAMETERS_H


#include "Template.h"
#include<map>



struct Parameters
{
    int xmin;

    int xmax;

    int nbVar;

    int RK;

    int Order;

    double dx;

    double gamma;

    int method;

    
};


struct Quadrature
{

    std::vector<double> Weights;

    std::vector<double>  Points;
    
};















#endif