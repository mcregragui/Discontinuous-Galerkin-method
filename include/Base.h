#ifndef BASE_H
#define BASE_H

#include "Parameters.h"

class Base
{
private:
    
public:
    Base(/* args */);
    ~Base();

    double getLegendre(int order, double x);

    double getCoeff(int order);
};

Base::Base(/* args */)
{
}

Base::~Base()
{
}












#endif