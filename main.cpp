#include "include/TimeScheme.h"



int main()
{
    Parameters* param=new Parameters();

    param->adv=1.0;
    param->dt=0.0001;
    param->dx=0.001;
    param->end=50;
    param->gamma=1.4;
    param->nbVar=1;
    param->Order=3;
    param->RK=1;
    param->xmin=-5.0;
    param->nbCells=10000;
    

    Quadrature* quad=new Quadrature();

    
    quad->Points.push_back(sqrt((7-2*sqrt(7))/21.0));
    quad->Points.push_back(-sqrt((7-2*sqrt(7))/21.0));
    quad->Points.push_back(sqrt((7+2*sqrt(7))/21.0));
    quad->Points.push_back(-sqrt((7+2*sqrt(7))/21.0));
    quad->Points.push_back(1.0);
    quad->Points.push_back(-1.0);

    quad->Weights.push_back((14+sqrt(7))/30.0);
    quad->Weights.push_back((14+sqrt(7))/30.0);
    quad->Weights.push_back((14-sqrt(7))/30.0);
    quad->Weights.push_back((14-sqrt(7))/30.0);
    quad->Weights.push_back(1.0/15.0);
    quad->Weights.push_back(1.0/15.0);

    TimeScheme* time=new TimeScheme(param,quad);
    
    time->advances();
    return 0;
}