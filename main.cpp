#include "include/TimeScheme.h"



int main()
{

    Parameters param;

    param.adv=1.0;
    param.dt=0.01;
    param.dx=0.1;
    param.end=50;
    param.gamma=1.4;
    param.nbVar=1;
    param.Order=1;
    param.RK=1;
    param.xmin=-5.0;
    param.nbCells=100;
    

    Quadrature quad;

    
    quad.Points.push_back(sqrt((7-2*sqrt(7))/21.0));
    quad.Points.push_back(-sqrt((7-2*sqrt(7))/21.0));
    quad.Points.push_back(sqrt((7+2*sqrt(7))/21.0));
    quad.Points.push_back(-sqrt((7+2*sqrt(7))/21.0));
    quad.Points.push_back(1.0);
    quad.Points.push_back(-1.0);

    quad.Weights.push_back((14+sqrt(7))/30.0);
    quad.Weights.push_back((14+sqrt(7))/30.0);
    quad.Weights.push_back((14-sqrt(7))/30.0);
    quad.Weights.push_back((14-sqrt(7))/30.0);
    quad.Weights.push_back(1.0/15.0);
    quad.Weights.push_back(1.0/15.0);


    TimeScheme* time=new TimeScheme();
    
    time->advances();

    return 0;
}