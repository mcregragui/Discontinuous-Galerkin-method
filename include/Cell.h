#ifndef CELL_H
#define CELL_H

#include "Parameters.h"



class Cell
{
private:
    int m_j;

    double m_dx;

    double m_rightPos;

    double m_leftPos;

    double m_pos;

    std::map<std::pair<int,int>,double> m_freedom;

    std::map<int,double> m_integral;

    Parameters* m_param;

    Quadrature* m_quad;

    std::map<int,double> m_leftBorder;

    std::map<int,double> m_rightBorder;



public:
    Cell(/* args */);
    ~Cell();

    void positions();

    void quadrature();

    std::map<int,double> getSolution(double x);

    std::map<int,double> getFunction(double x);

    void borders();

    double getLegendre(int order, double x);

    double getCoff(int order);


};

Cell::Cell(/* args */)
{
}

Cell::~Cell()
{
}






















#endif 