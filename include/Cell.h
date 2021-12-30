#ifndef CELL_H
#define CELL_H

#include "Parameters.h"
#include <algorithm>
#include<cmath>




class Cell
{
private:
    int m_j;

    double m_dx;

    double m_rightPos;

    double m_leftPos;

    double m_pos;

    std::map<std::pair<int,int>,double> m_freedom;

    std::map<std::pair<int,int>,double> m_integral;

    Parameters* m_param;

    Quadrature* m_quad;

    std::map<int,double> m_leftBorder;

    std::map<int,double> m_rightBorder;

    double m_leftEigen;

    double m_rightEigen;

    double m_maxEigen;


public:
    Cell(int j, double dx,Parameters* param,Quadrature* quad);
    ~Cell();

    void positions();

    void initial();

    std::map<int,double> getInit();

    void quadrature();

    std::map<int,double> getSolution(double x);

    std::map<int,double> getFunction(double x);

    std::map<int,double> getFU(std::map<int,double> sol);

    void borders();

    double getLegendre(int order, double x);

    double getCoff(int order);

    void eigens();

    void updateFreedom(std::map<std::pair<int,int>,double> freedom){m_freedom=freedom;};

    std::map<int,double> getLeftBorder(){return m_leftBorder;};
    std::map<int,double> getRightBorder(){return m_rightBorder;};
    double getLeftEigen(){return m_leftEigen;};
    double getRightEigen(){return m_rightEigen;};
    std::map<std::pair<int,int>,double> getIntegral(){return m_integral;};
    double getdx(){return m_dx;};
    std::map<std::pair<int,int>,double> getFreedom(){return m_freedom;};
    double getPos(){return m_pos;};

    


};























#endif 