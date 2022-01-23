#ifndef DOMAIN_H
#define DOMAIN_H

#include "Cell.h"


class Domain
{
private:
    int m_nbCells;

    std::vector<Cell*> m_cells;

    std::map<std::pair<int,int>, double> m_rightBorderMod;

    std::map<std::pair<int,int>, double> m_leftBorderMod;

    std::map<std::pair<int,int>, double> m_rightBorderCharc;

    std::map<std::pair<int,int>, double> m_leftBorderCharc;

    std::map<std::pair<int,int>, double> m_rightBorderTrue;

    std::map<std::pair<int,int>, double> m_leftBorderTrue;

    
    std::map<std::pair<int,int>, double> m_freeCharc;

    std::map<std::pair<int,int>, double> m_freeLeftCharc;
    std::map<std::pair<int,int>, double> m_freeRightCharc;
    std::map<std::pair<int,int>,std::vector<double>> m_leftEigen;

    
    std::map<std::pair<int,int>,double> m_leftEigenVal;

    std::map<std::pair<int,int>,std::vector<double>> m_EigenLeft;
    std::map<std::pair<int,int>,std::vector<double>> m_EigenRight;
    std::map<int,std::vector<double>> m_EigenValLeft;
    std::map<int,std::vector<double>> m_EigenValRight;
    std::map<std::pair<int,int>,std::vector<double>> m_invEigenRight;
    std::map<std::pair<int,int>,std::vector<double>> m_invEigenLeft;
    std::map<std::pair<int,int>,std::vector<double>> m_rightEigen;
    
    std::map<std::pair<int,int>,double> m_rightEigenVal;

    Parameters* m_param;

public:
    Domain(Parameters* param,Quadrature* quad);
    ~Domain();

    std::map<int, double> rightFlux(int i);

    std::map<int, double> leftFlux(int i);

    std::map<std::pair<int,int>,double> flux(int i);

    std::map<std::pair<int,int>,double> RHS(int i);

    std::vector<Cell*> getCells(){return m_cells;};

    void minmod();

    void minmodCharc();

    double mimod(double a, double b, double c);

    std::map<int,double> getLeftBorderMinmod(int cell);

    std::map<int,double> getRightBorderMinmod(int cell);

    void updateMod();

    //Characteristic projection

    void Eigenvector();

   // void rEigenvector();

    double Project(std::map<int,double>, int p);

    std::map<int,double> getLeftBorderCharc(int cell);

    std::map<int,double> getRightBorderCharc(int cell);

    std::map<int, double> rightFluxCharc(int i);

    std::map<int, double> leftFluxCharc(int i);

    void updateCharc();

    double project(std::vector<double> eigen, std::vector<double> v);
    double projectLeft(int cell, int i, std::vector<double> eigen, std::vector<double> v);

};















#endif