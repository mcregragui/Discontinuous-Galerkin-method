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

    double mimod(double a, double b, double c);

    std::map<int,double> getLeftBorderMinmod(int cell);

    std::map<int,double> getRightBorderMinmod(int cell);

    void updateMod();

    //Characteristic projection

    std::vector<double> lEigenvector(int i);

    std::vector<double> rEigenvector(int i);

    double Project(std::map<int,double>, int p);



};















#endif