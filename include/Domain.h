#ifndef DOMAIN_H
#define DOMAIN_H

#include "Cell.h"


class Domain
{
private:
    int m_nbCells;

    std::vector<Cell*> m_cells;

    Parameters* m_param;

public:
    Domain(int nbCells);
    ~Domain();

    std::map<int, double> rightFlux(int i);

    std::map<int, double> leftFlux(int i);

    std::map<std::pair<int,int>,double> flux(int i);

    std::map<std::pair<int,int>,double> RHS(int i);


};















#endif