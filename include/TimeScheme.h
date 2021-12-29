#ifndef TIME_SCHEME_H
#define TIME_SCHEME_H

#include "Domain.h"
#include<fstream>
#include<string>

class TimeScheme
{
private:
    int m_end;

    int m_nb;
    Parameters* m_param;
    Domain* m_domain;
public:
    TimeScheme();
    ~TimeScheme();

    void advance();

    void advances();

    void saveSol();
};





#endif