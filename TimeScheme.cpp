#include "include/TimeScheme.h"

TimeScheme::TimeScheme()
{
    m_end=m_param->end;
    m_nb=m_param->nbCells;
    m_domain=new Domain();
};

TimeScheme::~TimeScheme()
{
    delete m_domain;
};

void TimeScheme::advance()
{

    std::map<std::pair<int,int>,double> freedom0;
    std::map<std::pair<int,int>,double> freedom1;
    std::map<std::pair<int,int>,double> freedom2;

    if(m_param->RK==1)
    {
        
        for(int l=0;l<m_domain->getCells().size();l++)
        {
            for(int i=0;i<m_param->nbVar;i++)
            {
                for(int j=0;j<m_param->Order;j++)
                {
                    freedom0[std::make_pair(i,j)]=m_domain->getCells()[l]->getFreedom()[std::make_pair(i,j)]+m_param->dt*m_domain->RHS(l)[std::make_pair(i,j)];
                }
            }
            m_domain->getCells()[l]->updateFreedom(freedom0);
        } 
    }
    /*if(m_param->RK==2)
    {
        for(int l=0;l<m_domain->getCells().size();l++)
        {
            for(int i=0;i<m_param->nbVar;i++)
            {
                for(int j=0;j<m_param->Order;j++)
                {
                    freedom0[std::make_pair(i,j)]=m_domain->getCells()[l]->getFreedom()[std::make_pair(i,j)]+m_param->dt*m_domain->RHS(l)[std::make_pair(i,j)];
                }
            }
            m_domain->getCells()[l]->updateFreedom(freedom0);

            for(int i=0;i<m_param->nbVar;i++)
            {
                for(int j=0;j<m_param->Order;j++)
                {
                    freedom1[std::make_pair(i,j)]=m_domain->getCells()[l]->getFreedom()[std::make_pair(i,j)]+m_param->dt*m_domain->RHS(l)[std::make_pair(i,j)];
                }
            }
        }
    }*/
}

void TimeScheme::advances()
{
    int n=0;
    while (n<m_param->end)
    {
        advance();
        n++;
    } 
    saveSol();
}

void TimeScheme::saveSol()
{
    std::ofstream myfile;
    myfile.open ("sol.txt");

    double a=0.0;
    double rho=0.0;
    double u=0.0;
    double e=0.0;
    double p=0.0;
    for(int l=0;l<m_domain->getCells().size();l++)
    {
        a=m_domain->getCells()[l]->getPos();
        if(m_param->nbVar==1)
        {
            myfile <<a<<" "<<m_domain->getCells()[l]->getSolution(a)[0]<<std::endl;;
        }
        if(m_param->nbVar==3)
        {
            rho=m_domain->getCells()[l]->getSolution(a)[0];
            u=m_domain->getCells()[l]->getSolution(a)[1]/rho;
            e=m_domain->getCells()[l]->getSolution(a)[2];
            p=(m_param->gamma-1.0)*sqrt(e-0.5*rho*u*u);
            myfile <<a<<" "<<rho<<" "<<u<<" "<<e<<" "<<p<<std::endl;
        }
        
    }
    
    myfile.close();
}