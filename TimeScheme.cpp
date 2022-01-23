#include "include/TimeScheme.h"

TimeScheme::TimeScheme(Parameters* param,Quadrature* quad)
{
    m_param=param;
    m_end=m_param->end;
    m_nb=m_param->nbCells;
    m_domain=new Domain(param,quad);
};

TimeScheme::~TimeScheme()
{
   // delete m_domain;
};

void TimeScheme::advance()
{
    
    m_domain->Eigenvector();
    for(int l=0;l<m_nb;l++)
    {
        m_domain->getCells()[l]->borders();
    }
    m_domain->minmodCharc();
    m_domain->minmod();
    m_domain->updateCharc();
    for(int l=0;l<m_nb;l++)
    {
        m_domain->getCells()[l]->quadrature();
       // m_domain->getCells()[l]->borders();
        m_domain->getCells()[l]->eigens();
    }
    
    
    std::map<std::pair<int, int>, double> RHS;
    if(m_param->RK==1)
    {
        std::map<std::pair<int,int>,double> freedom0;
        std::vector<std::map<std::pair<int,int>,double>> freedom1;
        std::map<std::pair<int,int>,std::vector<double>> freedom2;
        for(int l=0;l<m_domain->getCells().size();l++)
        {
            RHS=m_domain->RHS(l);
            for(int i=0;i<m_param->nbVar;i++)
            {
                
                for(int j=0;j<m_param->Order;j++)
                {
                    
                    freedom0[std::make_pair(i,j)]=m_domain->getCells()[l]->getFreedom()[std::make_pair(i,j)]+m_param->dt*RHS[std::make_pair(i,j)];
                    //std::cout<<"pf3= "<<RHS[std::make_pair(i,j)]<<std::endl;
                }
            }
            m_domain->getCells()[l]->updateTemp(freedom0);
                //freedom1.push_back(freedom0);
            freedom0.clear();
            RHS.clear();
           
        } 
        for(int l=0;l<m_domain->getCells().size();l++)
        {
            //std::cout<<"pf3= "<<freedom1[l][std::make_pair(0,0)]<<std::endl;
            m_domain->getCells()[l]->updateFreedom();
        }

        freedom0.clear();

        freedom1.clear();
        
    }
    if(m_param->RK==2)
    {
        std::map<std::pair<int,int>,double> freedom;
        std::vector<std::map<std::pair<int,int>,double>> freedomb;
        std::map<std::pair<int,int>,double> freedom0;
        std::vector<std::map<std::pair<int,int>,double>> freedom1;
        std::map<std::pair<int,int>,double> freedom2;
        std::vector<std::map<std::pair<int,int>,double>> freedom3;

        for(int l=0;l<m_domain->getCells().size();l++)
        {
            for(int i=0;i<m_param->nbVar;i++)
            {
                
                for(int j=0;j<m_param->Order;j++)
                {                    
                    freedom[std::make_pair(i,j)]=m_domain->getCells()[l]->getFreedom()[std::make_pair(i,j)];
                }

                
                
            }
            freedomb.push_back(freedom);

            freedom.clear();
        }

        for(int l=0;l<m_domain->getCells().size();l++)
        {
            RHS=m_domain->RHS(l);
            for(int i=0;i<m_param->nbVar;i++)
            {
                
                for(int j=0;j<m_param->Order;j++)
                {
                    
                    freedom0[std::make_pair(i,j)]=m_domain->getCells()[l]->getFreedom()[std::make_pair(i,j)]+m_param->dt*RHS[std::make_pair(i,j)];
                }

                
                
            }
            freedom1.push_back(freedom0);

            freedom0.clear();
            RHS.clear();
        }

        for(int l=0;l<m_domain->getCells().size();l++)
        {
            m_domain->getCells()[l]->updateFree(freedom1[l]);
        }    

        for(int l=0;l<m_domain->getCells().size();l++)
        {
            RHS=m_domain->RHS(l);
            for(int i=0;i<m_param->nbVar;i++)
            {
                
                for(int j=0;j<m_param->Order;j++)
                {
                    
                    freedom2[std::make_pair(i,j)]=0.5*freedomb[l][std::make_pair(i,j)]+0.5*freedom1[l][std::make_pair(i,j)]+0.5*m_param->dt*RHS[std::make_pair(i,j)];
                }      
            }
            freedom3.push_back(freedom2);

            freedom0.clear();
        }

        for(int l=0;l<m_domain->getCells().size();l++)
        {
            m_domain->getCells()[l]->updateFree(freedom3[l]);
        }    
    }
    if(m_param->RK==3)
    {
        std::map<std::pair<int,int>,double> freedom;
        std::vector<std::map<std::pair<int,int>,double>> freedomb;
        std::map<std::pair<int,int>,double> freedom0;
        std::vector<std::map<std::pair<int,int>,double>> freedom1;
        std::map<std::pair<int,int>,double> freedom2;
        std::vector<std::map<std::pair<int,int>,double>> freedom3;
        std::map<std::pair<int,int>,double> freedom4;
        std::vector<std::map<std::pair<int,int>,double>> freedom5;

        for(int l=0;l<m_domain->getCells().size();l++)
        {
            for(int i=0;i<m_param->nbVar;i++)
            {
                
                for(int j=0;j<m_param->Order;j++)
                {                    
                    freedom[std::make_pair(i,j)]=m_domain->getCells()[l]->getFreedom()[std::make_pair(i,j)];
                }

                
                
            }
            freedomb.push_back(freedom);

            freedom.clear();
        }

        for(int l=0;l<m_domain->getCells().size();l++)
        {
            RHS=m_domain->RHS(l);
            for(int i=0;i<m_param->nbVar;i++)
            {
                
                for(int j=0;j<m_param->Order;j++)
                {
                    
                    freedom0[std::make_pair(i,j)]=m_domain->getCells()[l]->getFreedom()[std::make_pair(i,j)]+m_param->dt*RHS[std::make_pair(i,j)];
                }

                
                
            }
            freedom1.push_back(freedom0);

            freedom0.clear();
            RHS.clear();
        }

        for(int l=0;l<m_domain->getCells().size();l++)
        {
            m_domain->getCells()[l]->updateFree(freedom1[l]);
        }    

        for(int l=0;l<m_domain->getCells().size();l++)
        {
            RHS=m_domain->RHS(l);
            for(int i=0;i<m_param->nbVar;i++)
            {
                
                for(int j=0;j<m_param->Order;j++)
                {
                    
                    freedom2[std::make_pair(i,j)]=(3.0/4.0)*freedomb[l][std::make_pair(i,j)]+0.25*freedom1[l][std::make_pair(i,j)]+0.25*m_param->dt*RHS[std::make_pair(i,j)];
                }      
            }
            freedom3.push_back(freedom2);

            freedom2.clear();
        }

        for(int l=0;l<m_domain->getCells().size();l++)
        {
            m_domain->getCells()[l]->updateFree(freedom3[l]);
        }    
        for(int l=0;l<m_domain->getCells().size();l++)
        {
            RHS=m_domain->RHS(l);
            for(int i=0;i<m_param->nbVar;i++)
            {
                for(int j=0;j<m_param->Order;j++)
                {
                    freedom4[std::make_pair(i,j)]=(1.0/3.0)*freedomb[l][std::make_pair(i,j)]+(2.0/3.0)*freedom3[l][std::make_pair(i,j)]+(2.0/3.0)*m_param->dt*RHS[std::make_pair(i,j)];
                }
            }
            freedom5.push_back(freedom4);
            freedom4.clear();
        }
        for(int l=0;l<m_domain->getCells().size();l++)
        {
            m_domain->getCells()[l]->updateFree(freedom5[l]);
        }  
    }
};

void TimeScheme::advances()
{
    
   
    int n=0;
    while (n<m_param->end)
    {
        advance();
        n++;
        
    } 
    saveSol();
};

void TimeScheme::saveSol()
{
    std::ofstream myfile;
    myfile.open ("solMinmod.txt");

    double a=0.0;
    double rho=0.0;
    double u=0.0;
    double e=0.0;
    double p=0.0;
    std::map<int, double> sol;
    for(int l=0;l<m_domain->getCells().size();l++)
    {
        a=m_domain->getCells()[l]->getPos();
        if(m_param->nbVar==1)
        {
            myfile <<a<<" "<<m_domain->getCells()[l]->getSolution(a)[0]<<std::endl;;
        }
        if(m_param->nbVar==3)
        {
            sol=m_domain->getCells()[l]->getSolution(a);
            rho=sol[0];
            u=sol[1]/rho;
            e=sol[2];
            p=(m_param->gamma-1.0)*(e-0.5*rho*u*u);
           
            myfile <<a<<" "<<rho<<" "<<u<<" "<<e<<" "<<p<<std::endl;
        }
        
    }
    
    myfile.close();
};