#include "include/Domain.h"


Domain::Domain(Parameters* param,Quadrature* quad)
{
    m_param=param;
    m_nbCells=m_param->nbCells;
    for(int i=0;i<m_param->nbCells;i++)
    {
        if(i==0)
        {
            m_cells.push_back(new Cell(i,m_param->xmin,param,quad));
        }
        else
        {
            
            m_cells.push_back(new Cell(i,m_cells[i-1]->getRightPos(),param,quad));
        }
        
    }
};

Domain::~Domain()
{
   for(std::vector<Cell*>::iterator it=m_cells.begin();it!=m_cells.end();++it)
    {
        delete *it;
    }
};


std::map<int, double> Domain::leftFlux(int i)
{
    std::map<int, double> flux;

    if(m_param->nbVar==1)
    {
        if(i==0)
        {
            double fRight=m_cells[i]->getFU(m_cells[i]->getLeftBorder())[0];
            flux[0]=fRight;
        }
        else
        {
            double fLeft=m_cells[i-1]->getFU(m_cells[i-1]->getRightBorder())[0];
            double fRight=m_cells[i]->getFU(m_cells[i]->getLeftBorder())[0];
            double left=m_cells[i-1]->getRightBorder()[0];
            double right=m_cells[i]->getLeftBorder()[0];
            flux[0]=0.5*(fLeft+fRight-m_param->adv*(right-left));
        }
        return flux;      
    }
    if(m_param->nbVar==3)
    {
        if(i==0)
        {
            flux[0]=m_cells[i]->getFU(m_cells[i]->getLeftBorder())[0];
            flux[1]=m_cells[i]->getFU(m_cells[i]->getLeftBorder())[1];
            flux[2]=m_cells[i]->getFU(m_cells[i]->getLeftBorder())[2];
        }
        else
        {
            double lambda1=m_cells[i-1]->getRightEigen();
            double lambda2=m_cells[i]->getLeftEigen();
            double lambda=std::max(lambda1,lambda2);
           
            double fLeft=m_cells[i-1]->getFU(m_cells[i-1]->getRightBorder())[0];
            double fRight=m_cells[i]->getFU(m_cells[i]->getLeftBorder())[0];
          
            double left=m_cells[i-1]->getRightBorder()[0];
            double right=m_cells[i]->getLeftBorder()[0];
            
            flux[0]=0.5*(fLeft+fRight-lambda*(right-left));
            

            fLeft=m_cells[i-1]->getFU(m_cells[i-1]->getRightBorder())[1];
            fRight=m_cells[i]->getFU(m_cells[i]->getLeftBorder())[1];
            left=m_cells[i-1]->getRightBorder()[1];
            right=m_cells[i]->getLeftBorder()[1];
            flux[1]=0.5*(fLeft+fRight-lambda*(right-left));

            fLeft=m_cells[i-1]->getFU(m_cells[i-1]->getRightBorder())[2];
            fRight=m_cells[i]->getFU(m_cells[i]->getLeftBorder())[2];
            left=m_cells[i-1]->getRightBorder()[2];
            right=m_cells[i]->getLeftBorder()[2];
            flux[2]=0.5*(fLeft+fRight-lambda*(right-left));
        }
        return flux;     
    }
};


std::map<int, double> Domain::rightFlux(int i)
{
    std::map<int, double> flux;

    if(m_param->nbVar==1)
    {
        if(i==m_nbCells-1)
        {
            flux[0]=m_cells[i]->getRightBorder()[0];
        }
        else
        {
            double fLeft=m_cells[i]->getFU(m_cells[i]->getRightBorder())[0];
            double fRight=m_cells[i+1]->getFU(m_cells[i+1]->getLeftBorder())[0];
            double left=m_cells[i]->getRightBorder()[0];
            double right=m_cells[i+1]->getLeftBorder()[0];
            flux[0]=0.5*(fLeft+fRight-m_param->adv*(right-left));
        }
        return flux;      
    }
    if(m_param->nbVar==3)
    {
        if(i==m_nbCells-1)
        {
            flux[0]=m_cells[i]->getFU(m_cells[i]->getRightBorder())[0];
            flux[1]=m_cells[i]->getFU(m_cells[i]->getRightBorder())[1];
            flux[2]=m_cells[i]->getFU(m_cells[i]->getRightBorder())[2];
        }
        else
        {
            double lambda1=m_cells[i]->getRightEigen();
            double lambda2=m_cells[i+1]->getLeftEigen();
            double lambda=std::max(lambda1,lambda2);
          
            double fLeft=m_cells[i]->getFU(m_cells[i]->getRightBorder())[0];
            double fRight=m_cells[i+1]->getFU(m_cells[i+1]->getLeftBorder())[0];
            double left=m_cells[i]->getRightBorder()[0];
            double right=m_cells[i+1]->getLeftBorder()[0];
            flux[0]=0.5*(fLeft+fRight-lambda*(right-left));

            fLeft=m_cells[i]->getFU(m_cells[i]->getRightBorder())[1];
            fRight=m_cells[i+1]->getFU(m_cells[i+1]->getLeftBorder())[1];
            left=m_cells[i]->getRightBorder()[1];
            right=m_cells[i+1]->getLeftBorder()[1];
            flux[1]=0.5*(fLeft+fRight-lambda*(right-left));

            fLeft=m_cells[i]->getFU(m_cells[i]->getRightBorder())[2];
            fRight=m_cells[i+1]->getFU(m_cells[i+1]->getLeftBorder())[2];
            left=m_cells[i]->getRightBorder()[2];
            right=m_cells[i+1]->getLeftBorder()[2];
            flux[2]=0.5*(fLeft+fRight-lambda*(right-left));
        }
        return flux;     
    }
};


std::map<std::pair<int,int>,double> Domain::flux(int l)
{
    std::map<std::pair<int,int>,double> flux;
    
    std::map<int, double> rFlux=rightFlux(l);
    
    std::map<int, double> lFlux=leftFlux(l);
    for(int i=0;i<m_param->nbVar;i++)
    {
        
        for(int j=0;j<m_param->Order;j++)
        {
             flux[std::make_pair(i,j)]=0.0;
        }
    }


    for(int i=0;i<m_param->nbVar;i++)
    {
        
        for(int j=0;j<m_param->Order;j++)
        {
            
            if(j==0)
            {
                
                flux[std::make_pair(i,j)]=rFlux[i]-lFlux[i];
              
            }
            if(j==1)
            {
                flux[std::make_pair(i,j)]=rFlux[i]+lFlux[i];
            }
            if(j==2)
            {
                flux[std::make_pair(i,j)]=rFlux[i]-lFlux[i];
            }       
        }
    }
    
    return flux;

};

std::map<std::pair<int,int>,double> Domain::RHS(int l)
{
   
    std::map<std::pair<int,int>,double> rhs;
    
    std::map<std::pair<int,int>,double> f=flux(l);
    
   
    for(int i=0;i<m_param->nbVar;i++)
    {
        
        for(int j=0;j<m_param->Order;j++)
        {
            rhs[std::make_pair(i,j)]=0.0;
        }
    }

    double dx=m_cells[l]->getdx();
    
    for(int i=0;i<m_param->nbVar;i++)
    {
        
        for(int j=0;j<m_param->Order;j++)
        {
            if(j==0)
            {
                rhs[std::make_pair(i,j)]=-(1.0/dx)*f[std::make_pair(i,j)];
               
            }
            if(j==1)
            {
                rhs[std::make_pair(i,j)]=-(0.5/dx)*f[std::make_pair(i,j)]+(1.0/(dx*dx))*m_cells[l]->getIntegral()[std::make_pair(i,j)];
            }
            if(j==2)
            {
                rhs[std::make_pair(i,j)]=-(1.0/(6.0*dx))*f[std::make_pair(i,j)]+(2.0/(dx*dx*dx))*m_cells[l]->getIntegral()[std::make_pair(i,j)];
            }
            
        }
    }
    return rhs;
};


void Domain::minmod(int cell)
{
    if(cell==0)
    {
        
    }
};