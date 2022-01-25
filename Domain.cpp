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
        for(int j=0;j<m_param->nbVar;j++)
        {
            m_EigenLeft[std::make_pair(i,j)]={0.0,0.0,0.0};
            m_EigenRight[std::make_pair(i,j)]={0.0,0.0,0.0};
            
            m_invEigenRight[std::make_pair(i,j)]={0.0,0.0,0.0};
            m_invEigenLeft[std::make_pair(i,j)]={0.0,0.0,0.0};
            m_rightEigen[std::make_pair(i,j)]={0.0,0.0,0.0};
            m_leftEigen[std::make_pair(i,j)]={0.0,0.0,0.0};
        }

        m_EigenValLeft[i]={0.0,0.0,0.0};
        m_EigenValRight[i]={0.0,0.0,0.0};        
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
           
            flux[0]=m_cells[i]->getFU(getLeftBorderMinmod(i))[0];
           
            flux[1]=m_cells[i]->getFU(getLeftBorderMinmod(i))[1];
            flux[2]=m_cells[i]->getFU(getLeftBorderMinmod(i))[2];
        }
        else
        {
            
            double lambda1=m_cells[i-1]->getRightEigen();
            double lambda2=m_cells[i]->getLeftEigen();
            double lambda=std::max(lambda1,lambda2);
            
            double fLeft=m_cells[i-1]->getFU(getRightBorderMinmod(i-1))[0];
            double fRight=m_cells[i]->getFU(getLeftBorderMinmod(i))[0];
          
            double left=getRightBorderMinmod(i-1)[0];
            double right=getLeftBorderMinmod(i)[0];
            
            flux[0]=0.5*(fLeft+fRight-lambda*(right-left));
            

            fLeft=m_cells[i-1]->getFU(getRightBorderMinmod(i-1))[1];
            fRight=m_cells[i]->getFU(getLeftBorderMinmod(i))[1];
            left=getRightBorderMinmod(i-1)[1];
            right=getLeftBorderMinmod(i)[1];
            flux[1]=0.5*(fLeft+fRight-lambda*(right-left));

            fLeft=m_cells[i-1]->getFU(getRightBorderMinmod(i-1))[2];
            fRight=m_cells[i]->getFU(getLeftBorderMinmod(i))[2];
            left=getRightBorderMinmod(i-1)[2];
            right=getLeftBorderMinmod(i)[2];
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
            
            flux[0]=m_cells[i]->getFU(getRightBorderMinmod(i))[0];
            
            flux[1]=m_cells[i]->getFU(getRightBorderMinmod(i))[1];
            flux[2]=m_cells[i]->getFU(getRightBorderMinmod(i))[2];
        }
        else
        {
            double lambda1=m_cells[i]->getRightEigen();
            double lambda2=m_cells[i+1]->getLeftEigen();
            double lambda=std::max(lambda1,lambda2);
            
            double fLeft=m_cells[i]->getFU(getRightBorderMinmod(i))[0];
            double fRight=m_cells[i+1]->getFU(getLeftBorderMinmod(i+1))[0];
            double left=getRightBorderMinmod(i)[0];
            double right=getLeftBorderMinmod(i+1)[0];
            flux[0]=0.5*(fLeft+fRight-lambda*(right-left));

            fLeft=m_cells[i]->getFU(getRightBorderMinmod(i))[1];
            fRight=m_cells[i+1]->getFU(getLeftBorderMinmod(i+1))[1];
            left=getRightBorderMinmod(i)[1];
            right=getLeftBorderMinmod(i+1)[1];
            flux[1]=0.5*(fLeft+fRight-lambda*(right-left));

            fLeft=m_cells[i]->getFU(getRightBorderMinmod(i))[2];
            fRight=m_cells[i+1]->getFU(getLeftBorderMinmod(i+1))[2];
            left=getRightBorderMinmod(i)[2];
            right=getLeftBorderMinmod(i+1)[2];
            flux[2]=0.5*(fLeft+fRight-lambda*(right-left));
        }
        return flux;     
    }
};


std::map<std::pair<int,int>,double> Domain::flux(int l)
{
   
    std::map<std::pair<int,int>,double> flux;
    std::map<int, double> lFlux, rFlux;
    if(m_param->isCharac)
    {
        rFlux=rightFluxCharc(l);
        lFlux=leftFluxCharc(l);
    }
    else
    {
        rFlux=rightFlux(l);
        lFlux=leftFlux(l);
    }
    
    
    //std::map<int, double> 
    
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
                //std::cout<<"flux0= "<<flux[std::make_pair(i,j)]<<std::endl;
            }
            if(j==1)
            {
                flux[std::make_pair(i,j)]=rFlux[i]+lFlux[i];
              //  std::cout<<"flux1= "<<flux[std::make_pair(i,j)]<<std::endl;
            }
            if(j==2)
            {
                flux[std::make_pair(i,j)]=rFlux[i]-lFlux[i];
               // std::cout<<"flux2= "<<rFlux[i]<<std::endl;
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


void Domain::minmod()
{
    std::vector<double> Minmod(2,0.0);
    double left;
    double right;
    double deltaPlus;
    double deltaMinus;
    std::pair<int,int> key;
    for(int cell=0;cell<m_nbCells;cell++)
    {
        for(int j=0;j<m_param->nbVar;j++)
        {
            key=std::make_pair(cell,j);
            left= -m_cells[cell]->getLeftBorder()[j]+m_cells[cell]->getFreedom()[std::make_pair(j,0)];
            
            right=m_cells[cell]->getRightBorder()[j]-m_cells[cell]->getFreedom()[std::make_pair(j,0)];
            
            if(cell==0)
            {
            
                deltaPlus=m_cells[cell+1]->getFreedom()[std::make_pair(j,0)]-m_cells[cell]->getFreedom()[std::make_pair(j,0)];
                deltaMinus=0.0;
                m_leftBorderMod[key]=minmod2(left,deltaPlus);
                
                m_rightBorderMod[key]=minmod2(right,deltaPlus);
            }
            else if(cell==m_nbCells-1)
            {
                
                deltaPlus=0.0;
                deltaMinus=m_cells[cell]->getFreedom()[std::make_pair(j,0)]-m_cells[cell-1]->getFreedom()[std::make_pair(j,0)];

                m_leftBorderMod[key]=minmod2(left,deltaMinus);
                
                m_rightBorderMod[key]=minmod2(right,deltaMinus);
            
            }
            else
            {
            
                deltaPlus=m_cells[cell+1]->getFreedom()[std::make_pair(j,0)]-m_cells[cell]->getFreedom()[std::make_pair(j,0)];

                deltaMinus=m_cells[cell]->getFreedom()[std::make_pair(j,0)]-m_cells[cell-1]->getFreedom()[std::make_pair(j,0)];

                m_leftBorderMod[key]=mimod(left,deltaPlus,deltaMinus);

                m_rightBorderMod[key]=mimod(right,deltaPlus,deltaMinus);                
            }
        }       
    }
};

double Domain::mimod(double a, double b, double c)
{
    double M=1.0;
    if(fabs(a)<=M*m_param->dx*m_param->dx)
    {
        return a;
    }
    else if(a*b>0 && a*c>0)
    {
        if(a>=0)
        {
            return std::min({fabs(a),fabs(b),fabs(c)});
        }
        else
        {
            return -1.0*std::min({fabs(a),fabs(b),fabs(c)});
        }        
    }
    else
    {
        //std::cout<<"enter 2"<<std::endl;
        return 0.0;
    }
  
};

double Domain::minmod2(double a, double b)
{
    double M=1.0;
    if(fabs(a)<=M*m_param->dx*m_param->dx)
    {
        return a;
    }
    else if(a*b>0)
    {
        if(a>=0)
        {
            return std::min({fabs(a),fabs(b)});
        }
        else
        {
            return -1.0*std::min({fabs(a),fabs(b)});
        }        
    }
    else
    {
        //std::cout<<"enter 2"<<std::endl;
        return 0.0;
    }
};

std::map<int,double> Domain::getLeftBorderMinmod(int cell)
{
    std::map<int,double> leftBorderMinmod;
    
    for(int i=0;i<m_param->nbVar;i++)
    {
        leftBorderMinmod[i]=m_cells[cell]->getFreedom()[std::make_pair(i,0)]-m_leftBorderMod[std::make_pair(cell,i)];//minmod(cell,i)[0];
    }
  
    return leftBorderMinmod;
}
std::map<int,double> Domain::getRightBorderMinmod(int cell)
{
    std::map<int,double> rightBorderMinmod;
 
    for(int i=0;i<m_param->nbVar;i++)
    {
        rightBorderMinmod[i]=m_cells[cell]->getFreedom()[std::make_pair(i,0)]+m_rightBorderMod[std::make_pair(cell,i)];///minmod(cell,i)[1];
     
    }
  
    return rightBorderMinmod;
}

void Domain::updateMod()
{
    double a;
    for(int cell=0;cell<m_nbCells;cell++)
    {
        for(int i=0;i<m_param->nbVar;i++)
        {
            if(m_param->Order==2)
            {
                a=(1.0/6.0)*m_leftBorderMod[std::make_pair(cell,i)];
                m_cells[cell]->updateMod(i,1,a);
            }
            if(m_param->Order==3)
            {
                a=(1.0/12.0)*(m_leftBorderMod[std::make_pair(cell,i)]+m_rightBorderMod[std::make_pair(cell,i)]);
                m_cells[cell]->updateMod(i,1,a);
                a=(1.0/60.0)*(m_rightBorderMod[std::make_pair(cell,i)]-m_leftBorderMod[std::make_pair(cell,i)]);
                m_cells[cell]->updateMod(i,2,a);
            }
        }
    }
}

