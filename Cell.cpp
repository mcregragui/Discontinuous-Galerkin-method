#include<include/Cell.h>
#include<cassert>




//---------------------------CONSTRUCTOR-------------------------------------------------------------------------

Cell::Cell(int j, double dx,Parameters* param,Quadrature* quad)
{
    m_param=param;
    m_quad=quad;
    m_j=j;
    m_leftPos=dx;
    m_dx=m_param->dx;
    m_freedom.clear();
    m_temp.clear();
    positions();

    initial();

}

Cell::~Cell()
{
   delete m_param;
};


//----------------------------------------------------------------------------------------------------------



//-------------------------------------------DETERMINE POSITIONS: LEFT, RIGHT, MIDPOINTS-------------------------
void Cell::positions()
{
 
    m_rightPos=m_leftPos+m_dx;
    m_pos=(m_rightPos+m_leftPos)/2.0;
};
//---------------------------------------------------------------------------------------------------------------

//------------------------------------------INITIAL CONDITION IN THIS CELL---------------------------------------
void Cell::initial()
{
    assert(m_quad->Points.size()==m_quad->Weights.size());
    double h=(m_rightPos-m_leftPos)/2.0;
  
    for(int i=0;i<m_param->nbVar;i++)
    {
        for(int j=0;j<m_param->Order;j++)
        {
            m_freedom[std::make_pair(i,j)]=0.0;
        }
    }

    for(int i=0;i<m_param->nbVar;i++)
    {
        for(int j=0;j<m_param->Order;j++)
        {
            for(int l=0;l<m_quad->Points.size();l++)
            {
                m_freedom[std::make_pair(i,j)]+=m_quad->Weights[l]*getInit()[i]*getLegendre(j,m_pos+h*m_quad->Points[l]);
            }
            m_freedom[std::make_pair(i,j)]=(h/pow(m_dx,j+1))*m_freedom[std::make_pair(i,j)];
          
        }
    }
};

std::map<int,double> Cell::getInit()
{
    std::map<int,double> init;
    if(m_param->nbVar==1)
    {
        //init[0]=cos(m_pos);
        if(m_rightPos<=0.0)
        {
            init[0]=0.5;
        }
        else if(0<m_leftPos&&m_rightPos<2.0)
        {
            init[0]=3;
        }
        else
        {
            init[0]=0.5;
        }
        return init;
    }
    if(m_param->nbVar==3)
    {
        if(m_rightPos<=0.0)
        {
            init[0]=1.0;
            init[1]=0.0;
            init[2]=1.0/(m_param->gamma-1.0);

        }
        else
        {
            init[0]=0.125;
            init[1]=0.0;
            init[2]=0.1/(m_param->gamma-1.0);
        }

        return init;
    }
};
//---------------------------------------------------------------------------------------------------------------

//-------------------------------------------QUADRATURE-----------------------------------------------------------
void Cell::quadrature()
{
    if(m_param->Order==1)
    {
        for(int i=0;i<m_param->nbVar;i++)
        {
            for(int j=0;j<m_param->Order;j++)
            {
                m_integral[std::make_pair(i,j)]=0.0;
            }
            
        }
        
    }
    if(m_param->Order==2)
    {
        assert(m_quad->Points.size()==m_quad->Weights.size());

        for(int i=0;i<m_param->nbVar;i++)
        {
            for(int j=0;j<m_param->Order;j++)
            {
                m_integral[std::make_pair(i,j)]=0.0;
            }
            
        }

        double h=(m_rightPos-m_leftPos)/2.0;

        for(int i=0;i<m_param->nbVar;i++)
        {
            
            for(int j=0;j<m_quad->Points.size();j++)
            {
                m_integral[std::make_pair(i,1)]+=m_quad->Weights[j]*getFunction(m_pos+h*m_quad->Points[j])[i];
            }
            m_integral[std::make_pair(i,1)]=h*m_integral[std::make_pair(i,1)];
        
            
        }
    }

    if(m_param->Order==3)
    {
        assert(m_quad->Points.size()==m_quad->Weights.size());

        for(int i=0;i<m_param->nbVar;i++)
        {
            for(int l=0;l<m_param->Order;l++)
            {
                m_integral[std::make_pair(i,l)]=0.0;
            }
        }

        double h=(m_rightPos-m_leftPos)/2.0;

        for(int i=0;i<m_param->nbVar;i++)
        {
            for(int l=0;l<m_param->Order;l++)
            {
                if(l==1)
                {
                    for(int j=0;j<m_quad->Points.size();j++)
                    {
                        m_integral[std::make_pair(i,l)]+=m_quad->Weights[j]*getFunction(m_pos+h*m_quad->Points[j])[i];
                    }
                    m_integral[std::make_pair(i,l)]=h*m_integral[std::make_pair(i,l)];
                    
                }
                if(l==2)
                {
                    for(int j=0;j<m_quad->Points.size();j++)
                    {
                        m_integral[std::make_pair(i,l)]+=m_quad->Weights[j]*getFunction(m_pos+h*m_quad->Points[j])[i]*getLegendre(1,m_pos+h*m_quad->Points[j]);
                    }
                    m_integral[std::make_pair(i,l)]=h*m_integral[std::make_pair(i,l)];
                  
                }
            
            }
        }
    }
};

//-----------------------------------------------------------------------------------------------------------------


//------------------------------------------------------COMPUTE SOLUTION FROM FREEDOM DEGREES-----------------------

std::map<int,double> Cell::getSolution(double x)
{
    std::map<int,double> sol;

    for(int i=0;i<m_param->nbVar;i++)
    {
        sol[i]=0.0;
    }


    for(int i=0;i<m_param->nbVar;i++)
    {
        for(int j=0;j<m_param->Order;j++)
        {
            sol[i]=sol[i]+getCoff(j)*m_freedom[std::make_pair(i,j)]*getLegendre(j,x);
            
        }
        
    }
 
    return sol;
};

//-----------------------------------------------------------------------------------------------------------------


//--------------------------------------------------------COMPUTE FUNCTION FROM SOLUTION F(U)-----------------------

std::map<int,double> Cell::getFunction(double x)
{
    std::map<int,double> sol0=getSolution(x);
    std::map<int,double> sol;
    double p;
    
    if(m_param->nbVar==1)
    {
        sol[0]=m_param->adv*sol0[0];
        return sol;
    }
    if(m_param->nbVar==3)
    {
        p=(m_param->gamma-1.0)*(sol0[2]-0.5*sol0[1]*sol0[1]/sol0[0]);

        sol[0]=sol0[1];
        sol[1]=sol0[1]*sol0[1]/sol0[0]+p;

        sol[2]=sol0[1]*(sol0[2]+p)/sol0[0];

        return sol;
    }
};


std::map<int,double> Cell::getFU(std::map<int,double> sol0)
{
    std::map<int,double> sol;
    double p;
    switch (m_param->nbVar)
    {
    case 1:
        sol[0]=m_param->adv*sol0[0];
        return sol;

    case 3:
        p=(m_param->gamma-1.0)*(sol0[2]-0.5*sol0[1]*sol0[1]/sol0[0]);

        sol[0]=sol0[1];
        //std::cout<<"sol01= "<<sol[0]<<std::endl;
        sol[1]=sol0[1]*sol0[1]/sol0[0]+p;

        sol[2]=sol0[1]*(sol0[2]+p)/sol0[0];

        return sol;
            
    default:
        std::cout<<"No other function is defined"<<std::endl;
        break;
    }
};
//---------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------GET LEGENDRE POLYNOMES----------------------------------

double Cell::getLegendre(int order, double x)
{

    if(order==0)
    {
        return 1.0;
    }
    if(order==1)
    {
        return x-m_pos;
    }
    if(order==2)
    {
        return (x-m_pos)*(x-m_pos)-(1.0/12.0)*m_dx*m_dx;
    }
};
//------------------------------------------------------------------------------------------------------------

//---------------------------------------------------GET COEFFICIENT OF THE BASE-------------------------------

double Cell::getCoff(int order)
{
    if(order==0)
    {
        return 1.0;
    }
    if(order==1)
    {
        return 12.0/m_dx;
    }
    if(order==2)
    {
        return 180.0/(m_dx*m_dx);
    }
};

//-------------------------------------------------------------------------------------------------------------

//----------------------------------------CELL'S BORDERS TREATEMENT---------------------------------------------
void Cell::borders()
{

    //std::map<int,double> leftSol=getSolution();

    for(int i=0;i<m_param->nbVar;i++)
    {
        m_leftBorder[i]=0.0;
        m_rightBorder[i]=0.0;
    }
    if(m_param->Order==1)
    {
        for(int i=0;i<m_param->nbVar;i++)
        {
            m_leftBorder[i]=m_freedom[std::make_pair(i,0)];
            m_rightBorder[i]=m_freedom[std::make_pair(i,0)];
        }
    }
    if(m_param->Order==2)
    {
        for(int i=0;i<m_param->nbVar;i++)
        {
            m_leftBorder[i]=m_freedom[std::make_pair(i,0)]-6.0*m_freedom[std::make_pair(i,1)];
            m_rightBorder[i]=m_freedom[std::make_pair(i,0)]+6.0*m_freedom[std::make_pair(i,1)];
        }
    }
    if(m_param->Order==3)
    {
        for(int i=0;i<m_param->nbVar;i++)
        {
            m_leftBorder[i]=m_freedom[std::make_pair(i,0)]-6.0*m_freedom[std::make_pair(i,1)]+30.0*m_freedom[std::make_pair(i,2)];
            m_rightBorder[i]=m_freedom[std::make_pair(i,0)]+6.0*m_freedom[std::make_pair(i,1)]+30.0*m_freedom[std::make_pair(i,2)];
        }
    }
    
    
};
//---------------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------EIGENS--------------------------------------------------------


void Cell::eigens()
{
    if(m_param->nbVar==1)
    {
        m_leftEigen=m_param->adv;
        m_rightEigen=m_param->adv;
        m_maxEigen=m_param->adv;
    }
    if(m_param->nbVar==3)
    {
        m_leftEigen=0.0;
        m_rightEigen=0.0;
        m_maxEigen=0.0;


        double u=0.0;
        double c=0.0;
        double p=0.0;
        double max=0.0;
        //--------------------------------------------------Left----------------------
        std::map<int,double> sol0=getSolution(m_leftPos);
        
        u=sol0[1]/sol0[0];
        
        p=fabs((m_param->gamma-1.0)*(sol0[2]-0.5*sol0[1]*sol0[1]/sol0[0]));
        
        c=sqrt(fabs(m_param->gamma*p/sol0[0]));
        
       // std::cout<<"c= "<<c<<std::endl;
        max=std::max(fabs(u),fabs(u-c));
        m_leftEigen=std::max(max,fabs(u+c));
        
        //--------------------------------------------------Right-------------------------

        sol0=getSolution(m_rightPos);
        
        u=sol0[1]/sol0[0];
      
        p=((m_param->gamma-1.0)*(sol0[2]-0.5*sol0[1]*sol0[1]/sol0[0]));
      
        c=sqrt(fabs(m_param->gamma*p/sol0[0]));
       
        max=std::max(fabs(u),fabs(u-c));
        m_rightEigen=std::max(max,fabs(u+c));
        
        //---------------------------------------------------Max eigen---------------------------------

        m_maxEigen=std::max(fabs(m_rightEigen),fabs(m_leftEigen));
        
    }
};

void Cell::updateFreedom()
{
    m_freedom=m_temp;

    m_temp.clear();
};
