#include<include/Cell.h>
#include<cassert>



//-------------------------------------------DETERMINE POSITIONS: LEFT, RIGHT, MIDPOINTS-------------------------
void Cell::positions()
{
    m_rightPos=m_param->xmin+m_j*m_dx;
    m_leftPos=m_rightPos+m_dx;
    m_pos=(m_rightPos+m_leftPos)/2.0;
};
//---------------------------------------------------------------------------------------------------------------

//-------------------------------------------QUADRATURE-----------------------------------------------------------
void Cell::quadrature()
{
    if(m_param->Order==1)
    {
        for(int i=0;i<m_param->nbVar;i++)
        {
            m_integral[i]=0.0;
        }
        
    }
    if(m_param->Order==2)
    {
        assert(m_quad->Points.size()==m_quad->Weights.size());

        for(int i=0;i<m_param->nbVar;i++)
        {
            m_integral[i]=0.0;
        }

        double h=(m_rightPos-m_leftPos)/2.0;

        for(int i=0;i<m_param->nbVar;i++)
        {
            for(int j=0;j<m_quad->Points.size();j++)
            {
                m_integral[i]+=m_quad->Points[j]*getFunction(m_pos+h*m_quad->Weights[j])[i];
            }
            m_integral[i]=h*m_integral[i];
        }
    }

    if(m_param->Order==3)
    {
        assert(m_quad->Points.size()==m_quad->Weights.size());

        for(int i=0;i<m_param->nbVar;i++)
        {
            m_integral[i]=0.0;
        }

        double h=(m_rightPos-m_leftPos)/2.0;

        for(int i=0;i<m_param->nbVar;i++)
        {
            for(int j=0;j<m_quad->Points.size();j++)
            {
                m_integral[i]+=m_quad->Points[j]*getFunction(m_pos+h*m_quad->Weights[j])[i]*getLegendre(1,m_pos+h*m_quad->Weights[j]);
            }
            m_integral[i]=h*m_integral[i];
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
            sol[i]=+getCoff(j)*m_freedom[std::make_pair(i,j)]*getLegendre(j,x);
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

    switch (m_param->nbVar)
    {
    case 1:
        return sol0;

    case 3:
        double p=0.0;

        p=(m_param->gamma-1.0)*(sol0[2]-0.5*sol0[1]*sol0[1]/sol0[0]);

        sol[0]=sol0[1];

        sol[1]=sol0[1]*sol0[1]/sol0[0]+p;

        sol[2]=sol0[1]*(sol0[2]+p)/sol0[0];
            
    default:
        std::cout<<"No other function is defined"<<std::endl;
        break;
    }
};
//---------------------------------------------------------------------------------------------------------------

//--------------------------------------------------------GET LEGENDRE POLYNOMES----------------------------------

double Cell::getLegendre(int order, double x)
{
    switch (order)
    {
    case 0:
        return 1.0;
    
    case 1:
        return x-m_pos;

    case 2:
        return (x-m_pos)*(x-m_pos)-(1.0/12.0)*m_dx*m_dx;
    
    default:
        std::cout<<"higher ordre is not available, you want to compute the ordre: "<<order<<std::endl; 
        break;
    }
};
//------------------------------------------------------------------------------------------------------------

//---------------------------------------------------GET COEFFICIENT OF THE BASE-------------------------------

double Cell::getCoff(int order)
{
    switch (order)
    {
    case 0:
        return 1.0;
    
    case 1:
        return 12.0/m_dx;

    case 2:
        return 180.0/(m_dx*m_dx);
    
    default:
        std::cout<<"higher ordre is not available, you want to compute the ordre: "<<order<<std::endl;
        break;
    }
};

//-------------------------------------------------------------------------------------------------------------


void Cell::borders()
{

    //std::map<int,double> leftSol=getSolution();

    for(int i=0;i<m_param->nbVar;i++)
    {
        m_leftBorder[i]=0.0;
        m_rightBorder[i]=0.0;
    }
    switch (m_param->Order)
    {
    case 1:
        for(int i=0;i<m_param->nbVar;i++)
        {
            m_leftBorder[i]=m_freedom[std::make_pair(i,0)];
            m_rightBorder[i]=m_freedom[std::make_pair(i,0)];
        }
    
    case 2:
        for(int i=0;i<m_param->nbVar;i++)
        {
            m_leftBorder[i]=m_freedom[std::make_pair(i,0)]-6*m_freedom[std::make_pair(i,1)];
            m_rightBorder[i]=m_freedom[std::make_pair(i,0)]+6*m_freedom[std::make_pair(i,1)];
        }

    case 3:
        for(int i=0;i<m_param->nbVar;i++)
        {
            m_leftBorder[i]=m_freedom[std::make_pair(i,0)]-6*m_freedom[std::make_pair(i,1)]+30*m_freedom[std::make_pair(i,2)];
            m_rightBorder[i]=m_freedom[std::make_pair(i,0)]+6*m_freedom[std::make_pair(i,1)]+30*m_freedom[std::make_pair(i,2)];
        }

    default:
        break;
    }
}