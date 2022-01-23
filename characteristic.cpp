/*#include "include/Domain.h"


void Domain::Eigenvector()
{
    double rhoL, rhoR;
    double uL, uR;
    double hL, hR;
    double u, h, a;
    double b1, b2;
    for(int cell=0;cell<m_nbCells;cell++)
    {
                   
        if(cell==m_nbCells-1)
        {
            rhoL=m_cells[cell]->getFreedom()[std::make_pair(0,0)];
            rhoR=m_cells[cell]->getFreedom()[std::make_pair(0,0)];
            
            uL=m_cells[cell]->getFreedom()[std::make_pair(1,0)];
            uR=m_cells[cell]->getFreedom()[std::make_pair(1,0)];

            hL=m_param->gamma*m_cells[cell]->getFreedom()[std::make_pair(2,0)]/rhoL-((m_param->gamma-1.0)/2.0)*uL*uL;
            hR=hL;
            
        }
        else
        {
            rhoL=m_cells[cell]->getFreedom()[std::make_pair(0,0)];
            rhoR=m_cells[cell+1]->getFreedom()[std::make_pair(0,0)];

            uL=m_cells[cell]->getFreedom()[std::make_pair(1,0)];
            uR=m_cells[cell+1]->getFreedom()[std::make_pair(1,0)];

            hL=m_param->gamma*m_cells[cell]->getFreedom()[std::make_pair(2,0)]/rhoL-((m_param->gamma-1.0)/2.0)*uL*uL;

            hR=m_param->gamma*m_cells[cell+1]->getFreedom()[std::make_pair(2,0)]/rhoR-((m_param->gamma-1.0)/2.0)*uR*uR;
        }
        
        u=(sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR));
        h=(sqrt(rhoL)*hL+sqrt(rhoR)*hR)/(sqrt(rhoL)+sqrt(rhoR));
        a=sqrt((m_param->gamma-1)*(h-u*u/2.0));
        b1=(m_param->gamma-1)*0.5*(u/a)*(u/a);
        b2=(m_param->gamma-1)/(a*a);
        //P=1
        //std::cout<<"b1= "<<b1<<std::endl;
        //std::cout<<"b2= "<<b2<<std::endl;
        //std::cout<<"uR= "<<rhoR<<std::endl;
        m_leftEigen[std::make_pair(cell,0)][0]=1.0;
        m_leftEigen[std::make_pair(cell,0)][1]=(u-a);
        m_leftEigen[std::make_pair(cell,0)][2]=(h-u*a);
        m_leftEigenVal[std::make_pair(cell,0)]=u-a;

        m_rightEigen[std::make_pair(cell,1)][0]=0.5*(b1+u/a);
        m_rightEigen[std::make_pair(cell,1)][1]=-0.5*(u*b2+1.0/a);
        m_rightEigen[std::make_pair(cell,1)][2]=b2*0.5;
        //std::cout<<"n1= "<<0.5*(b1+u/a)-0.5*(u*b2+1.0/a)*(u-a)+(h-u*a)*b2*0.5<<std::endl;
        //P=2
       
        m_leftEigen[std::make_pair(cell,1)][0]=1.0;
        m_leftEigen[std::make_pair(cell,1)][1]=u;
        m_leftEigen[std::make_pair(cell,1)][2]=(u*u/2.0);
        m_leftEigenVal[std::make_pair(cell,1)]=u;

        m_rightEigen[std::make_pair(cell,1)][0]=1.0-b1;;
        m_rightEigen[std::make_pair(cell,1)][1]=b2*u;
        m_rightEigen[std::make_pair(cell,1)][2]=-b2;

        //std::cout<<"n2= "<<1.0-b1+b2*u*(u)-(u*u/2.0)*b2<<std::endl;
        //P=3
        
        m_leftEigen[std::make_pair(cell,2)][0]=1.0;
        m_leftEigen[std::make_pair(cell,2)][1]=(u+a);
        m_leftEigen[std::make_pair(cell,2)][2]=(h+u*a);
        m_leftEigenVal[std::make_pair(cell,2)]=(u+a);

        m_rightEigen[std::make_pair(cell,1)][0]=0.5*(b1-u/a);
        m_rightEigen[std::make_pair(cell,1)][1]=-0.5*(u*b2-1.0/a);
        m_rightEigen[std::make_pair(cell,1)][2]=b2*0.5;
        //std::cout<<"n3= "<<0.5*(b1-u/a)-0.5*(u*b2-1.0/a)*(u+a)+(h+u*a)*b2*0.5<<std::endl;
    }
};

void Domain::minmodCharc()
{
    std::vector<double> Minmod(2,0.0);
    double left;
    double l1, l2, l3;
    double r1, r2,r3;
    double p1,p2,p3;
    double m1,m2,m3;
    double right;
    double deltaPlus;
    double deltaMinus;
    double f1,f2,f3, f;
    
    std::pair<int,int> key;
    for(int cell=0;cell<m_nbCells;cell++)
    {
        
        
        l1=-m_cells[cell]->getLeftBorder()[0]+m_cells[cell]->getFreedom()[std::make_pair(0,0)];
        l2=-m_cells[cell]->getLeftBorder()[1]+m_cells[cell]->getFreedom()[std::make_pair(1,0)];
        l3=-m_cells[cell]->getLeftBorder()[2]+m_cells[cell]->getFreedom()[std::make_pair(2,0)];
        
        r1=m_cells[cell]->getRightBorder()[0]-m_cells[cell]->getFreedom()[std::make_pair(0,0)];
        r2=m_cells[cell]->getRightBorder()[1]-m_cells[cell]->getFreedom()[std::make_pair(1,0)];
        r3=m_cells[cell]->getRightBorder()[2]-m_cells[cell]->getFreedom()[std::make_pair(2,0)];

        f1=m_cells[cell]->getFreedom()[std::make_pair(0,0)];
        f2=m_cells[cell]->getFreedom()[std::make_pair(1,0)];
        f3=m_cells[cell]->getFreedom()[std::make_pair(2,0)];
        
        if(cell==0)
        {
            p1=m_cells[cell+1]->getFreedom()[std::make_pair(0,0)]-m_cells[cell]->getFreedom()[std::make_pair(0,0)];
            p2=m_cells[cell+1]->getFreedom()[std::make_pair(1,0)]-m_cells[cell]->getFreedom()[std::make_pair(1,0)];
            p3=m_cells[cell+1]->getFreedom()[std::make_pair(2,0)]-m_cells[cell]->getFreedom()[std::make_pair(2,0)];
            
            m1=0.0;
            m2=0.0;
            m3=0.0;
            
        }
        else if(cell==m_nbCells-1)
        {
            p1=0.0;
            p2=0.0;
            p3=0.0;
            

            m1=m_cells[cell]->getFreedom()[std::make_pair(0,0)]-m_cells[cell-1]->getFreedom()[std::make_pair(0,0)];
            m2=m_cells[cell]->getFreedom()[std::make_pair(1,0)]-m_cells[cell-1]->getFreedom()[std::make_pair(1,0)];
            m3=m_cells[cell]->getFreedom()[std::make_pair(2,0)]-m_cells[cell-1]->getFreedom()[std::make_pair(2,0)];
            
        }
        else
        {
            p1=m_cells[cell+1]->getFreedom()[std::make_pair(0,0)]-m_cells[cell]->getFreedom()[std::make_pair(0,0)];
            p2=m_cells[cell+1]->getFreedom()[std::make_pair(1,0)]-m_cells[cell]->getFreedom()[std::make_pair(1,0)];
            p3=m_cells[cell+1]->getFreedom()[std::make_pair(2,0)]-m_cells[cell]->getFreedom()[std::make_pair(2,0)];
            

            m1=m_cells[cell]->getFreedom()[std::make_pair(0,0)]-m_cells[cell-1]->getFreedom()[std::make_pair(0,0)];
            m2=m_cells[cell]->getFreedom()[std::make_pair(1,0)]-m_cells[cell-1]->getFreedom()[std::make_pair(1,0)];
            m3=m_cells[cell]->getFreedom()[std::make_pair(2,0)]-m_cells[cell-1]->getFreedom()[std::make_pair(2,0)];
        }


        for(int j=0;j<m_param->nbVar;j++)
        {
            key=std::make_pair(cell,j);
            left=l1*m_leftEigen[key][0]+l2*m_leftEigen[key][1]+l3*m_leftEigen[key][2];
            right=r1*m_leftEigen[key][0]+r2*m_leftEigen[key][1]+r3*m_leftEigen[key][2];
            deltaPlus=p1*m_leftEigen[key][0]+p2*m_leftEigen[key][1]+p3*m_leftEigen[key][2];
            deltaMinus=m1*m_leftEigen[key][0]+m2*m_leftEigen[key][1]+m3*m_leftEigen[key][2];

            m_leftBorderCharc[key]=mimod(left,deltaPlus,deltaMinus);

            m_rightBorderCharc[key]=mimod(right,deltaPlus,deltaMinus);
            //std::cout<<"r0= "<<m_rightBorderCharc[key]<<std::endl;

            f=f1*m_leftEigen[key][0]+f2*m_leftEigen[key][1]+f3*m_leftEigen[key][2];

            m_freeCharc[key]=f;
        }
        for(int j=0;j<m_param->nbVar;j++)
        {
            key=std::make_pair(cell,j);
            m_leftBorderTrue[key]=m_leftBorderCharc[std::make_pair(cell,0)]*m_rightEigen[std::make_pair(cell,0)][j]+m_leftBorderCharc[std::make_pair(cell,1)]*m_rightEigen[std::make_pair(cell,1)][j]+m_leftBorderCharc[std::make_pair(cell,2)]*m_rightEigen[std::make_pair(cell,2)][j];

            m_rightBorderTrue[key]=m_rightBorderCharc[std::make_pair(cell,0)]*m_rightEigen[std::make_pair(cell,0)][j]+m_rightBorderCharc[std::make_pair(cell,1)]*m_rightEigen[std::make_pair(cell,1)][j]+m_rightBorderCharc[std::make_pair(cell,2)]*m_rightEigen[std::make_pair(cell,2)][j];
        }

        
    }
    
};

std::map<int,double> Domain::getLeftBorderCharc(int cell)
{

    std::map<int,double> leftBorderCharc;
    std::map<int,double> leftCharc;
    std::pair<int,int> key;
    double l1, l2 ,l3;

    double l;

    l1=m_cells[cell]->getFreedom()[std::make_pair(0,0)];
    l2=m_cells[cell]->getFreedom()[std::make_pair(1,0)];
    l3=m_cells[cell]->getFreedom()[std::make_pair(2,0)];
    for(int i=0;i<m_param->nbVar;i++)
    {
        key=std::make_pair(cell,i);
        l=l1*m_leftEigen[key][0]+l2*m_leftEigen[key][1]+l3*m_leftEigen[key][2];
        leftBorderCharc[i]=l-m_leftBorderCharc[std::make_pair(cell,i)];//minmod(cell,i)[0];
    }
    for(int i=0;i<m_param->nbVar;i++)
    {
        leftCharc[i]=leftBorderCharc[0]*m_rightEigen[std::make_pair(cell,0)][i]+leftBorderCharc[1]*m_rightEigen[std::make_pair(cell,1)][i]+leftBorderCharc[2]*m_rightEigen[std::make_pair(cell,2)][i];
    }
  
    return leftCharc;
}


std::map<int,double> Domain::getRightBorderCharc(int cell)
{
    std::map<int,double> rightBorderCharc;
    std::map<int,double> rightCharc;
    std::pair<int,int> key;
    double r1, r2 ,r3;

    double r;

    r1=m_cells[cell]->getFreedom()[std::make_pair(0,0)];
    r2=m_cells[cell]->getFreedom()[std::make_pair(1,0)];
    r3=m_cells[cell]->getFreedom()[std::make_pair(2,0)];
    for(int i=0;i<m_param->nbVar;i++)
    {
        key=std::make_pair(cell,i);
        r=r1*m_leftEigen[key][0]+r2*m_leftEigen[key][1]+r3*m_leftEigen[key][2];
        rightBorderCharc[i]=r+m_rightBorderCharc[std::make_pair(cell,i)];//minmod(cell,i)[0];
        //std::cout<<"r= "<<m_rightBorderCharc[std::make_pair(cell,i)]<<std::endl;
    }
    for(int i=0;i<m_param->nbVar;i++)
    {
        rightCharc[i]=rightBorderCharc[0]*m_rightEigen[std::make_pair(cell,0)][i]+rightBorderCharc[1]*m_rightEigen[std::make_pair(cell,1)][i]+rightBorderCharc[2]*m_rightEigen[std::make_pair(cell,2)][i];
    }
  
    return rightCharc;
}


std::map<int, double> Domain::rightFluxCharc(int i)
{
    std::map<int, double> flux;
    std::map<int, double> fluxChar;
    std::pair<int,int> key;
    double lambda1,lambda2,lambda, left, right;
    if(i==m_nbCells-1)
    {
        
        flux[0]=m_cells[i]->getFU(getRightBorderCharc(i))[0];
        
        flux[1]=m_cells[i]->getFU(getRightBorderCharc(i))[1];
        flux[2]=m_cells[i]->getFU(getRightBorderCharc(i))[2];
    }
    else
    {

        double fLeft=m_cells[i]->getFU(getRightBorderCharc(i))[0];
        double fRight=m_cells[i+1]->getFU(getLeftBorderCharc(i+1))[0];
       
        flux[0]=0.5*(fLeft+fRight);

        fLeft=m_cells[i]->getFU(getRightBorderCharc(i))[1];
        fRight=m_cells[i+1]->getFU(getLeftBorderCharc(i+1))[1];
       
        flux[1]=0.5*(fLeft+fRight);

        fLeft=m_cells[i]->getFU(getRightBorderCharc(i))[2];
        fRight=m_cells[i+1]->getFU(getLeftBorderCharc(i+1))[2];
       
        flux[2]=0.5*(fLeft+fRight);
    }
    
    for(int j=0;j<m_param->nbVar;j++)
    {
        key=std::make_pair(i,j);
        fluxChar[i]=flux[0]*m_leftEigen[key][0]+flux[1]*m_leftEigen[key][1]+flux[2]*m_leftEigen[key][2];

        if(i==m_nbCells-1)
        {
            lambda=m_rightEigenVal[key];
            right=0;
            left=0;
        }
        else
        {
            lambda1=m_rightEigenVal[key];
            lambda2=m_rightEigenVal[std::make_pair(i+1,j)];;
            lambda=std::max(lambda1,lambda2);
            //right=getLeftBorderCharc(i+1)[j];
            right=m_freeCharc[std::make_pair(i+1,j)]-m_leftBorderCharc[std::make_pair(i+1,j)];
            left=m_freeCharc[std::make_pair(i,j)]+m_leftBorderCharc[std::make_pair(i,j)];
            //left=getRightBorderCharc(i)[j];
        }
        
        
        fluxChar[i]=fluxChar[i]-0.5*lambda*(right-left);
    }
    //std::cout<<"f1= "<<flux[0]<<std::endl;
    std::map<int, double> finalFlux;
    for(int j=0;j<m_param->nbVar;j++)
    {
        if(i==m_nbCells-1)
        {
            finalFlux[j]=fluxChar[0]*m_rightEigen[std::make_pair(i,0)][j]+fluxChar[1]*m_rightEigen[std::make_pair(i,1)][j]+fluxChar[2]*m_rightEigen[std::make_pair(i,2)][j];
        }
        else
        {
            finalFlux[j]=fluxChar[0]*m_rightEigen[std::make_pair(i+1,0)][j]+fluxChar[1]*m_rightEigen[std::make_pair(i+1,1)][j]+fluxChar[2]*m_rightEigen[std::make_pair(i+1,2)][j];
        }
        //std::cout<<"f2= "<<m_rightEigen[std::make_pair(i,2)][j]<<std::endl;
    }
    //std::cout<<"f2= "<<finalFlux[0]<<std::endl;
    return finalFlux;

}




std::map<int, double> Domain::leftFluxCharc(int i)
{
    std::map<int, double> flux;
    std::map<int, double> fluxChar;
    std::pair<int,int> key;
    double lambda1,lambda2,lambda, left, right;
    if(i==0)
    {
        
        flux[0]=m_cells[i]->getFU(getLeftBorderCharc(i))[0];
        
        flux[1]=m_cells[i]->getFU(getLeftBorderCharc(i))[1];

        flux[2]=m_cells[i]->getFU(getLeftBorderCharc(i))[2];
    }
    else
    {
        
        double fLeft=m_cells[i-1]->getFU(getRightBorderCharc(i-1))[0];
        double fRight=m_cells[i]->getFU(getLeftBorderCharc(i))[0];
        
        flux[0]=0.5*(fLeft+fRight);
        

        fLeft=m_cells[i-1]->getFU(getRightBorderCharc(i-1))[1];
        fRight=m_cells[i]->getFU(getLeftBorderCharc(i))[1];
       
        flux[1]=0.5*(fLeft+fRight);

        fLeft=m_cells[i-1]->getFU(getRightBorderCharc(i-1))[2];
        fRight=m_cells[i]->getFU(getLeftBorderCharc(i))[2];
       
        flux[2]=0.5*(fLeft+fRight);
    }
    for(int j=0;j<m_param->nbVar;j++)
    {
        key=std::make_pair(i,j);
        fluxChar[i]=flux[0]*m_rightEigen[key][0]+flux[1]*m_rightEigen[key][1]+flux[2]*m_rightEigen[key][2];
        if(i==0)
        {
            lambda=m_rightEigenVal[key];
            right=0;
            left=0;
        }
        else
        {
            lambda1=m_rightEigenVal[key];
            lambda2=m_rightEigenVal[std::make_pair(i+1,j)];;
            lambda=std::max(lambda1,lambda2);
            //right=getLeftBorderCharc(i+1)[j];
            right=m_freeCharc[std::make_pair(i+1,j)]-m_leftBorderCharc[std::make_pair(i,j)];
            left=m_freeCharc[std::make_pair(i,j)]+m_leftBorderCharc[std::make_pair(i-1,j)];
            //left=getRightBorderCharc(i)[j];
        }
        
        
        fluxChar[i]=fluxChar[i]-0.5*lambda*(right-left);
       
    }
    std::map<int, double> finalFlux;
    for(int j=0;j<m_param->nbVar;j++)
    {
        if(i==m_nbCells-1)
        {
            finalFlux[j]=fluxChar[0]*m_rightEigen[std::make_pair(i,0)][j]+fluxChar[1]*m_rightEigen[std::make_pair(i,1)][j]+fluxChar[2]*m_rightEigen[std::make_pair(i,2)][j];
        }
        else
        {
            finalFlux[j]=fluxChar[0]*m_leftEigen[std::make_pair(i+1,0)][j]+fluxChar[1]*m_leftEigen[std::make_pair(i+1,1)][j]+fluxChar[2]*m_leftEigen[std::make_pair(i+1,2)][j];
        }
    }
    
    return finalFlux;
}

void Domain::updateCharc()
{
    double a;
    for(int cell=0;cell<m_nbCells;cell++)
    {
        for(int i=0;i<m_param->nbVar;i++)
        {
            if(m_param->Order==2)
            {
                a=(1.0/6.0)*m_leftBorderTrue[std::make_pair(cell,i)];
                m_cells[cell]->updateMod(i,1,a);
            }
            if(m_param->Order==3)
            {
                a=(1.0/12.0)*(m_leftBorderTrue[std::make_pair(cell,i)]+m_rightBorderTrue[std::make_pair(cell,i)]);
                m_cells[cell]->updateMod(i,1,a);
                a=(1.0/60.0)*(m_rightBorderTrue[std::make_pair(cell,i)]-m_leftBorderTrue[std::make_pair(cell,i)]);
                m_cells[cell]->updateMod(i,2,a);
            }
        }
    }
}
*/