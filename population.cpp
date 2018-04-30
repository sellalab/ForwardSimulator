// Class automatically generated by Dev-C++ New Class wizard
//#include "mtrand.h"
#include "population.h" // class's header file
#include "BRand.hpp"
#include <math.h>
#include <vector>
#include <iostream>
#include <cstdlib>

#define PI 3.141592654

double population::s,population::h,population::hs;

void population::initialize(double sel, double dom)
{
s=sel;
h=dom;
hs=h*s;                           
}

// class constructor
population::population(int N)
{
    alleleholders[0]=N;
    alleleholders[1]=0;
    alleleholders[2]=0;
    size=N;
}

void population::mutateup(int n)
{
     for (int i=0;i<n;i++)
     {
         if ((alleleholders[0]!=0)&&(BRand::Controller.nextOpened()>((0.5*alleleholders[1])/alleleholders[0])))
         {alleleholders[0]-=1;
         alleleholders[1]+=1;}
         else if (alleleholders[1]!=0)
         {alleleholders[1]-=1;
         alleleholders[2]+=1;}
     }         

}

void population::mutatedown(int n)
{
     for (int i=0;i<n;i++)
     {
         if ((alleleholders[2]!=0)&&(BRand::Controller.nextOpened()>((0.5*alleleholders[1])/alleleholders[2])))
         {alleleholders[2]-=1;
         alleleholders[1]+=1;}
         else if (alleleholders[1]!=0)
         {alleleholders[1]-=1;
         alleleholders[0]+=1;}
     }          

}

void population::clear()
{
    alleleholders[0]=size;
    alleleholders[1]=0;
    alleleholders[2]=0;
}

void population::fix()
{
    alleleholders[0]=0;
    alleleholders[1]=0;
    alleleholders[2]=size;
}

double population::prob()
{
       return (alleleholders[0]+(1.-hs)*alleleholders[1]*0.5)/(alleleholders[0]+(1.-hs)*alleleholders[1]+(1.-s)*alleleholders[2]);
}

void population::populate_from(population &p, int N)
{
    populate_from(p.prob(),N);
}

void population::populate_from(double prob, int N)
{
    alleleholders[0]=0;
    alleleholders[1]=0;
    alleleholders[2]=0;    
    if (N!=0) size=N;
    if (prob<0.5) //drawing a binomial variate works faster when p<0.5
       {alleleholders[0]=binom(size,prob*prob);
       alleleholders[1]=binom(size-alleleholders[0],2*prob/(1+prob));
       alleleholders[2]=size-alleleholders[0]-alleleholders[1];
       }
    else  //drawing a binomial variate works faster when p<0.5
       {prob=1-prob;
       alleleholders[2]=binom(size,prob*prob);
       alleleholders[1]=binom(size-alleleholders[2],2*prob/(1+prob));
       alleleholders[0]=size-alleleholders[2]-alleleholders[1];
       }        
}

int population::allelenum()
{
return (alleleholders[1]+2*alleleholders[2]);
}

/*int population::binom(int n, double p)
{int x=0;
double y=0,c=log(1-p);
if (c==0)
   {//std::cout<<"c is zero!\n";
   return 0;}
while (1) 
    {
    y+=ceil(log(BRand::Controller.nextOpened())/c);
    if (y>n) return x;    
    x++;    
    } 
}*/
int population::binom(int n, double pp)
{
int j;
static int nold=(-1);
double am,em,g,angle,p,bnl,sq,t,y;
static double pold=(-1.0),pc,plog,pclog,en,oldg;
p=(pp <= 0.5 ? pp : 1.0-pp);
am=n*p;
if (n < 25) 
  {bnl=0.0;
  for (j=1;j<=n;j++)
  if (BRand::Controller.nextOpened() < p) ++bnl; } 
 else if (am < 1.0) //Note to self: tried to change this 1.0 limit. Simulation perofrmance not senstive to small changes, eg 3.0. 
  {g=exp(-am);
   t=1.0;
   for (j=0;j<=n;j++) 
      {t *= BRand::Controller.nextOpened();
      if (t < g) break;}
   bnl=(j <= n ? j : n); } 
else 
  {if (n != nold) 
     { en=n;
       oldg=lgamma(en+1.0); //was gammaln
     nold=n;} 
  if (p != pold) 
     {pc=1.0-p;
     plog=log(p);
     pclog=log(pc);
     pold=p;} 
  sq=sqrt(2.0*am*pc); 
  do 
     {do 
        {angle=PI*BRand::Controller.nextOpened();
        y=tan(angle);
        em=sq*y+am;
        } 
     while (em < 0.0 || em >= (en+1.0));
     em=floor(em); 
     t=1.2*sq*(1.0+y*y)*exp(oldg-lgamma(em+1.0)-lgamma(en-em+1.0)+em*plog+(en-em)*pclog);//was gammaln
     // exp(oldg-gammln(em+1.0)-gammln(en-em+1.0)+em*plog+(en-em)*pclog) is p(em|n,p) and acts as the target PDF for em which a continuous variable. Therefore floor(em) is distributed binomially.
     // sq*(1.0+y*y) is the actual PDF from which we're drawing em and we're using acception/rejection to get the target PDF.
   } 
   while ( BRand::Controller.nextOpened() > t); 
   bnl=em;
   }
if (p != pp) bnl=n-bnl; 
return (int) bnl;
}

 
population& population::operator= (const population &Source)
{
    // do the copy
    size = Source.size;
    alleleholders[0]=Source.alleleholders[0];
    alleleholders[1]=Source.alleleholders[1];
    alleleholders[2]=Source.alleleholders[2];    
    // return the existing object
    return *this;
}

int population::choose_allele()
{double d=BRand::Controller.nextClosed();
if (d<(double(alleleholders[0])/size))
   return 0;
else if (d>(1-double(alleleholders[2])/size))
   return 2;
else
   return 1;

}

population::~population()
{
	std::cout<<"destructor called\n";

}
