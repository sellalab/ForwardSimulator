//This code was used to run the simulations from Amorim et al. (2017) - The Population Genetics of Human Disease: The Case of Recessive Lethal Mutations
//Modified by Eduardo Amorim (guerraamorim@gmail.com) from the original of Yuval Simons (see reference 1, Simons et al. 2014).

#include <chrono>
#include <random>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "population.h"
#include "BRand.hpp"
#include <stdint.h>
#include <map>
#include <time.h>
#include <sstream>
#include <fstream>
#include <string>
#define tau 2040 
#define taujump 5920 
#define Urate 0.000000012 
#define ratio 1  

double lognormal();
double Uvar;

int nfinal,RUNS,Ntau=14474; 
double sel,DOM;

using namespace std;

struct freqs 
{ 
  double freq0; 
  double freq1; 

  freqs(const double a=0,const double b=0) :  
    freq0(a),freq1(b) {}   
};

class valueComp 
{
public:
  bool operator()(const freqs& A,
          const freqs& B)
    const
    { if (A.freq0!=B.freq0) 
         return A.freq0<B.freq0;
      else
         return A.freq1<B.freq1; }
};

char* filename(char *);
char* seriesfilename(char * str);
int popsize(int gen, int i);
int poi(double l);
int poicond1(double l);
int binom(int n, double p);
int binom1(int n, double p);
double lognormal();

int main(int argc, char *argv[])
{
    int initgen=150000,stopover;  
    double m=0.00015;  
    if (argc==4) // User may run the script with "./a.out RUNS sel DOM", where RUNS = number of runs, sel = selective coefficient (in this work always set to 1), and DOM = heterozygote effect (in this work set to 0 or 1%)  
       {RUNS=atoi(argv[1]);
       sel=atof(argv[2]);
       DOM=atof(argv[3]); 
       }
    else         // User may run the script with "./a.out" and the above mentioned parameters will be asked in the command line prompt.
        {
        cout<<"enter number of runs:";
        cin>>RUNS;
        cout<<"enter selective coefficient:";
        cin>>sel;
        cout<<"enter dominance:";
        cin>>DOM;
        }
    stopover=RUNS;
    cout<<"RUNS:"<<RUNS<<" tau:"<<tau<<" sel:"<<sel<<" initgen:"    
    <<initgen<<" stopover:"<<stopover<<" dominace:"<<DOM<<" m:"<<m<<"\n";
 
    BRand::Controller.seed(time(NULL));
    
    population::initialize(sel,DOM); 
    population* pops[2]; 
    pops[0]= new population(popsize(0,0));
    pops[1]= new population(popsize(0,0));

    int euroflag=0;  
    double halftheta=2*popsize(initgen,0)*Urate; 
    double p1=1./(1.+ratio*exp((2*popsize(initgen,0)-1)*sel)*(1-(2*DOM-1)*(2*popsize(initgen,0)-1)*sel*sel/(12*popsize(initgen,0))) ); 
    map<int,int> count[2],counts[2][2],countd[2],taucount,taucountd;
    map<freqs,int,valueComp> joint,joints[2],jointd;

    int gen,thispop;
    double baseup=(1./halftheta); 
    double basedown=(1./(ratio*halftheta));
    double oneovertwoU=(1./(2*Urate));
    double oneovertwoUratio=(1./(2*Urate*ratio));
    int initialstate,zeroruns=0,oneruns=0;
    int skip=0,stretch;
    time_t tt;
    struct tm *tim;
    char library[70],comm[200],totalfile[200],totalres[200];
    double Ep,s1,S2N,Ep2,Ex,Exderived,s12,Ex2,Ex4,Dx,Dx2;
    double europ, euros1, eurox, euroS2N;
    int tauallele,before, after;
    freqs f;
    
    ofstream myfile,results,afrDistrib,eurDistrib,mutDistribution;
    
    for(int run=0;run<RUNS;run++)  
            {					  
            gen=initgen;
            double Uvar = lognormal();
            double halfthetaVAR = 2*popsize(initgen,0)*Uvar;
            double baseupVAR = (1./halfthetaVAR);
            
            pops[0]->size=popsize(gen,0); //Initializing the population size for the beginning of each run with deleterious allele absent from the population

                pops[0]->clear();
                initialstate=0;
                zeroruns++;
            
            euroflag=0; //A flag indicating if the African-European population split has occured (0 = it hasn't occurred)

            while (gen>0)
                    {
		      if (gen==tau) //tau is the generation of the out of Africa exodus (gen=2040) according to Tennessen et al (2012).
                    {
                        euroflag=1;
                        *pops[1]=*pops[0]; 
                        taucount[pops[0]->allelenum()]++; 
                        if (initialstate==0)
			  			taucountd[pops[0]->allelenum()]++; 
                        else 
                    	taucountd[2*pops[0]->size-pops[0]->allelenum()]++;
                    }
		      if (gen==920) m=0.000025; 
                    
		      for(int i=0;i<(1+euroflag);i++) 
                    { 
		      if (((pops[i]->allelenum()>0)&&(pops[i]->allelenum()<(2*pops[i]->size)))||(gen<=taujump)) //If the allele is segregating or we are in recent history (last 5920 generations), introduce mutations
                        {    
                            pops[i]->mutateup(poi( Uvar*(2*pops[i]->size-pops[i]->allelenum()) )) ; 
                        }
		      else if (pops[i]->allelenum()==0) //If the ancestral allele is fixed calculate when will the next deleterious allele appear and advance population state to that point in time (unless that point in time is in recent history; then just set time to initial history).
                        {
                             gen-=int(-log( BRand::Controller.nextOpened())*baseupVAR); 
                             if (gen<=taujump) 
                                gen=taujump+1; 
                             else
                                 {pops[i]->size=popsize(gen,i);
                                 pops[i]->clear(); 
                                 pops[i]->mutateup(poicond1(2*Uvar*pops[i]->size)); } 

                        }    
                        else if (pops[i]->allelenum()==(2*pops[i]->size))//If the deleterious allele is fixed calculate when will the next deleterious allele appear and advance population state to that point in time (unless that point in time is in recent history and then just set time to initial history). Note that with s = 1 this cannot happen
                        {
                             gen-=int(-log( BRand::Controller.nextOpened())*basedown);
                             if (gen<=taujump) 
                                gen=taujump+1;
                             else
                                 {pops[i]->size=popsize(gen,i);
                                 pops[i]->fix();
                                 }
                        }       
                                       
		      if (euroflag==0) //If a split between Af and Eur populations has occured apply migration
                           pops[i]->populate_from(pops[i]->prob(),popsize(gen-1,i));
                        else
                           pops[i]->populate_from((1-m)*(pops[i]->prob())+m*(pops[1-i]->prob()),popsize(gen-1,i)); 
                    }
                    
                                    
		      gen--; //Next generation please
}
	    //End of a run. Time for some statitics.
	    
          for(int i=0;i<2;i++) 
	      {count[i][pops[i]->allelenum()]++; 
                if (initialstate==0)
                   {counts[i][initialstate][pops[i]->allelenum()]++; 
                   countd[i][pops[i]->allelenum()]++;} 
                else 
		  {counts[i][initialstate][2*pops[i]->size-pops[i]->allelenum()]++;
                   countd[i][2*pops[i]->size-pops[i]->allelenum()]++;}
               }
            f.freq0=double(pops[0]->allelenum())/(2*pops[0]->size); //Allele frequency in Africans
            f.freq1=double(pops[1]->allelenum())/(2*pops[1]->size); //Allele frequency in Europeans
            joint[f]++; 
            if (initialstate==0) 
                 {joints[0][f]++;
                 jointd[f]++;}
            else 
                 {f.freq0=1-f.freq0;
                 f.freq1=1-f.freq1;
                 joints[1][f]++;
                 jointd[f]++;}
            
           	afrDistrib.open("afrDistrib.txt", ios::app);
            afrDistrib<<f.freq0;
            afrDistrib<<"\n";
            afrDistrib.close();
            
			eurDistrib.open("eurDistrib.txt", ios::app);
            eurDistrib<<f.freq1;
            eurDistrib<<"\n";
            eurDistrib.close();
            
			mutDistribution.open("mutDistribution.txt", ios::app);
            mutDistribution<<Uvar;
            mutDistribution<<"\n";
            mutDistribution.close();
            
            if (((100*(run+1))%RUNS)==0) cout<<((100.0*(run+1))/RUNS)<<"%\n";              

            if (((run+1)%stopover==0)||((run+1)==RUNS))  
            {

	      time(&tt); //We use a time stamp to differentiate files in parallel computing
                tim=localtime(&tt);
                sprintf(library, "Tennessen/regular");
                sprintf(comm, "mkdir %s",library);
                system(comm);
                sprintf(library, "%s/res%d_%d_%d",library, tim->tm_mday,tim->tm_mon+1,tim->tm_year+1900);
                sprintf(comm, "mkdir %s",library);
                system(comm);
                sprintf(library, "%s/dom%f",library, DOM);
                sprintf(comm, "mkdir %s",library);
                system(comm);

		//Printing statistics and data
                sprintf(totalfile, "%s/Tennessendata_sel%f_dom%f_runs%d_%d_%d_%d.csv", library,log10(sel),DOM,run+1,tim->tm_hour,tim->tm_min,tim->tm_sec);
                sprintf(totalres, "%s/Tennessen_sel%f_dom%f_runs%d_%d_%d_%d.csv", library,log10(sel),DOM,run+1,tim->tm_hour,tim->tm_min,tim->tm_sec);
                cout<<"opening "<<totalfile<<"\n";
                myfile.open (totalfile); //Data
                cout<<"opened "<<totalfile<<"\n";
                cout<<"opening "<<totalres<<"\n";
                results.open (totalres); //Statistics
                cout<<"opened "<<totalres<<"\n";
                results<<"Ntau="<<Ntau<<",tau="<<tau<<",halftheta="<<halftheta<<",p1="<<p1<<",Urate="<<Urate<<",ratio="<<ratio<<",initgen="<<initgen<<",sel="<<sel<<",DOM="<<DOM<<",runs="<<(run+1)<<"\n";
                myfile<<"Ntau="<<Ntau<<",tau="<<tau<<",halftheta="<<halftheta<<",p1="<<p1<<",Urate="<<Urate<<",ratio="<<ratio<<",initgen="<<initgen<<",sel="<<sel<<",DOM="<<DOM<<",runs="<<(run+1)<<"\n";
                for(int i=0;i<2;i++)                
                    {Ep=0; Ep2=0;
                    s1=0; s12=0;
                    Ex=0;Ex2=0;Ex4=0;Dx=0;Dx2=0;
                    for(map<int,int>::iterator it=count[i].begin();it!=count[i].end();it++) {Ep+=(double(it->first)/(2*pops[i]->size))*double(it->second)/(run+1);
                                                                                            Ep2+=pow(double(it->first)/(2*pops[i]->size),2)*double(it->second)/(run+1);
                                                                                            }
                    S2N=1-double(count[i][0])/(run+1)-double(count[i][2*popsize(0,i)])/(run+1);
                    for(map<int,int>::iterator it=countd[i].begin();it!=countd[i].end();it++) {s1+=(double(it->first)/(2*pops[i]->size))*double(it->second)/(run+1);
                                                                                               s12+=pow(double(it->first)/(2*pops[i]->size),2)*double(it->second)/(run+1);
                                                                                               if ((it->first!=0)&&(it->first!=pops[i]->size))
                                                                                                  {Ex+=(double(it->first)/(2*pops[i]->size))*double(it->second)/(run+1);
                                                                                                  Ex2+=pow(double(it->first)/(2*pops[i]->size),2)*double(it->second)/(run+1);
                                                                                                  Ex4+=pow(double(it->first)/(2*pops[i]->size),4)*double(it->second)/(run+1);
                                                                                                  }
                                                                                               }

    
                    Ex=Ex/S2N;
                    Ex2=Ex2/S2N;
                    Ex4=Ex4/S2N;                
                    Dx=sqrt((Ex2-pow(Ex,2))/(S2N*(run+1)));
                    Dx2=sqrt((Ex4-pow(Ex2,2))/(S2N*(run+1)));                
                    
                    results<<"Ex,Dx,Ex2,Dx2,s1,Ds1,Ep,Dp,S2N\n";//myfile << "Writing this to a file.\n";    
                    results<<Ex<<","<<Dx<<","<<Ex2<<","<<Dx2<<","<<s1<<","<<sqrt((s12-pow(s1,2))/(run+1))<<","<<Ep<<","<<sqrt((Ep2-pow(Ep,2))/(run+1))<<","<<S2N<<"\n";//myfile << "Writing this to a file.\n";
                    myfile.precision(10);
                    myfile<<"deleterious allele frequencies\n";
                    for(map<int,int>::iterator it=count[i].begin();it!=count[i].end();it++) myfile<<(double(it->first)/(2*popsize(0,i)))<<",";
                    myfile<<"\n";
                    for(map<int,int>::iterator it=count[i].begin();it!=count[i].end();it++) myfile<<double(it->second)/(run+1)<<",";
                    myfile<<"\n";
                    myfile<<"derived allele frequencies conditional on zero\n";
                    for(map<int,int>::iterator it=counts[i][0].begin();it!=counts[i][0].end();it++) myfile<<(double(it->first)/(2*popsize(0,i)))<<",";
                    myfile<<"\n";
                    for(map<int,int>::iterator it=counts[i][0].begin();it!=counts[i][0].end();it++) myfile<<double(it->second)/zeroruns<<",";
                    myfile<<"\n";
                    if (oneruns!=0)
                        {myfile<<"derived allele frequencies conditional on one\n";
                        for(map<int,int>::iterator it=counts[i][1].begin();it!=counts[i][1].end();it++) myfile<<(double(it->first)/(2*popsize(0,i)))<<",";
                        myfile<<"\n";
                        for(map<int,int>::iterator it=counts[i][1].begin();it!=counts[i][1].end();it++) myfile<<double(it->second)/oneruns<<",";
                        myfile<<"\n";}
                    myfile<<"derived allele frequencies\n";
                    for(int j=0;j<(2*popsize(0,i)+1);j++) if ((counts[i][0][j]+counts[i][1][j])!=0) myfile<<double(j)/(2*popsize(0,i))<<",";    
                    myfile<<"\n";
                    for(int j=0;j<(2*popsize(0,i)+1);j++) if ((counts[i][0][j]+counts[i][1][j])!=0) myfile<<double(counts[i][0][j]+counts[i][1][j])/(run+1)<<",";    
                    myfile<<"\n";
                    }
                    
                    
                    myfile<<"joint deleterious allele frequencies\n";
                    for(map<freqs,int,valueComp>::iterator it=joint.begin();it!=joint.end();it++) myfile<<(it->first).freq0<<",";
                    myfile<<"\n";
                    for(map<freqs,int,valueComp>::iterator it=joint.begin();it!=joint.end();it++) myfile<<(it->first).freq1<<",";
                    myfile<<"\n";
                    for(map<freqs,int,valueComp>::iterator it=joint.begin();it!=joint.end();it++) myfile<<double(it->second)/(run+1)<<",";
                    myfile<<"\n";
                    myfile<<"joint derived allele frequencies conditional on zero\n";
                    for(map<freqs,int,valueComp>::iterator it=joints[0].begin();it!=joints[0].end();it++) myfile<<(it->first).freq0<<",";
                    myfile<<"\n";
                    for(map<freqs,int,valueComp>::iterator it=joints[0].begin();it!=joints[0].end();it++) myfile<<(it->first).freq1<<",";
                    myfile<<"\n";
                    for(map<freqs,int,valueComp>::iterator it=joints[0].begin();it!=joints[0].end();it++) myfile<<double(it->second)/zeroruns<<",";
                    myfile<<"\n";
                                if (oneruns!=0)
                        {myfile<<"joint derived allele frequencies conditional on one\n";
                        for(map<freqs,int,valueComp>::iterator it=joints[1].begin();it!=joints[1].end();it++) myfile<<(it->first).freq0<<",";
                        myfile<<"\n";
                        for(map<freqs,int,valueComp>::iterator it=joints[1].begin();it!=joints[1].end();it++) myfile<<(it->first).freq1<<",";
                        myfile<<"\n";
                        for(map<freqs,int,valueComp>::iterator it=joints[1].begin();it!=joints[1].end();it++) myfile<<double(it->second)/oneruns<<",";
                        myfile<<"\n";}
                    myfile<<"joint derived allele frequencies\n";
                    for(map<freqs,int,valueComp>::iterator it=jointd.begin();it!=jointd.end();it++) myfile<<(it->first).freq0<<",";
                    myfile<<"\n";
                    for(map<freqs,int,valueComp>::iterator it=jointd.begin();it!=jointd.end();it++) myfile<<(it->first).freq1<<",";
                    myfile<<"\n";
                    for(map<freqs,int,valueComp>::iterator it=jointd.begin();it!=jointd.end();it++) myfile<<double(it->second)/(run+1)<<",";
                    myfile<<"\n";
                    
                    myfile<<"split deleterious allele frequencies\n";
                    for(map<int,int>::iterator it=taucount.begin();it!=taucount.end();it++) myfile<<(double(it->first)/(2*Ntau))<<",";
                    myfile<<"\n";
                    for(map<int,int>::iterator it=taucount.begin();it!=taucount.end();it++) myfile<<double(it->second)/(run+1)<<",";
                    myfile<<"\n";
                    myfile<<"split derived allele frequencies\n";
                    for(map<int,int>::iterator it=taucountd.begin();it!=taucountd.end();it++) myfile<<(double(it->first)/(2*Ntau))<<",";
                    myfile<<"\n";
                    for(map<int,int>::iterator it=taucountd.begin();it!=taucountd.end();it++) myfile<<double(it->second)/(run+1)<<",";
                    myfile<<"\n";
                //myfile<<"check deletrious allele frequencies\n";
                //for(int i=0;i<(2*popsize(0)+1);i++) if ((counts[0][i]+counts[1][2*popsize(0)-i])!=0) myfile<<double(i)/(2*popsize(0))<<",";    
                //myfile<<"\n";
                //for(int i=0;i<(2*popsize(0)+1);i++) if ((counts[0][i]+counts[1][2*popsize(0)-i])!=0) myfile<<double(counts[0][i]+counts[1][2*popsize(0)-i])/RUNS<<",";    
                //myfile<<"\n";        
                myfile.close();
                cout<<"closed myfile\n";
                results.close();
                cout<<"closed results\n";
		

		//Printing results sub-sampled to be comparable with experimental data. Also includes admixture to simulate the difference between AF and Af-Am
                sprintf(totalfile, "%s/TennessenBinomed_sel%f_dom%f_runs%d_%d_%d_%d.csv", library,log10(sel),DOM,run+1,tim->tm_hour,tim->tm_min,tim->tm_sec);
                cout<<"opening "<<totalfile<<"\n";
                myfile.open (totalfile);
                cout<<"opened "<<totalfile<<"\n";
                myfile<<"Ntau="<<Ntau<<",tau="<<tau<<",halftheta="<<halftheta<<",p1="<<p1<<",Urate="<<Urate<<",ratio="<<ratio<<",initgen="<<initgen<<",sel="<<sel<<",DOM="<<DOM<<",runs="<<(run+1)<<"\n";
                map<freqs,int,valueComp> binomed;
                for(map<freqs,int,valueComp>::iterator it=jointd.begin();it!=jointd.end();it++)
                                                       for(int j=0; j<(it->second);j++)
                                                               {f.freq0=double(binom(435,(it->first).freq1)+binom(1741,(it->first).freq0))/2176; 
                                                               f.freq1=double(binom(2176,(it->first).freq1))/2176;
                                                               binomed[f]++;}
                s1=0; s12=0;S2N=0;
                Ex=0;Ex2=0;Ex4=0;Dx=0;Dx2=0;
                for(map<freqs,int,valueComp>::iterator it=binomed.begin();it!=binomed.end();it++) {s1+=((it->first).freq0)*(it->second);
                                                                                                 
                                                                                                  s12+=pow((it->first).freq0,2)*(it->second);
                                                                                                  if ((((it->first).freq0)!=0)&&(((it->first).freq0)!=1)) {Ex+=((it->first).freq0)*(it->second);
                                                                                                                                                          S2N+=(it->second);
                                                                                                                                                          Ex2+=pow((it->first).freq0,2)*(it->second);
                                                                                                                                                          Ex4+=pow((it->first).freq0,4)*(it->second);
                                                                                                                                                          }
                                                                                                  }
                s1=s1/(run+1);s12=s12/(run+1); S2N=S2N/(run+1);
                Ex=Ex/(run+1);Ex2=Ex2/(run+1);Ex4=Ex4/(run+1);
                Ex=Ex/S2N;Ex2=Ex2/S2N;Ex4=Ex4/S2N;                
                Dx=sqrt((Ex2-pow(Ex,2))/(S2N*(run+1)));
                Dx2=sqrt((Ex4-pow(Ex2,2))/(S2N*(run+1)));                 
                myfile<<"African Binomed Results\n";
                myfile<<"Ex, Dx, Ex2, Dx2, s1, Ds1,S2N\n";
                myfile<<Ex<<","<< Dx<<","<< Ex2<<","<< Dx2<<","<< s1<<","<< sqrt((s12-pow(s1,2))/(run+1))<<","<< S2N<<"\n";
                
                s1=0; s12=0;S2N=0;
                Ex=0;Ex2=0;Ex4=0;Dx=0;Dx2=0;
                for(map<freqs,int,valueComp>::iterator it=binomed.begin();it!=binomed.end();it++) {s1+=((it->first).freq1)*(it->second);
                                                                                                  //myfile<<"European,"<<(it->first).freq0<<","<<(it->first).freq1<<","<<it->second<<"\n";
                                                                                                  s12+=pow((it->first).freq1,2)*(it->second);
                                                                                                  if ((((it->first).freq1)!=0)&&(((it->first).freq1)!=1)) {Ex+=((it->first).freq1)*(it->second);
                                                                                                                                                          S2N+=(it->second);
                                                                                                                                                          Ex2+=pow((it->first).freq1,2)*(it->second);
                                                                                                                                                          Ex4+=pow((it->first).freq1,4)*(it->second);
                                                                                                                                                          }
                                                                                                  }
                s1=s1/(run+1);s12=s12/(run+1); S2N=S2N/(run+1);
                Ex=Ex/(run+1);Ex2=Ex2/(run+1);Ex4=Ex4/(run+1);
                Ex=Ex/S2N;Ex2=Ex2/S2N;Ex4=Ex4/S2N;  
                Dx=sqrt((Ex2-pow(Ex,2))/(S2N*(run+1)));
                Dx2=sqrt((Ex4-pow(Ex2,2))/(S2N*(run+1)));                 
                myfile<<"European Binomed Results\n";
                myfile<<"Ex, Dx, Ex2, Dx2, s1, DS1,S2N\n";
                myfile<<Ex<<","<< Dx<<","<< Ex2<<","<< Dx2<<","<< s1<<","<< sqrt((s12-pow(s1,2))/(run+1))<<","<< S2N<<"\n";
                
                S2N=0;
                Ex=0;Ex2=0;Ex4=0;Dx=0;Dx2=0;
                double Eurx=0, Eurx2=0, Eurx4=0, Durx=0, Durx2=0;
                for(map<freqs,int,valueComp>::iterator it=binomed.begin();it!=binomed.end();it++) if ((((it->first).freq0+(it->first).freq1)!=0)&&(((it->first).freq0+(it->first).freq1)!=2)) {Ex+=((it->first).freq0)*(it->second);

                                                                                                                                                                                              S2N+=(it->second);
                                                                                                                                                                                              Ex2+=pow((it->first).freq0,2)*(it->second);
                                                                                                                                                                                              Ex4+=pow((it->first).freq0,4)*(it->second);
                                                                                                                                                                                              Eurx+=((it->first).freq1)*(it->second);
                                                                                                                                                                                              Eurx2+=pow((it->first).freq1,2)*(it->second);
                                                                                                                                                                                              Eurx4+=pow((it->first).freq1,4)*(it->second);
                                                                                                                                                                                              }
                                                                                                  
                S2N=S2N/(run+1);
                Ex=Ex/(run+1);Ex2=Ex2/(run+1);Ex4=Ex4/(run+1);
                Eurx=Eurx/(run+1);Eurx2=Eurx2/(run+1);Eurx4=Eurx4/(run+1);
                Ex=Ex/S2N;Ex2=Ex2/S2N;Ex4=Ex4/S2N;
                Eurx=Eurx/S2N;Eurx2=Eurx2/S2N;Eurx4=Eurx4/S2N;
                Dx=sqrt((Ex2-pow(Ex,2))/(S2N*(run+1)));
                Dx2=sqrt((Ex4-pow(Ex2,2))/(S2N*(run+1)));                 
                Durx=sqrt((Eurx2-pow(Eurx,2))/(S2N*(run+1)));
                Durx2=sqrt((Eurx4-pow(Eurx2,2))/(S2N*(run+1)));                 
                myfile<<"joined Binomed Results\n";
                myfile<<"Ex, Dx, Ex2, Dx2,S2N\n";
                myfile<<Ex<<","<< Dx<<","<< Ex2<<","<< Dx2<<","<< S2N<<"\n";
                myfile<<"Eurx, Durx, Eurx2, Durx2,S2N\n";
                myfile<<Eurx<<","<< Durx<<","<< Eurx2<<","<< Durx2<<","<< S2N<<"\n";
                
                
                myfile<<"joint binomed derived allele frequencies\n";
                for(map<freqs,int,valueComp>::iterator it=binomed.begin();it!=binomed.end();it++) myfile<<(it->first).freq0<<",";
                myfile<<"\n";
                for(map<freqs,int,valueComp>::iterator it=binomed.begin();it!=binomed.end();it++) myfile<<(it->first).freq1<<",";
                myfile<<"\n";
                for(map<freqs,int,valueComp>::iterator it=binomed.begin();it!=binomed.end();it++) myfile<<double(it->second)/(run+1)<<",";
                myfile<<"\n";
                 
                myfile.close();
                cout<<"closed file\n";               
                
            }
    }

    delete pops[0];
    delete pops[1];
    cout<<"deleted pops\n";

    return EXIT_SUCCESS;
}


int popsize(int gen, int i)  // This version of the code generates an executable for running the simulations according to the Tennessen et al (2012) demography.
{							 // In our paper, there were cases where we used different population sizes and demography. This is the place where changes were made to accomodate the differences between scenarios.
    static double lam1=log(9300.0/1032)/(920-205);
    static double lam2=log(512000.0/9300)/205;  
    //static double lam2=log(1024000.0/9300)/205;
    //static double lam2=log(2048000.0/9300)/205;  
    //static double lam2=log(5120000.0/9300)/205;  
    static double lambda=log(424000.0/Ntau)/205;
   
    if (gen>tau)
       if (gen>taujump)
          return 7310;
       else
          return Ntau;
    else if (i==0)
         if (gen>205) 
            return Ntau;
         else
            return (int)round( Ntau*exp(lambda *(205-gen))); 
    else 
        if (gen>920) 
           return 1861;
        else if (gen>205)
           return (int)round( 1032*exp(lam1 *(920-gen)));
        else 
           return (int)round( 9300*exp(lam2 *(205-gen))); 

}

double lognormal()
{
  
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  
  std::default_random_engine generator (seed);

  // The following lines have to be uncommented or commented prior to compilation according to the type of mutation that is going to be simulated.
  // There are 4 types (CpGti, CpGtv, nonCpGti and nonCpGtv) and the genomic average, based on Kong et al. (2012) - See Methods for detail
  // We also apply the correction from Harpak et al. (2016) - See Methods for detail.
  
  //std::normal_distribution<double> distribution(-7.324836926,0.57); //CpGti 
  //std::normal_distribution<double> distribution(-8.392236341,0.57); //CpGtv
  //std::normal_distribution<double> distribution(-8.583066473,0.57); //nonCpGti 
  //std::normal_distribution<double> distribution(-8.798867103,0.57); //nonCpGtv 
  std::normal_distribution<double> distribution(-7.920818754,0.57); //average u 
  
  double u = distribution(generator);

  return pow(10,u);
}

int poi(double l) //Regular Poisson random variate
{double p=1,expl=exp(-l);
int k=0;
while (p>=expl)
{k++;
p=p*BRand::Controller.nextClosed();
}
return k-1;
}

int poicond1(double l) //Poisson random variate conditional on the result being at least 1
{double expl=exp(-l);
double p=expl+BRand::Controller.nextClosed()*(1-expl);
int k=1;
while (p>=expl)
{k++;
p=p*BRand::Controller.nextClosed();
}
return k-1;
}

int binom(int n, double p) //binomial random variate generator
{ if (p<0.5) return binom1(n,p);
else return  n-binom1(n,1-p);
}

int binom1(int n, double p) //binomial random variate generator
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
}
