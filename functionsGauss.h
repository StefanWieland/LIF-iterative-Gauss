#pragma once
#include "constants.h"

MTRand mtrand1;
const double Pi=M_PI;
std::fstream data;

double rea[tWin]={};
fftw_complex com[tWin/2+1]={};
double Sxx[tWin/2+1]={};
double SxxTemp[tWin/2+1]={};

double noiseterm[tMax]={};//noise vector
fftw_complex eta[tMax/2+1]={};//input vector for FT backtransform to noise vector

class neuronClass
{
    public:
	double state[2];
	long int lastSpike;
        neuronClass()
	{
		state[0]=vR+(vTh-vR)*mtrand1();//track membrane potential
		state[1]=0;//can be any additional variable
		lastSpike=-int(rp/dt);
	};
        ~neuronClass(){};
};

class spikeT
{
    public:
	std::vector<long int> data;
        spikeT(){};
        ~spikeT(){};
};

void getSxx(double rateGlobal, double *SxxTemp, spikeT *spiketrains)//computing power spectra, consider dimensionless spiketrains
{	
	fftw_plan forward = fftw_plan_dft_r2c_1d(tWin, rea, com, FFTW_ESTIMATE);
	for(long int run=0;run<runs;run++)
	{
	        for(int v=0;v<tWin;v++)//initialize
	   	    	rea[v]=-rateGlobal*dt;//actually /dt
		for(long int hi=0; hi<spiketrains[run].data.size(); hi++)//record
			rea[spiketrains[run].data.at(hi)-tStart]+=1;//actually /dt
		fftw_execute(forward);
		for(long int v=0;v<tWin/2+1;v++)
			SxxTemp[v]+=pow(com[v][0],2)+pow(com[v][1],2);//actually *dt^2 (but skipped /dt at input), actually /tWin/dt/runs (but done only when averaging)
	};
	fftw_destroy_plan(forward);
	fftw_cleanup();
};

void writeStats(int gen, spikeT *spiketrains)
{
	std::stringstream ss;
	ss<<gen;
	std::string str=ss.str();
	std::fstream dS((type+"Sxx"+str+".dat").c_str(), std::ios::out);

	for(long int v=0;v<tWin/2+1;v++)//in Hertz
		dS<<double(v)/dt/tWin*1000<<" "<<Sxx[v]*1000<<"\n";
	dS.close();


        std::fstream ST((type+"Spikes"+str+".dat").c_str(), std::ios::out);
        std::fstream RH((type+"RateHisto"+str+".dat").c_str(), std::ios::out);
	for(long int run=0;run<runs;run++)
	{
		ST<<run<<" ";
		for(long int hi=0; hi<spiketrains[run].data.size(); hi++)
			ST<<spiketrains[run].data.at(hi)<<" ";
		ST<<"\n";

		RH<<run<<" "<<1e3*double(spiketrains[run].data.size())/tWin/dt<<"\n";
	};
	ST.close();
	RH.close();

//--------------CV and firing rate
	double ISIav[2]={};//
	double ISIarray[runs][2]={};//record local ISI moments
	double rateAv=0;//record averages
	double rateArray[runs]={};//record local rates
	double CVav=0;//record averages
	double CVarray[runs]={};//record local CVs
	long int zeroIntervals=0;
	double spikeCount[2]={};

	double ISIhis[tWin+1]={};//ISI histogram, all other statistics cross-neuronal and written straight to file
	double nISIOverall=0;//total number of ISIs

	for(long int nPost=0; nPost<runs; nPost++)
	{
		double ISIlocal=0,ISI2local=0;// for each neuron
		long int STsize=spiketrains[nPost].data.size();
		spikeCount[0]+=STsize;
		spikeCount[1]+=pow(STsize,2);
		if(STsize>1)
			nISIOverall+=STsize-1;//ISI count
		for(long int hi=1; hi<STsize; hi++)//browse through spike train
		{
			long int dummyISI=spiketrains[nPost].data.at(hi)-spiketrains[nPost].data.at(hi-1);
			ISIlocal+=dummyISI;
			ISI2local+=dummyISI*dummyISI;
			ISIhis[dummyISI]++;
		};
		if(STsize>1)
		{
			ISIlocal=ISIlocal/(STsize-1);
			ISI2local=ISI2local/(STsize-1);
			if(ISI2local-ISIlocal*ISIlocal<0)
				CVarray[nPost]=0;
			else
				CVarray[nPost]=sqrt(ISI2local-ISIlocal*ISIlocal)/ISIlocal;//sonst 0
			CVav+=CVarray[nPost];
		}
		else
		{
			CVarray[nPost]+=-1;
			zeroIntervals++;
		};
		rateArray[nPost]=double(STsize)/tWin;//still dimensionless	
		rateAv+=double(STsize)/tWin;
		ISIarray[nPost][0]=ISIlocal;
		ISIarray[nPost][1]=ISI2local;
		ISIav[0]+=ISIlocal;
		ISIav[1]+=ISI2local;
	};
	rateAv=rateAv/runs;
	CVav=double(CVav)/(runs-zeroIntervals);
	ISIav[0]=ISIav[0]/(runs-zeroIntervals);
	ISIav[1]=ISIav[1]/(runs-zeroIntervals);
	spikeCount[0]=spikeCount[0]/runs;
	spikeCount[1]=spikeCount[1]/runs;
			
//------------------SCCs across ensemble-----
	double rhoAv[kmax]={};
	double rhoMix[kmax]={};
	long int intervalsM[kmax]={};
	double rhoAv2[kmax]={};//through strict ensemble average
	long int zeroIntervals2[kmax]={};
	for(long int nPost=0; nPost<runs; nPost++)
	{
		double rhoMix2[kmax]={};
		long int ST=spiketrains[nPost].data.size();
		for(int k=0; k<kmax; k++)
		{
			if(k+2<=ST)
			{	
				intervalsM[k]+=ST-k-1;			
				for(long int nj=1; nj<ST-k; nj++)
				{
					rhoMix[k]+=((spiketrains[nPost].data.at(nj)-spiketrains[nPost].data.at(nj-1))-ISIav[0])*((spiketrains[nPost].data.at(nj+k)-spiketrains[nPost].data.at(nj+k-1))-ISIav[0]);
					rhoMix2[k]+=((spiketrains[nPost].data.at(nj)-spiketrains[nPost].data.at(nj-1))-ISIarray[nPost][0])*((spiketrains[nPost].data.at(nj+k)-spiketrains[nPost].data.at(nj+k-1))-ISIarray[nPost][0]);
				};
				rhoAv2[k]+=(rhoMix2[k]/(ST-k-1))/(ISIarray[nPost][1]-pow(ISIav[0],2));
			}
			else
				zeroIntervals2[k]++;
		};
	};
	for(int k=0; k<kmax; k++)
	{
		rhoAv[k]=(rhoMix[k]/intervalsM[k])/(ISIav[1]-pow(ISIav[0],2));
		rhoAv2[k]=rhoAv2[k]/(runs-zeroIntervals2[k]);
	};
//write data for given generation
	double FanoFactor=(spikeCount[1]-pow(spikeCount[0],2))/spikeCount[0];
	data<<"gen "<<gen<<" rate="<<rateAv/dt<<" CV="<<CVav<<" Fano factor="<<Sxx[0]/(rateAv/dt)<<"="<<FanoFactor<<"\n";
	data<<"|rho|"<<"\n";
	for(int k=0; k<kmax; k++)
		data<<"k="<<k<<" : "<<rhoAv[k]<<" "<<rhoAv2[k]<<"\n";
	data<<"FRACTION OF ZERO SPIKING: "<<double(zeroIntervals)/runs<<"\n";
	data<<"\n";
}

void writeParameters()
{
	std::fstream param((type+"parameters.dat").c_str(), std::ios::out);
	param<<"leak="<<Leak<<"\n";
	param<<"Ce="<<Ce<<"\n";
	param<<"Ci="<<Ci<<"\n\n";

	param<<"g="<<g<<"\n";
	param<<"vTh="<<vTh<<"\n";
	param<<"mu="<<mu<<"\n";
	param<<"tauM="<<tau<<"\n";
	param<<"tauS="<<tauS<<"\n";
	param<<"tRef="<<rp<<"\n";
	param<<"J="<<J<<"\n";
	param<<"rate of initial Poisson process="<<r0<<"\n";
	param<<"initial power at S(0), decay to r0="<<r0b<<"\n\n";

	param<<"#generations="<<genMax<<"\n";
	param<<"#runs="<<runs<<"\n";
	param<<"dt="<<dt<<"\n";
	param<<"tStart="<<tStart<<"\n";
	param<<"tMax="<<tMax<<"\n";
	param<<"max SCC lag="<<kmax<<"\n";
	param.close();
};

void boxmuller(double *gaussR, double *gaussI)//generate two Gaussian variables with unit variance
{
	static double V1, V2, S;
	do
	{
//		double U1 = (double)rand() / RAND_MAX;
//		double U2 = (double)rand() / RAND_MAX;
		double U1=mtrand1();	
		double U2=mtrand1();	
		V1 = 2 * U1 - 1;
		V2 = 2 * U2 - 1;
		S = V1 * V1 + V2 * V2;
	}
	while(S >= 1 || S == 0);
	*gaussR = V1 * sqrt(-2 * log(S) / S);
	*gaussI = V2 * sqrt(-2 * log(S) / S);
}; 
	
double generateNoise(int gen,double *noiseterm)
{
	fftw_plan back=fftw_plan_dft_c2r_1d(tMax, eta, noiseterm, FFTW_ESTIMATE);
	double gaussR=0,gaussI=0;//real and imaginary part
//0th and (tMax/2)-th bin must be real (Hermitian symmetry)
	long int v2=0;
	boxmuller(&gaussR,&gaussI);
	eta[0][0]=gaussR*sqrt(Sxx[v2]*(Ce+Ci*g*g)*tMax*dt)*J*tau;//sole real part has to account for all the variance
	v2=(long int)(tWin/2+0.5);
	eta[tMax/2][0]=gaussI*sqrt(Sxx[v2]*(Ce+Ci*g*g)*tMax*dt)*J*tau;//sole real part has to account for all the variance

	for(long int v=1;v<tMax/2;v++)// otherwise two random numbers for each frequency bin
	{
		boxmuller(&gaussR,&gaussI);//... both normally distributed around 0 with SD 1
		v2=(long int)(v*double(tWin)/(tMax)+0.5);//interpolate binning; after backtransform we want to have time series of tMax steps
		eta[v][0]=gaussR*sqrt(Sxx[v2]*(Ce+Ci*g*g)/2*tMax*dt)*J*tau;
		eta[v][1]=gaussI*sqrt(Sxx[v2]*(Ce+Ci*g*g)/2*tMax*dt)*J*tau;
	};
//backtransform
	fftw_execute(back);
	fftw_destroy_plan(back);
	fftw_cleanup();
        for(long int tStep=0;tStep<tMax;tStep++) //in ms
		noiseterm[tStep]=noiseterm[tStep]/tMax/dt;//rescaling with 1/tMax/dt=df necessary for DFT backtransform to emulate df in integration
};
//#endif
