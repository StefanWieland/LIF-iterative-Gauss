#include "functionsGauss.h"

int main()
{
        time_t t1,t2;(void)time(&t1);	//initializing time
	writeParameters();		//writing system parameters
	data.open((type+".dat").c_str(), std::ios::out);//general data
	std::fstream RM((type+"Sxx0ReturnMap.dat").c_str(), std::ios::out);

	for(int v=0;v<tWin/2+1;v++)//input 0th generation, S(0) included!
		Sxx[v]=2*(r0b-r0)/(1+exp(10*double(v-1)/dt/tWin))+r0;			
	
	double ReturnMapS0=Sxx[0];//first coordinate of Sxx(0) return map
	double ratePrev=0;//rate of previous generation
	double totalSpikes=r0*runs*dt*tWin;//expected total spike count over all runs

	for(int gen=1;gen<=genMax;gen++)//loop over generations
	{
		for(int v=0;v<tWin/2+1;v++)//initializing Sxx
			SxxTemp[v]=0;

		ratePrev=totalSpikes/runs/dt/tWin;//averaged firing rate of gen i-1 -> mean input current of gen i
		totalSpikes=0;//initializing to measure current average firing rate	

		spikeT spiketrains[runs]={};//initializing	
		for(long int run=0;run<runs;run++)//loop over runs
		{		
			generateNoise(gen,noiseterm);//generate noise vector with tMax time bins
			double gaussR=0,gaussI=0;
			boxmuller(&gaussR,&gaussI);

			neuronClass neuron;//initializing neuron

			for(long int tD=0;tD<tMax;tD++)//integrating LIF
			{
				if(neuron.state[0]>=vTh)//threshold condition	
				{
					if(tStart <= tD)//if transient over...
						spiketrains[run].data.push_back(tD);//...recording spike
					neuron.state[0]=vR;//setting voltage back to vR
					neuron.lastSpike=tD;//recording spiking time
				};

				if((tD-neuron.lastSpike)*dt>=rp)//refractory condition
					neuron.state[0]=neuron.state[0]+dt/tau*(-Leak*neuron.state[0]+J*(Ce-g*Ci)*tau*ratePrev+mu+noiseterm[tD]);
			};//end of time window
			totalSpikes+=double(spiketrains[run].data.size());//increasing global spike count 
			
		};//end of run

		getSxx(double(totalSpikes)/runs/dt/tWin,SxxTemp,spiketrains);//generate new PSD to sample from in next generation
		for(long int v=0;v<tWin/2+1;v++)//in ms
			Sxx[v]=SxxTemp[v]/runs/tWin/dt;
		std::cout.precision(15);
		RM<<ReturnMapS0<<" "<<Sxx[0]<<"\n";
		ReturnMapS0=Sxx[0];

		if(gen%10==0)
			writeStats(gen,spiketrains);
		std::cout<<"generation:"<<gen<<" rate="<<1e3*totalSpikes/runs/dt/tWin<<"\n";
	
 	};
	data.close();
	RM.close();
}
