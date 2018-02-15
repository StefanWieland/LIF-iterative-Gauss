#pragma once
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include "MersenneTwister.h"
#include <string>
#include <fftw3.h>
#include <sstream>
#include <complex>

std::string type="test";//data path 

const long int Ce=1000;//number of excitatory presynaptic neurons
const long int Ci=250;//number of inhibitory presynaptic neurons

const double Leak=1;//size of leak term in membran voltage dynamics
const double g=4;//relative strength of inhibitory synaptic amplitude in mV
const double vTh=20;//firing threshold value in mV
const double vR=10; //voltage reset value in mV
const double mu=1.5*20; //constant input
const double tau=20;//membrane constant in ms
const double tauS=0;//synaptic time constant in ms
const double rp=2;//refractory period in ms
const double J=0.1; //excitatory synaptic amplitude in mV
const double r0=0.1;//power of bandpass-limited Gaussian white noise at first generation
const double r0b=0.1;//Sxx(0) in synaptic input of initial layer for different starting points of Sxx(0) return map

const int genMax=20;//number of generations/layers in the feed-forward network
const long int runs=1000;//number of averages taken
const double dt=1e-1;//width of time step of Euler scheme in ms
const long int tStart=2e3; //system transient
const long int tMax=7e3;//number of time steps
const long int tWin=tMax-tStart;//time window of sampling

const int kmax=11;//maximum lag of serial correlation coefficient





