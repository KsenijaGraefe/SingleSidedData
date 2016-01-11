clear all
close all
%system matrix
SM_particle =  hdf5read('Data\SystemResponse.h5','systemResponseFrequencies');
%system matrix withou SPION
SM_empty_1 =   hdf5read('Data\SystemResponse_empty.h5','systemResponseFrequencies');

%loading measuremnts
particle2 = hdf5read('Data\Measurement.h5','frequencies');
% phantom measurement
p = 6;
% empty meausrement
e = 7;

%S: system matrix, u: voltag signal, c: particle concentration
addpath ('Functions') 
[S,u,c]=Reko(SM_particle,SM_empty_1,particle2,particle2,e,p);