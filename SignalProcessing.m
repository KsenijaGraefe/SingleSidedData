clear all
close all
%Laden der Systemmatrix, die mit einer Partikelprobe vermessen wurde
SM_particle =  hdf5read('dropbox\home\BetterSiSi\Documents\data\SystemResponse\SystemResponse.h5','systemResponseFrequencies');
%Laden der Systemmatrix, die ohne eine Partikelprobe vermessen wurde
SM_empty_1 =   hdf5read('dropbox\home\BetterSiSi\Documents\data\SystemResponse\SystemResponse_empty.h5','systemResponseFrequencies');

%Laden einer Phantommessung
particle2 = hdf5read('dropbox\home\BetterSiSi\Documents\data\SystemResponse\Measurement.h5','frequencies');
% falls mehrere Phantommessungen hintereinander abgespeichert sind:
% bestimmen welche benutzt wird
p=6;
%Laden einer Leermessung
particle3 = hdf5read('dropbox\home\BetterSiSi\Documents\data\SystemResponse\Measurement_empty.h5','frequencies');
% falls mehrere Messungen hintereinander abgespeichert sind:
% bestimmen welche benutzt wird
e=1;

%S: Systemmatrix, u: Spannungssignal, c: rekonstruierte
%Partikelkonzentration
[S,u,c]=Reko(SM_particle,SM_empty_1,particle2,particle3,e,p);