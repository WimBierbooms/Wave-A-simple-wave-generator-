function H=waveheight(eta)
% syntax: function H=waveheight(eta)
% Holthuijsen: wave defined by zero crossings
% H: wave height (m) 
% eta: surface elevation (m)  

N=1:length(eta)-1;

% indices of zero upcrossings
N0=find(eta(N)<0 & eta(N+1)>0);

% wave heights (from trough to crest)
for i=1:length(N0)-1
   H(i)=max(eta(N0(i):N0(i+1)))-min(eta(N0(i):N0(i+1)));
end
