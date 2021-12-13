 function [t,eta]=wave(Hs,Tz,N,deltat,seed)
% syntax: function [t,eta]=wave(Hs,Tz,N,deltat,seed)
% Wave simulation based on Fourier coefficients
% Output:
% t: time (s)
% eta: surface elevation (m)
% Input:
% Hs: significant wave height (m)
% Tz: zero crossing period (s)
% N: number of time points (including zero); N must be a power of 2
% deltat: time step (s)
% seed: used to initialize random number generator

% time vector
t=[0:N-1]'*deltat;

% period
T=N*deltat;

% frequency step
deltaf=1/T;

% discretized frequencies
k=[1:N/2-1]';
f=k.*deltaf;

% autopower spectral density (one-sided)
Sa=autopow3(f,Hs,Tz);

% (half) diagonal of the weight matrix
dW=sqrt(Sa/T);
% diagonal of weight matrix
W=[dW;dW];

% vector of unit variance normal random numbers
rng(seed);
r=randn(N-2,1);

% multiplication: W*r
% Note: the matrix multiplication is replaced by an inner product
%       since W is a vector containing the diagonal of the matrix
abk=W.*r;

% complex notation
i=sqrt(-1);
C=abk(1:N/2-1)-i*abk(N/2:N-2);
C=1/2*[0;C;0;rot90(C')];

% inverse FFT
eta=N.*ifft(C);

if any(abs(imag(eta)) >= 1e-7*abs(eta) & abs(imag(eta)) >= 1e-12)
  max(abs(eta))
  max(imag(eta))
  error('imag too large')
end
eta=real(eta);
