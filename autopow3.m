function Sa=autopow3(f,Hs,Tz)
% syntax: function Sa=autopow3(f,Hs,Tz)
% Autopower spectral density function of a random sea
% Input:
%   f: frequency (Hz)
%   Hs: significant wave height (m)
%   Tz: zero crossing period (s)
% Output:
%   Sa: autopower spectral density (m^2*s)

% Pierson-Moskowitz spectrum
Sa=Hs.^(2)./(4.*pi.*Tz.^(4).*f.^(5)).*exp(-1./pi.*(f.*Tz).^(-4));
