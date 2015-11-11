function dj = dbesselj(n, x)
%DBESSELJ Bessel function derivative of the first kind.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DBESSELJ.M - 15/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n==0
    dj = -besselj(1,x);
else
    dj = (besselj(n-1,x) - besselj(n+1,x))/2;

end
