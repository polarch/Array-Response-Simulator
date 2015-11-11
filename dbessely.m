function dj = dbessely(n, x)
%DBESSELY Bessel function derivative of the second kind.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DBESSELY.M - 15/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n==0
    dj = -bessely(1,x);
else
    dj = (bessely(n-1,x) - bessely(n+1,x))/2;
end
