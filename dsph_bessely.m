function dj = dsph_bessely(n, x)
%DSPH_BESSELY Spherical bessel function derivative of the second kind.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DSPH_BESSELY.M - 15/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dj = 1/(2*n+1)*(n*sph_bessely(n-1,x) - (n+1)*sph_bessely(n+1,x));

end

