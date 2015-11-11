function dj = dsph_besselj(n, x)
%DSPH_BESSELJ Spherical bessel function derivative of the first kind.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DSPH_BESSELJ.M - 15/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dj = 1/(2*n+1)*(n*sph_besselj(n-1,x) - (n+1)*sph_besselj(n+1,x));
% Teutch:
% dj = sph_besselj(n-1,x) - (n+1)/x*sph_besselj(n,x);

end
