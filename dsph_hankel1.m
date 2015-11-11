function dj = dsph_hankel1(n, x)
%DSPH_HANKEL1 Spherical hankel function derivative of the first kind.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DSPH_HANKEL1.M - 15/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dj = dsph_besselj(n,x) + 1i*dsph_bessely(n,x);

end
