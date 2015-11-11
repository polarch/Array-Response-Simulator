function dj = dsph_hankel2(n, x)
%DSPH_HANKEL2 Spherical hankel function derivative of the second kind.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DSPH_HANKEL2.M - 15/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dj = dsph_besselj(n,x) - 1i*dsph_bessely(n,x);

end
