function j = sph_hankel1(n, x)
%SPH_HANKEL1 Spherical hankel function of the first kind.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPH_HANKEL1.M - 15/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = sph_besselj(n,x) + 1i*sph_bessely(n,x);

end

