function j = sph_besselj(n, x)
%SPH_BESSELJ Spherical bessel function of the first kind.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPH_BESSELJ.M - 15/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find zeros in the argument
idx_zero = (x==0);

j = sqrt(pi./(2*x)).*besselj(n+0.5, x);
if n==0
    j(idx_zero) = 1;
else
    j(idx_zero) = 0;
end

end

