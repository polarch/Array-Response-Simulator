function b_N = cylModalCoeffs(N, kr, arrayType)
%CYLMODALCOEFFS Modal coefficients for rigid or open cylindrical array
%
%   N: maximum order
%   kr: wavenumber-radius product
%   arrayType: {'open','rigid'} open for open array of omnidirectional
%              sensors, rigid for sensors mounted on a rigid baffle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CYLMODALCOEFFS.M - 15/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b_N = zeros(length(kr), N+1);

for n=0:N
    
    if (isequal(arrayType, 'open'))
        b_N(:, n+1) = 1i^n * besselj(n, kr);
        
    elseif (isequal(arrayType,'rigid'))
        jn = besselj(n, kr);
        jnprime = dbesselj(n, kr);
        hn = besselh(n, 2, kr);
        hnprime = dhankel2(n, kr);
        
        temp = 1i^n * (jn-(jnprime./hnprime).*hn);
        if n==0
            temp(kr==0) = 1;
        else
            temp(kr==0) = 0;
        end
        b_N(:, n+1) = temp;
        
    else
        error('Wrong array type')
    end
end

% Avoid NaNs for very high orders, instead of very small values
b_N(isnan(b_N)) = 0;
