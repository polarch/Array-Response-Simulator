function f_x = sph_function(n, x, funcName)
%SPH_FUNCTION Spherical Bessel and Hankel function

switch funcName
    case 'besselj'
        if n==0
            f_x = sinc(x/pi);
        else
            idx_zero = (x==0);
            f_x = sqrt(pi./(2*x)).*besselj(n+0.5, x);
            f_x(idx_zero) = 0;
        end

    case 'bessely'
        f_x = sqrt(pi./(2*x)).*bessely(n+0.5, x);
        if any(x==0)
            warning('Zero argument for the Bessel function of the second kind results to -Inf')
        end
        
    case 'hankel1'
        f_x = sqrt(pi./(2*x)).*besselh(n+0.5, 1, x);
        
    case 'hankel2'
        f_x = sqrt(pi./(2*x)).*besselh(n+0.5, 2, x);
end

end
