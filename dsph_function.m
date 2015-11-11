function df_x = dsph_function(n, x, funcName)
%DSPH_FUNCTION Derivative of spherical Bessel and Hankel functions

switch funcName
    case 'besselj'
    df_x = (n./x).*sph_besselj(n, x) - sph_besselj(n+1, x);

    case 'bessely'
    df_x = (n./x).*sph_bessely(n, x) - sph_bessely(n+1, x);        
        
    case 'hankel1'
    df_x = (n./x).*sph_hankel1(n, x) - sph_hankel1(n+1, x);
        
    case 'hankel2'
    df_x = (n./x).*sph_hankel2(n, x) - sph_hankel2(n+1, x);
end

end
