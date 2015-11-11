function [h_mic, H_mic] = sphericalScatterer(mic_dirs_rad, src_dirs_rad, R, N_order, N_filt, fs)
%SPHERICALSCATTERER Compute the pressure due to a spherical scatterer
%   SPHERICALSCATTERER computes the impulse responses of the pressure
%   measured at some points in the field with a spherical rigid scatterer
%   centered at the origin and due to incident plane waves.
%
%   mic_dirs_rad =  [mic_azi1 mic_elev1 mic_dist1; ...] directions and
%                   distances of the measurement positions (in rads,
%                   meters)
%
%   src_dirs_rad =  [src_azi1 src_elev1; ...] DOAs of plane waves for
%                   evaluation (in rads)
%
%   R            =  scatterer radius in m
%   N_order      =  maximum order of spherical approximation
%   fs           =  sample rate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPHERICALSCATTERER.M - 12/5/2014
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the frequency-dependent part of the microphone responses (radial
% dependence)
f = (0:N_filt/2)'*fs/N_filt;
c = 343;
kR = 2*pi*f*R/c;
N_mic = size(mic_dirs_rad,1);
N_pw = size(src_dirs_rad,1);
if any(mic_dirs_rad(:,3)<R)
    error('the distance of the measurement point cannot be less than the radius')
end
% check if all microphones at same radius
ALLATSAMERAD = all((mic_dirs_rad(:,3) - mic_dirs_rad(1,3))==0);
% if all microphones ar same radius, simplify
if ALLATSAMERAD
    b_N = zeros(N_filt/2+1, N_order+1);
    r = mic_dirs_rad(1, 3);
    kr = 2*pi*f*r/c;
    
    for n=0:N_order
        
        jn = sph_besselj(n, kr);
        jnprime = dsph_besselj(n, kR);
        hn = sph_hankel2(n, kr);
        hnprime = dsph_hankel2(n, kR);
        
        b_N(:, n+1) = (2*n+1) * 1i^n * (jn-(jnprime./hnprime).*hn);
    end
else
    
    b_N = zeros(N_filt/2+1, N_order+1, N_mic);
    for nm=1:N_mic
        r = mic_dirs_rad(nm, 3);
        kr = 2*pi*f*r/c;
        
        for n=0:N_order
            
            jn = sph_besselj(n, kr);
            jnprime = dsph_besselj(n, kR);
            hn = sph_hankel2(n, kr);
            hnprime = dsph_hankel2(n, kR);
            
            b_N(:, n+1, nm) = (2*n+1) * 1i^n * (jn-(jnprime./hnprime).*hn);
        end
    end
end
% Avoid NaNs for very high orders, instead of (very) very small values
b_N(isnan(b_N)) = 0;

% Compute angular-dependent part of the microphone responses
% unit vectors of DOAs and microphones
H_mic = zeros(N_filt/2+1, N_mic, N_pw);
for np=1:N_pw
    azi0 = src_dirs_rad(np,1);
    elev0 = src_dirs_rad(np,2);
    azi = mic_dirs_rad(:,1);
    elev = mic_dirs_rad(:,2);
    cosAlpha = sin(elev)*sin(elev0)+cos(elev)*cos(elev0).*cos(azi-azi0);
    
    P_N = zeros(N_order+1, N_mic);
    for n=0:N_order
        tempLeg = legendre(n, cosAlpha);
        P_N(n+1, :) = tempLeg(1,:);
    end
    % Accumulate across orders
    if ALLATSAMERAD
        H_mic(:,:,np) = b_N * P_N;
    else
        for nm=1:N_mic
            H_mic(:,nm,np) = b_N(:,:,nm) * P_N(:,nm);
        end
    end
end


%handle Nyquist for real impulse response
tempH_mic = H_mic;
tempH_mic(end,:,:) = abs(H_mic(end,:,:));
% conjugate ifft and fftshift for causal IR
tempH_mic = [tempH_mic; conj(tempH_mic(end-1:-1:2,:,:))];
h_mic = fftshift(ifft(tempH_mic),1);

end
