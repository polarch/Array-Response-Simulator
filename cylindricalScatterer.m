function [h_mic, H_mic] = cylindricalScatterer(mic_dirs_rad, src_azis_rad, R, N_order, N_filt, fs)
%CYLINDRICALSCATTERER Compute the pressure due to a spherical scatterer
%   CYLINDRICALSCATTERER computes the impulse responses of the pressure
%   measured at some points in the field with a cylindrical rigid scatterer
%   centered at the origin and due to incident plane waves.
%
%   mic_dirs_rad =  [mic_azi1 mic_dist1; ...] directions and
%                   distances of the measurement positions (in rads,
%                   meters)
%
%   src_azi_rad =  [src_azi1; ...] DOAs of plane waves for evaluation (in rads)
%
%   R            =  scatterer radius in m
%   N_order      =  maximum order of spherical approximation
%   fs           =  sample rate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CYLINDRICALSCATTERER.M - 12/5/2014
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the frequency-dependent part of the microphone responses (radial
% dependence)
f = (0:N_filt/2)'*fs/N_filt;
c = 343;
kR = 2*pi*f*R/c;
N_mic = size(mic_dirs_rad,1);
N_pw = size(src_azis_rad,1);
if any(mic_dirs_rad(:,2)<R)
    error('the distance of the measurement point cannot be less than the radius')
end
% check if all microphones at same radius
ALLATSAMERAD = all((mic_dirs_rad(:,2) - mic_dirs_rad(1,2))==0);
% if all microphones ar same radius, simplify
if ALLATSAMERAD
    b_N = zeros(N_filt/2+1, N_order+1);
    r = mic_dirs_rad(1, 2);
    kr = 2*pi*f*r/c;
    
    for n=0:N_order
        
        jn = besselj(n, kr);
        jnprime = dbesselj(n, kR);
        hn = besselh(n, 2, kr);
        hnprime = dhankel2(n, kR);
        
        b_N(:, n+1) = 1i^n * (jn-(jnprime./hnprime).*hn);
    end
else
    
    b_N = zeros(N_filt/2+1, N_order+1, N_mic);
    for nm=1:N_mic
        r = mic_dirs_rad(nm, 2);
        kr = 2*pi*f*r/c;
        
        for n=0:N_order
            
            jn = besselj(n, kr);
            jnprime = dbesselj(n, kR);
            hn = besselh(n, 2, kr);
            hnprime = dhankel2(n, kR);
            
            b_N(:, n+1, nm) = 1i^n * (jn-(jnprime./hnprime).*hn);
        end
    end
end
% Avoid NaNs for very high orders, instead of (very) very small values
b_N(isnan(b_N)) = 0;

% Compute angular-dependent part of the microphone responses
% unit vectors of DOAs and microphones
H_mic = zeros(N_filt/2+1, N_mic, N_pw);
for np=1:N_pw
    azi0 = src_azis_rad(np);
    azi = mic_dirs_rad(:,1);
    angle = azi - azi0;
    
    C = zeros(N_order+1, N_mic);
    for n=0:N_order
        % circular expansion        
        if n==0
            C(n+1,:) = ones(size(angle));
        else
            C(n+1,:) = 2*cos(n*angle);
        end
    end
    % Accumulate across orders
    if ALLATSAMERAD
        H_mic(:,:,np) = b_N * C;
    else
        for nm=1:N_mic
            H_mic(:,nm,np) = b_N(:,:,nm) * C(:,nm);
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
