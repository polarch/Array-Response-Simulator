function [h_mic, H_mic] = simulateCylArray(N_filt, mic_dirs_rad, src_dirs_rad, arrayType, R, N_order, fs)
%SIMULATECYLARRAY Simulate the impulse responses of a cylindrical array.
%
%   SIMULATECYLARRAY computes the impulse responses of the microphones of a
%   cylindrical microphone array for the given directions of incident plane
%   waves. The array type can be either 'open' for omnidirectional
%   microphones in an open setup, or 'rigid' for omnidirectional
%   microphones mounted on a cylinder.
%
%   mic_dirs_rad :  [mic_azi1; ...] directions of microphone capsules
%
%   src_dirs_rad :  [src_azi1; ...] DOAs of plane waves for evaluation
%
%   arrayType    :  'open' or 'rigid'
%   R            :  array radius in m
%   N_order      :  maximum order of cylindrical expansion
%   fs           :  sample rate
%
%
% H_mic:    Output the frequency responses for M directions for the N 
%           microphones. Only half the FFT spectrum is returned (up to Nyquist).
%           Last dimension is the DOA.
% h_mic:    Output the IRs for M directions for the N microphones. Last 
%           dimension is the DOA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SIMULATECYLARRAY.M - 15/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the radial part of the microphone responses
f = (0:N_filt/2)'*fs/N_filt;
c = 343;
kR = 2*pi*f*R/c;
b_N = cylModalCoeffs(N_order, kR, arrayType);
%handle Nyquist for real impulse response
temp = b_N;
temp(end,:) = real(temp(end,:));
b_Nt = fftshift(ifft([temp; conj(temp(end-1:-1:2,:))]), 1);

% Compute angular-dependent part of the microphone responses
% unit vectors of DOAs and microphones
N_doa = size(src_dirs_rad,1);
N_mic = size(mic_dirs_rad,1);
h_mic = zeros(N_filt, N_mic, N_doa);
H_mic = zeros(N_filt/2+1, N_mic, N_doa);
for i=1:N_doa
    angle = mic_dirs_rad - src_dirs_rad(i);
    
    C = zeros(N_order+1, N_mic);
    for n=0:N_order
        % Jacobi-Anger expansion
        if n==0
            C(n+1,:) = ones(size(angle));
        else
            C(n+1,:) = 2*cos(n*angle);
        end
    end
    h_mic(:,:,i) = b_Nt * C;
    
    if nargout>1
        H_mic(:,:,i) = b_N * C;
    end
end

end
