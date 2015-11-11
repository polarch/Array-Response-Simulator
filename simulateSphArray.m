function [h_mic, H_mic] = simulateSphArray(N_filt, mic_dirs_rad, src_dirs_rad, arrayType, R, N_order, fs, dirCoeff)
%SIMULATESPHARRAY Simulate the impulse responses of a spherical array.
%
%   SIMULATESPHARRAY computes the impulse responses of the microphones of a
%   spherical microphone array for the given directions of incident plane
%   waves. The array type can be either 'open' for omnidirectional
%   microphones in an open setup, or 'rigid' for omnidirectional
%   microphones mounted on a sphere.
%
%   mic_dirs_rad :  [mic_azi1 mic_elev1; ...] directions of microphone
%                   capsules
%
%   src_dirs_rad :  [src_azi1 src_elev1; ...] DOAs of plane waves for
%                   evaluation
%
%   arrayType    :  'open', 'rigid', or 'directional'
%   R            :  array radius in m
%   N_order      :  maximum order of spherical approximation
%   fs           :  sample rate
%   dirCoeff     :  if array consists of directional microphones, then
%                   dirCoeff is the directivity coefficient from 0 to 1 (1
%                   for omni, 0.5 cardioid, 0 dipole)
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
% SIMULATESPHARRAY.M - 15/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<8
    dirCoeff = [];
end

% Compute the frequency-dependent part of the microphone responses (radial
% dependence)
f = (0:N_filt/2)'*fs/N_filt;
c = 343;
kR = 2*pi*f*R/c;
b_N = sphModalCoeffs(N_order, kR, arrayType, dirCoeff);
%handle Nyquist for real impulse response
temp = b_N;
temp(end,:) = real(temp(end,:));
b_Nt = fftshift(ifft([temp; conj(temp(end-1:-1:2,:))]), 1);

% Compute angular-dependent part of the microphone responses
% unit vectors of DOAs and microphones
N_doa = size(src_dirs_rad,1);
N_mic = size(mic_dirs_rad,1);
[U_mic(:,1), U_mic(:,2), U_mic(:,3)] = sph2cart(mic_dirs_rad(:,1), mic_dirs_rad(:,2), 1);
[U_doa(:,1), U_doa(:,2), U_doa(:,3)] = sph2cart(src_dirs_rad(:,1), src_dirs_rad(:,2), 1);
h_mic = zeros(N_filt, N_mic, N_doa);
H_mic = zeros(N_filt/2+1, N_mic, N_doa);
for i=1:N_doa
    cosangle = dot(U_mic, ones(N_mic,1)*U_doa(i,:), 2);
    
    P = zeros(N_order+1, N_mic);
    for n=0:N_order
        % The Legendre polynomial gives the angular dependency
        Pn = legendre(n, cosangle);
        P(n+1,:) = (2*n+1)/(4*pi) * Pn(1,:);
    end
    h_mic(:,:,i) = b_Nt * P;
    
    if nargout>1
        H_mic(:,:,i) = b_N * P;
    end
end

end
