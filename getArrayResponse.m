function [micIRs, micTFs] = getArrayResponse(U_doa, R_mic, U_orient, fDir_handle, Lfilt, fs)
%GETARRAYRESPONSE Return array response of directional sensors.
%
%   GETARRAYRESPONSE computes the impulse responses of the microphones of
%   an open array of directional microphones, located at R_mic and with 
%   orientations U_orient, for the directions-of-incidence U_doa. Each
%   sensors directivity in defined by a function handle in the cell array 
%   fDir_handle.
%
% U_doa:    Mx3 matrix of vectors specifying the directions to evaluate the
%           impulse response
% R_mic:    Nx3 matrix of vectors specifying the positions of the microphones
% U_orient: Nx3 matrix of unit vectors specifying the orientation of the
%           microphones, or 1x3 unit vector if the microphones all have the
%           same orientation. If left empty or not defined then the
%           orientation of the microphones is assumed to be radial from the
%           origin.
% fDir_hadle:   1xN cell array of function handles describing the directivity
%               of each sensor with respect the sensor's orientation and the 
%               DOA. The function considers only axisymmetric directivities,
%               and the function should be of the form f(angle), where
%               angle is the angle between the sensors orientation and the
%               DOA. If a single function handle is provided, the same is
%               applied to all sensors.
% Lfilt:    Length of filter in samples.
% fs:       Sampling rate. If empty or undefined, 48kHz are used.
%
%
% micTFs:   Output the frequency responses for M directions for the N 
%           microphones. Only half the FFT spectrum is returned (up to Nyquist).
%           Last dimension is the DOA.
% micIRs:   Output the IRs for M directions for the N microphones. Last 
%           dimension is the DOA.
%
%   EXAMPLE - Simulate the response of a 3-microphone array of an
%   omnidirectional microphone, a first-order cardioid and a second-order
%   cardioid, with random locations and orientations, for front and side
%   incidence:
%   
%   Nmic = 3;
%   U_doa = [1 0 0; 0 1 0];
%   R_mic = rand(Nmic,3);
%   U_orient = rand(Nmic,3);
%   U_orient = U_orient./(sqrt(sum(U_orient.^2,2))*ones(1,3));   % convert to unit vectors
%
%   fdir_omni = @(angle) ones(size(angle));
%   fdir_card = @(angle) (1/2)*(1 + cos(angle));
%   fdir_card2 = @(angle) (1/2)^2 * (1 + cos(angle)).^2;
%   fDir_handle = {fdir_omni, fdir_card, fdir_card2};
%
%   [micIRs, micTFs] = getArrayResponse(U_doa, R_mic, U_orient, fDir_handle);
%
%   plot(micIRs(:,:,1))
%   legend('omni','1st-order cardioid','2nd-order cardioid')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GETARRAYRESPONSE.M - 15/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% DEFAULT ARGUMENTS %%%%%%%%%%
if nargin < 2
    error('Not enough input arguments.')
elseif nargin > 6
    error('Too many input arguments.')
end

Nmics = size(R_mic,1);
if size(R_mic, 2) ~= 3
    error(['Size of R_mic should be Nmicsx3, where the columns are the cartesian ' ...
        'coordinates of the vectors'])
end

Ndoa = size(U_doa,1);
if size(U_doa, 2) ~= 3
    error(['Size of U_src should be Ndoax3, where the columns are the cartesian ' ...
        'coordinates of the vectors'])
end

if ~exist('fs','var') || isempty(fs)
    fs = 48000;
end

% if no directivity coefficient is defined assume omnidirectional sensors
if ~exist('fDir_handle','var') || isempty(fDir_handle)
    fDir_handle = @(angle) 1;
end
if length(fDir_handle) == 1
    if iscell(fDir_handle)
        fDir_handle = repmat(fDir_handle,1,Nmics);
    elseif isa(fDir_handle,'function_handle')
        fDir_handle = {fDir_handle};
    	fDir_handle = repmat(fDir_handle,1,Nmics);
    end
elseif length(fDir_handle) ~= Nmics
    error('Size of fDir_handle should be 1xNmics')
end

% compute unit vectors of the microphone positions
normR_mic = sqrt(sum(R_mic.^2,2));
U_mic = R_mic./(normR_mic*ones(1,3));

% if no orientation is defined then assume that the microphones are
% oriented radially, similar to U_mic
if ~exist('U_orient','var') || isempty(U_orient)
    U_orient = U_mic;
elseif size(U_orient, 2) ~= 3
    error(['Size of U_orient should be 1x3, if all the micorphones are oriented '...
        'on the same direction, or Nmicsx3, where the columns are the cartesian '...
        'coordinates of the vectors'])
elseif isrow(U_orient)
    U_orient = U_orient*ones(Nmics,1);
elseif size(U_orient, 1) ~= Nmics
    error(['The number of rows of the orientation vectors should be equal '...
    'to the number of microphones'])
end


%%%%%%%%%% COMPUTATIONS %%%%%%%%%%

% speed of sound
c = 343;
% frequency vector
Nfft = Lfilt;
f = (0:Nfft/2+1)*fs/Nfft;
K = Nfft/2+1;

% unit vectors pointing to the evaluation points
U_eval = zeros(Ndoa, Nmics, 3);
U_eval(:,:,1) = U_doa(:,1) * ones(1, Nmics);
U_eval(:,:,2) = U_doa(:,2) * ones(1, Nmics);
U_eval(:,:,3) = U_doa(:,3) * ones(1, Nmics);

% computation of time delays and attenuation for each evaluation point to
% microphone, measured from the origin
tempR_mic(:,:,1) = ones(Ndoa, 1) * R_mic(:,1)';
tempR_mic(:,:,2) = ones(Ndoa, 1) * R_mic(:,2)';
tempR_mic(:,:,3) = ones(Ndoa, 1) * R_mic(:,3)';

tempU_orient(:,:,1) = ones(Ndoa, 1) * U_orient(:,1)';
tempU_orient(:,:,2) = ones(Ndoa, 1) * U_orient(:,2)';
tempU_orient(:,:,3) = ones(Ndoa, 1) * U_orient(:,3)';

% cos-angles between DOAs and sensor orientations
cosAngleU = dot(U_eval, tempU_orient, 3);
% d*cos-angles between DOAs and sensor positions
dcosAngleU = dot(U_eval, tempR_mic, 3);

% attenuation due to directionality of the sensors
B = zeros(Ndoa, Nmics);
for nm=1:Nmics
    B(:,nm) = fDir_handle{nm}(acos(cosAngleU(:,nm)));
end

% create TFs for each microphone
micTFs = zeros(K, Nmics, Ndoa);
for kk = 1:K
    omega = 2*pi*f(kk);
    tempTF = B .* exp(1i*(omega/c)*dcosAngleU);
    micTFs(kk,:,:) = tempTF.';
end

% create IRs for each microphone
micIRs = zeros(Nfft, Nmics, Ndoa);
for nd = 1:Ndoa
    tempTF = micTFs(:,:,nd);
    tempTF(end,:) = abs(tempTF(end,:));
    tempTF = [tempTF; conj(tempTF(end-1:-1:2,:))];
    micIRs(:,:,nd) = ifft(tempTF);
    micIRs(:,:,nd) = fftshift(micIRs(:,:,nd), 1);
end
