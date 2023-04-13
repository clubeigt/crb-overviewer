function W_Delta = generate_W_Delta(code, Fs, Fc, Delta_tau, Delta_b)
%--------------------------------------------------------------------------
%         USAGE: generate_W_Delta(code, Fs, Fc, Delta_tau, Delta_b)
%
%        AUTHOR: Corentin Lubeigt
%       CREATED: 05/11/2020
%   
%   DESCRIPTION: This function generates the interference matrix W_Delta
%   that is involved in the Fisher Information Matrix of the collection
%   with one antenna of a signal and a single reflection in a Gaussian
%   case and under band-limited assumption.
%   
%         NOTE: From this general function, it is easy to get the W matrix
%         which corresponds to the case with a unique signal [1]. The W matrix
%         is also used to compute the Fisher Information matrix. To obtain
%         the W matrix from this function, just set Delta_tau and Delta_b
%         to 0.
%
%        INPUTS: code       []   1xnSamp vector containing the code used
%                Fs         [Hz] sampling frequency (can be a reduced frequency)
%                Fc         [Hz] carrier frequency (can be a reduced frequency)
%                Delta_tau  [s]  time delay difference between path 0 and
%                                path 1 (can be a reduced time (in terms of
%                                C/A chips for instance))
%                Delta_b    []   doppler drift difference between path 0
%                                and path 1
%
%       OUTPUTS: W_Delta    []   a 3x3 matrix that characterises the
%                                interference between a signal and itself 
%
%    REFERENCES: [1] 2020 [Medina et al] Compact CRB for delay, Doppler, and phase estimation - application to GNSS SPP and RTK performance characterisation
%                [2] 2020 [Lubeigt et al] Joint Delay-Doppler Estimation Performance in a Dual Source Context
%--------------------------------------------------------------------------

%% I - Initialisation

nSamp = length(code);
wc = 2*pi*Fc;
Ts = 1/Fs;

%% II - Priori calculations

% 1) Vectors involved in the W_Delta components

% 1a- time vector
dt = -(1-nSamp:nSamp-1) + Delta_tau/Ts;

% 1b- 0 derivative
v0 = sinc(dt);

% 1c- first derivative
v1 = (cos(pi*dt) - v0);
v1(not(not(dt))) = v1(not(not(dt)))./dt(not(not(dt)));

% 1d- second derivative
v2 = 2*v1;
v2(not(not(dt))) = v2(not(not(dt)))./dt(not(not(dt)));
v2 = (pi^2)*v0 + v2;
v2(not(dt)) = pi^2/3;

% 2) series of matrices used to compute the W_Delta components

% 2a- Matrix D (due to first derivatives)
diagD = 1:nSamp;

% 2b- Matrix U (due to different doppler between signals)
if (Delta_b ~= 0)
    diagU = exp(-1j*2*pi*Ts*Fc*Delta_b*(1:nSamp));
else
    diagU = ones(1,nSamp);
end

% 3) Intermediate combination of vectors and matrices to avoid huge
% matrices handling

WD_cD    = code'.*diagD;
WD_Dc    = WD_cD';
if ( Delta_b ~= 0)
    WD_cU    = code'.*diagU;
    WD_cDU   = WD_cD.*diagU;
else
    WD_cU    = code';
    WD_cDU   = WD_cD;    
end

WD_VD0c = fftshift(ifft(fft(code,length(v0)).*(fft(v0)')));
WD_VD0c = WD_VD0c(1:length(code));

WD_VD0Dc = fftshift(ifft(fft(WD_Dc,length(v0)).*(fft(v0)')));
WD_VD0Dc = WD_VD0Dc(1:length(WD_Dc));

WD_VD1c = fftshift(ifft(fft(code,length(v1)).*(fft(v1)')));
WD_VD1c = WD_VD1c(1:length(code));

WD_VD1Dc = fftshift(ifft(fft(WD_Dc,length(v1)).*(fft(v1)')));
WD_VD1Dc = WD_VD1Dc(1:length(WD_Dc));

WD_VD2c = fftshift(ifft(fft(code,length(v2)).*(fft(v2)')));
WD_VD2c = WD_VD2c(1:length(code));

%% III - Evaluation of each components of W_Delta

WD11 = (1/Fs)*(WD_cU*WD_VD0c);
WD12 = (1/Fs^2)*(WD_cDU*WD_VD0c);
WD13 = -(WD_cU*WD_VD1c) + (1j*wc*Delta_b/Fs)*(WD_cU*WD_VD0c);
WD21 = (1/Fs^2)*(WD_cU*WD_VD0Dc);
WD22 = (1/Fs^3)*(WD_cDU*WD_VD0Dc);
WD23 = -(1/Fs)*(WD_cU*WD_VD1Dc) + (1j*wc*Delta_b/Fs^2)*(WD_cU*WD_VD0Dc);
WD31 = WD_cU*WD_VD1c;
WD32 = (1/Fs)*(WD_cDU*WD_VD1c);
WD33 = Fs*(WD_cU*WD_VD2c) + 1j*wc*Delta_b*(WD_cU*WD_VD1c);

%% IV - Computation of the matrix

W_Delta = [WD11, WD12, WD13;
           WD21, WD22, WD23;
           WD31, WD32, WD33];

% Note: here one still needs to multiply W_Delta with the complex
% exponential that contains the phase difference.

end % function generate_W_Delta