function W_Delta = generate_W_Delta(code, fs, fc, Delta_tau, Delta_b)
%--------------------------------------------------------------------------
%        AUTHOR: Corentin Lubeigt
%       CREATED: 05/11/2020
%   
%   DESCRIPTION: This function generates the interference matrix W_Delta
%   that is involved in the Fisher Information Matrix of the collection
%   with one antenna of a signal and a single reflection in an additive 
%   white Gaussian case and under band-limited assumption.
%   
%         NOTE: From this general function, it is easy to get the W matrix
%   which corresponds to the case with a unique signal [1]. The W matrix is
%   also used to compute the Fisher Information matrix. To obtain the W 
%   matrix from this function, just set Delta_tau and Delta_b to 0.
%
%        INPUTS: code       [] Nx1 vector containing the code used
%                fs         [] reduced sampling frequency
%                fc         [] reduced carrier frequency
%                Delta_tau  [] reduced time delay difference between path 0 and
%                              path 1
%                Delta_b    [] Doppler dilatation difference between path 0
%                              and path 1
%
%       OUTPUTS: W_Delta    [] a 3x3 matrix that characterises the
%                              interference between a signal and a
%                              time-delay, frequency shifted replica of it
%
%    REFERENCES: [1] 2020 [Medina et al] Compact CRB for delay, Doppler, 
%                    and phase estimation - application to GNSS SPP and RTK
%                    performance characterisation
%                    DOI: 10.1049/iet-rsn.2020.0168
%                [2] 2020 [Lubeigt et al] Joint Delay-Doppler Estimation 
%                    Performance in a Dual Source Context
%                    DOI:10.3390/rs12233894
%--------------------------------------------------------------------------

%% I - Initialisation

N = length(code);
wc = 2*pi*fc;
Ts = 1/fs;

%% II - Priori calculations

% 1- Vectors involved in the W_Delta components

% time vector
dt = -(1-N:N-1) + Delta_tau/Ts;

% 0 derivative
v0 = sinc(dt);

% first derivative
v1 = (cos(pi*dt) - v0);
v1(not(not(dt))) = v1(not(not(dt)))./dt(not(not(dt)));

% second derivative
v2 = 2*v1;
v2(not(not(dt))) = v2(not(not(dt)))./dt(not(not(dt)));
v2 = (pi^2)*v0 + v2;
v2(not(dt)) = pi^2/3;

% 2- series of matrices used to compute the W_Delta components

% Matrix D (due to first derivatives)
diagD = 1:N;

% Matrix U (due to different Doppler between signals)
if (Delta_b ~= 0)
    diagU = exp(-1j*2*pi*Ts*fc*Delta_b*(1:N));
else
    diagU = ones(1,N);
end

% 3- Intermediate combination of vectors and matrices to avoid huge
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

WD11 = (1/fs)*(WD_cU*WD_VD0c);
WD12 = (1/fs^2)*(WD_cDU*WD_VD0c);
WD13 = -(WD_cU*WD_VD1c) + (1j*wc*Delta_b/fs)*(WD_cU*WD_VD0c);
WD21 = (1/fs^2)*(WD_cU*WD_VD0Dc);
WD22 = (1/fs^3)*(WD_cDU*WD_VD0Dc);
WD23 = -(1/fs)*(WD_cU*WD_VD1Dc) + (1j*wc*Delta_b/fs^2)*(WD_cU*WD_VD0Dc);
WD31 = WD_cU*WD_VD1c;
WD32 = (1/fs)*(WD_cDU*WD_VD1c);
WD33 = fs*(WD_cU*WD_VD2c) + 1j*wc*Delta_b*(WD_cU*WD_VD1c);

%% IV - Computation of the matrix

W_Delta = [WD11, WD12, WD13;
           WD21, WD22, WD23;
           WD31, WD32, WD33];

% Note: here one still needs to multiply W_Delta with the complex
% exponential that contains the phase difference.

end % function generate_W_Delta