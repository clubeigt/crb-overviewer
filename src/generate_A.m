function A = generate_A(code, fs, fc, epsilon)
%--------------------------------------------------------------------------
%        AUTHOR: Corentin Lubeigt
%       CREATED: 02/21/2022
%   
%   DESCRIPTION: This function generates the A matrix as it is defined
%   in [1]. This matrix is involved in the computation of the misspecified 
%   Cramér-Rao bound for the estimation of the time-delay and the Doppler
%   frequency of a single source collected by a single antenna but in 
%   presence of a single multipath.
% 
%        INPUTS: code    [] Nx1 vector containing the considered code
%                fs      [] reduced sampling frequency
%                fc      [] reduced carrier frequency
%                epsilon [] 8x1 vector of parameters
%                           [tau_0, b_0, rho_0, phi_0, tau_1, b_1, rho_1, phi_1]^T
%
%       OUTPUTS: A       [] 8x8 matrix that corresponds to the matrix A
%                           involved in the computation of the Misspecified
%                           Cram'ér-Rao bound associated to the estimation 
%                           of a subset of the vector of parameters epsilon
%                           in the case of an additive white Gaussian 
%                           circular centered noise
%
%    REFERENCES: [1] 2023 [Lubeigt et al] Untangling first and second order
%                    statistics contributions in multipath scenarios
%                    DOI: 10.1016/j.sigpro.2022.108868
%--------------------------------------------------------------------------

%% I - Initialisation

wc = 2*pi*fc; % reduced pulsation

tau_0 = epsilon(1);
b_0   = epsilon(2);
rho_0 = epsilon(3);
phi_0 = epsilon(4);

tau_1 = epsilon(5);
b_1   = epsilon(6);
rho_1 = epsilon(7);
phi_1 = epsilon(8);

tau_pt = epsilon(9);
b_pt   = epsilon(10);
rho_pt = epsilon(11);
phi_pt = epsilon(12);

alpha_0  = rho_0*exp(1j*phi_0);
alpha_1  = rho_1*exp(1j*phi_1);
alpha_pt = rho_pt*exp(1j*phi_pt);

Delta_tau_0   = (tau_0 - tau_pt)/fs; % [chips]
Delta_tau_1   = (tau_1 - tau_pt)/fs; % [chips]

Delta_b_0     = b_0 - b_pt;
Delta_b_1     = b_1 - b_pt;

%% II - Q matrices

Q_tau = [    -alpha_pt*wc^2*b_pt^2,                  0, 0, -1j*2*alpha_pt*wc*b_pt,               0, alpha_pt;
                    1j*alpha_pt*wc, alpha_pt*wc^2*b_pt, 0,                      0,  1j*alpha_pt*wc,        0;
         1j*exp(1j*phi_pt)*wc*b_pt,                  0, 0,        -exp(1j*phi_pt),               0,        0;
                 -alpha_pt*wc*b_pt,                  0, 0,           -1j*alpha_pt,               0,        0];

Q_b   = [1j*alpha_pt*wc,    alpha_pt*wc^2*b_pt,              0, 0, 1j*alpha_pt*wc, 0;
                      0,                     0, -alpha_pt*wc^2, 0,              0, 0;
                      0, -1j*exp(1j*phi_pt)*wc,              0, 0,              0, 0;
                      0,           alpha_pt*wc,              0, 0,              0, 0];

Q_rho = [1j*exp(1j*phi_pt)*wc*b_pt,                      0, 0, -exp(1j*phi_pt), 0, 0;
                                  0, -1j*exp(1j*phi_pt)*wc, 0,               0, 0, 0;
                                  0,                     0, 0,               0, 0, 0;
                  1j*exp(1j*phi_pt),                     0, 0,               0, 0, 0];

Q_phi = [-alpha_pt*wc*b_pt,           0, 0, -1j*alpha_pt, 0, 0;
                         0, alpha_pt*wc, 0,            0, 0, 0;
         1j*exp(1j*phi_pt),           0, 0,            0, 0, 0;
                 -alpha_pt,           0, 0,            0, 0, 0];

%% III - W_m matrices

W_m_0  = generate_W_m(code, fs, fc, Delta_tau_0, Delta_b_0, b_0);
W_m_1  = generate_W_m(code, fs, fc, Delta_tau_1, Delta_b_1, b_1);
W_m_pt = generate_W_m(code, fs, fc, 0, 0, b_pt);

W_m = conj(alpha_0)*W_m_0 + conj(alpha_1)*W_m_1 - conj(alpha_pt)*W_m_pt;

%% IV - Computation of the A matrix

A = [Q_tau, Q_b, Q_rho, Q_phi]*kron(eye(4),W_m);

end % generate_A function

% --- ADDITIONAL FUNCTIONS ------------------------------------------------

function W_m = generate_W_m(code, fs, fc, Delta_tau, Delta_b, b)
% -------------------------------------------------------------------------
%   DESCRIPTION: This function generates the matrix W_m that is involved
%   in the A Matrix needed to compute the Misspecified Cramér-Rao bound for
%   the estimation of the time-daly and Doppler frequency of a signal in 
%   presence of a single multipath considering an additive white Gaussian
%   noise and under band-limited assumption.
% -------------------------------------------------------------------------
%% I - Initialisation

N = length(code);

wc = 2*pi*fc;
ts = 1/fs;

%% II - Priori calculations

% 1- Vectors involved in the W_Delta components

% time vector
dt = -(1-N:N-1) + (Delta_tau)/ts;
Nt = length(dt);

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

% Matrix U (due to different doppler between signals)
if (Delta_b ~= 0)
    diagU = exp(-1j*2*pi*ts*fc*(Delta_b)*(1:N));
else
    diagU = ones(1,N);
end

% 3- Intermediate combination of vectors and matrices to avoid huge
% matrices handling

WD_cD    = code'.*diagD;
WD_cDD   = WD_cD.*diagD;

if ( Delta_b ~= 0)
    WD_cU    = code'.*diagU;
    WD_cDU   = WD_cD.*diagU;
    WD_cDDU  = WD_cDD.*diagU;
else
    WD_cU    = code';
    WD_cDU   = WD_cD;    
    WD_cDDU  = WD_cDD;    
end

WD_VD0c = fftshift(ifft(fft(code,Nt).*(fft(v0)')));
WD_VD0c = WD_VD0c(1:N);

WD_VD1c = fftshift(ifft(fft(code,Nt).*(fft(v1)')));
WD_VD1c = WD_VD1c(1:N);

WD_VD2c = fftshift(ifft(fft(code,Nt).*(fft(v2)')));
WD_VD2c = WD_VD2c(1:N);

%% III - Evaluation of each components of W_Delta

% compute the conjugate of each terms
term_1 = (1/fs)*(WD_cU*WD_VD0c);
term_2 = (1/fs^2)*(WD_cDU*WD_VD0c);
term_3 = (1/fs^3)*(WD_cDDU*WD_VD0c);
term_4 = -(WD_cU*WD_VD1c) + (1j*wc*Delta_b/fs)*(WD_cU*WD_VD0c); 
term_5 = (-1/fs)*(WD_cU*WD_VD0c) - (1/fs)*(WD_cDU*WD_VD1c) + (1j*wc*Delta_b/fs^2)*(WD_cDU*WD_VD0c);
term_6 = (-fs)*(WD_cU*WD_VD2c) - (1j*2*wc*Delta_b)*(WD_cU*WD_VD1c) - (wc^2*Delta_b^2/fs)*(WD_cU*WD_VD0c);

%% IV - Computation of the matrix

W_m = exp(-1j*wc*b*Delta_tau)*[term_1', term_2', term_3', term_4', term_5', term_6'].';

end % generate_W_m function