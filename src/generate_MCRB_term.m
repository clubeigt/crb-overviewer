function MCRB_term = generate_MCRB_term(epsilon,code, Fs, Fc)


%% I - 
wc = 2*pi*Fc;

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

Delta_tau_0   = (tau_0 - tau_pt)/Fs; %[chips]
Delta_tau_1   = (tau_1 - tau_pt)/Fs; %[chips]

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

%% III - W^m matrices

W_m_0  = generate_W_m(code, Fs, Fc, Delta_tau_0, Delta_b_0, b_0);
W_m_1  = generate_W_m(code, Fs, Fc, Delta_tau_1, Delta_b_1, b_1);
W_m_pt = generate_W_m(code, Fs, Fc, 0, 0, b_pt);

W_m = conj(alpha_0)*W_m_0 + conj(alpha_1)*W_m_1 - conj(alpha_pt)*W_m_pt;

MCRB_term = [Q_tau, Q_b, Q_rho, Q_phi]*kron(eye(4),W_m);
end

function W_m = generate_W_m(code, Fs, Fc, Delta_tau, Delta_b, b)

%% I - Initialisation

nSamp = length(code);
wc = 2*pi*Fc;
Ts = 1/Fs;

%% II - Priori calculations

% 1) Vectors involved in the W_Delta components

% 1a- time vector
dt = -(1-nSamp:nSamp-1) + (Delta_tau)/Ts;

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
    diagU = exp(-1j*2*pi*Ts*Fc*(Delta_b)*(1:nSamp));
else
    diagU = ones(1,nSamp);
end

% 3) Intermediate combination of vectors and matrices to avoid huge
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

WD_VD0c = fftshift(ifft(fft(code,length(v0)).*(fft(v0)')));
WD_VD0c = WD_VD0c(1:length(code));

WD_VD1c = fftshift(ifft(fft(code,length(v1)).*(fft(v1)')));
WD_VD1c = WD_VD1c(1:length(code));

WD_VD2c = fftshift(ifft(fft(code,length(v2)).*(fft(v2)')));
WD_VD2c = WD_VD2c(1:length(code));

%% III - Evaluation of each components of W_Delta

% compute the conjugate of each terms
term_1 = (1/Fs)*(WD_cU*WD_VD0c);
term_2 = (1/Fs^2)*(WD_cDU*WD_VD0c);
term_3 = (1/Fs^3)*(WD_cDDU*WD_VD0c);
term_4 = -(WD_cU*WD_VD1c) + (1j*wc*Delta_b/Fs)*(WD_cU*WD_VD0c); 
term_5 = (-1/Fs)*(WD_cU*WD_VD0c) - (1/Fs)*(WD_cDU*WD_VD1c) + (1j*wc*Delta_b/Fs^2)*(WD_cDU*WD_VD0c);
term_6 = (-Fs)*(WD_cU*WD_VD2c) - (1j*2*wc*Delta_b)*(WD_cU*WD_VD1c) - (wc^2*Delta_b^2/Fs)*(WD_cU*WD_VD0c);

%% IV - Computation of the matrix

W_m = exp(-1j*wc*b*Delta_tau)*[term_1', term_2', term_3', term_4', term_5', term_6'].';

end