function FIM = generate_FIM(code, fs, fc, epsilon)
%--------------------------------------------------------------------------
%        AUTHOR: Corentin Lubeigt
%       CREATED: 05/11/2020
%   
%   DESCRIPTION: This function generates the Fisher Information Matrix for 
%   unknown deterministic parameters considering the collection with one
%   antenna of a signal and a single reflection in a Gaussian case and 
%   under band-limited assumption. This FIM is for the estimation of the 
%   vector of parameters epsilon defined by
%   epsilon = [tau_0, b_0, rho_0, phi_0, tau_1, b_1, rho_1, phi_1]^T 
%   where tau is the time-delay, b is the doppler dilatation coefficient, 
%   rho the signal amplitude (real-valued) and phi its phase. signal 0 
%   could be a direct signal and signal 1 its reflected replica.
%
%        INPUTS: code       [] 1xnSamp vector containing the code used
%                fs         [] reduced sampling frequency
%                fc         [] reduced carrier frequency 
%                epsilon    [] 8x1 vector of parameters 
%                              
%       OUTPUTS: FIM        [] 8x8 matrix that corresponds to the FIM of
%                              the vector of parameters epsilon in the 
%                              case of a Gaussian circular centered noise
%
%    REFERENCES: [1] 2019 [Das et al] On the Accuracy Limit of Time-delay 
%                    Estimation with a Band-limited Signal
%                [2] 2020 [Lubeigt et al] Joint Delay-Doppler Estimation 
%                    Performance in a Dual Source Context
%--------------------------------------------------------------------------

%% I - Initialisation

wc = 2*pi*fc; % reduced pulsation

tau_0   = epsilon(1);
b_0     = epsilon(2);
alpha_0 = epsilon(3);
phi_0   = epsilon(4);

tau_1   = epsilon(5);
b_1     = epsilon(6);
alpha_1 = epsilon(7);
phi_1   = epsilon(8);

Delta_tau   = (tau_1 - tau_0)/fs; %[chips]
Delta_b     = b_1 - b_0;
Delta_phi   = phi_1 - phi_0;      %[rad]

%% II - Evaluation of Matrices W, W_Delta and Q

% 1- W matrices
W = generate_W_Delta(code, fs, fc, 0, 0);
W_Delta = generate_W_Delta(code, fs, fc, Delta_tau, Delta_b);

DD = zeros(2*size(W));

DD(1:length(W(:,1)),1:length(W(1,:))) = W;
DD(length(W(:,1))+1:end,length(W(1,:))+1:end) = W;

DD(1:length(W_Delta(:,1)),length(W_Delta(1,:))+1:end) = (W_Delta*exp(1j*(Delta_phi+wc*Delta_tau*b_1)))';
DD(length(W_Delta(:,1))+1:end,1:length(W_Delta(1,:))) = W_Delta*exp(1j*(Delta_phi+wc*Delta_tau*b_1));

% 2- Q matrices
Q_0 = [1j*alpha_0*wc*b_0,              0, -alpha_0;
                       0, -1j*alpha_0*wc,        0;
                       1,              0,        0;
              1j*alpha_0,              0,        0];

Q_1 = [1j*alpha_1*wc*b_1,              0, -alpha_1;
                       0, -1j*alpha_1*wc,        0;
                       1,              0,        0;
              1j*alpha_1,              0,        0];

Q = zeros(2*size(Q_0));

Q(1:length(Q_0(:,1)),1:length(Q_0(1,:))) = Q_0;
Q(length(Q_0(:,1))+1:end,length(Q_0(1,:))+1:end) = Q_1;

%% III- Computation of the FIM

FIM = Q*DD*Q';

% Note: here one still needs to take the real part and multiply by 2*fs/sigma^2 (see [2])

end % function generate_FIM