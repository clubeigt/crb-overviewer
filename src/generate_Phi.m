function PHI = generate_Phi(code, Fs, Fc, epsilon)
%--------------------------------------------------------------------------
%         USAGE: generate_Phi(code, Fs, Fc, epsilon)
%
%        AUTHOR: Corentin Lubeigt
%       CREATED: 05/11/2020
%   
%   DESCRIPTION: This function generates the Phi matrix as it is defined
%   in [3]. This matrix is involved in the computation of the Fisher
%   Information Matrix for the estimation of the time delay and the Doppler
%   of two sources that impige a single antenna. In this approach, one considers
%   considers the amplitude and phase of both sources as nuisances.
%
%        INPUTS: code       []   1xnSamp vector containing the code used
%                Fs         [Hz] sampling frequency (can be a reduced frequency)
%                Fc         [Hz] carrier frequency (can be a reduced frequency)
%                epsilon    []   1x8 vector of parameters 
%                                [tau_0, b_0, alpha_0, phi_0, tau_1, b_1, alpha_1, phi_1]^T
%
%       OUTPUTS: FIM        []   nEstxnEst matrix that corresponds to the Fisher 
%                                Information of the vector of parameters 
%                                epsilon in the case of a Gaussian circular 
%                                centered noise 
%
%    REFERENCES: [1] 2019 [Das et al] On the Accuracy Limit of Time-delay Estimation with a Band-limited Signal
%                [2] 2020 [Lubeigt et al] Joint Delay-Doppler Estimation Performance in a Dual Source Context
%                [3] 2022 [Lubeigt et al] Clean-to-Composite Bound Ratio: A New Multipath Criterion for GNSS Signal Design and Analysis
%--------------------------------------------------------------------------

%% I - Initialisation
wc = 2*pi*Fc;

tau_0   = epsilon(1);
b_0     = epsilon(2);

tau_1   = epsilon(5);
b_1     = epsilon(6);

Delta_tau   = (tau_1 - tau_0)/Fs;
Delta_b     = b_1 - b_0;
    
%% II -  Evaluation of Matrices W, W_Delta and Q involved

% 1) Q matrices
Q_0 = [1j*wc*b_0,      0, -1;
               0, -1j*wc,  0];

Q_1 = [1j*wc*b_1,      0, -1;
               0, -1j*wc,  0];

Q = zeros(2*size(Q_0));

Q(1:length(Q_0(:,1)),1:length(Q_0(1,:))) = Q_0;
Q(length(Q_0(:,1))+1:end,length(Q_0(1,:))+1:end) = Q_1;


% 2) W matrices

% 2a- first term of Phi (DD1)
W = generate_W_Delta(code, Fs, Fc, 0, 0);

W_Delta = generate_W_Delta(code, Fs, Fc, Delta_tau, Delta_b);
W_Delta = W_Delta*exp(1j*wc*b_1*Delta_tau);

DD1 = zeros(2*size(W));

DD1(1:length(W(:,1)),1:length(W(1,:))) = W;
DD1(length(W(:,1))+1:end,length(W(1,:))+1:end) = W;

DD1(1:length(W_Delta(:,1)),length(W_Delta(1,:))+1:end) = W_Delta';
DD1(length(W_Delta(:,1))+1:end,1:length(W_Delta(1,:))) = W_Delta;

% 2b- second term of Phi (DD2)
w1 = W(1,1);
W_Delta_11_exp = W_Delta(1,1);

w = W(:,1);
w_Delta_1dot = W_Delta(1,:);
w_Delta_dot1 = W_Delta(:,1);

DD2 = zeros(2*size(W));

DD2(1:length(W(:,1)),1:length(W(1,:)))         = w*w' + w_Delta_1dot'*w_Delta_1dot - 2*real(w*w_Delta_1dot*(W_Delta_11_exp/w1)');
DD2(length(W(:,1))+1:end,length(W(1,:))+1:end) = w*w' + w_Delta_dot1*w_Delta_dot1' - 2*real(w_Delta_dot1*w'*(W_Delta_11_exp/w1)');

W_temp = w*w_Delta_1dot + w_Delta_dot1*w' - w*w'*(W_Delta_11_exp/w1) - w_Delta_dot1*w_Delta_1dot*(W_Delta_11_exp/w1)'; 

DD2(1:length(W_Delta(:,1)),length(W_Delta(1,:))+1:end) = W_temp';
DD2(length(W_Delta(:,1))+1:end,1:length(W_Delta(1,:))) = W_temp;

DD2 = (w1/(w1^2 - abs(W_Delta_11_exp)^2)) * DD2;

%% III - Computation of the matrix Phi

PHI = Q*(DD1 - DD2)*Q';

% Note: here one still needs to take the real part and multiply by 2*Fs/sigma^2 (see [3])

end % function generate_FIM