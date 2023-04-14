function PHI = generate_Phi(code, fs, fc, epsilon)
%--------------------------------------------------------------------------
%        AUTHOR: Corentin Lubeigt
%       CREATED: 05/11/2020
%   
%   DESCRIPTION: This function generates the PHI matrix as it is defined
%   in [2]. This matrix is involved in the computation of the Fisher
%   Information Matrix for the estimation of the time-delay and the Doppler
%   frequency of two sources collected by a single antenna. In this 
%   approach, one considers considers the amplitude and phase of both 
%   sources as unknown nuisance parameters.
%
%        INPUTS: code       [] Nx1 vector containing the considered code
%                fs         [] reduced sampling frequency
%                fc         [] reduced carrier frequency
%                epsilon    [] 8x1 vector of parameters
%                              [tau_0, b_0, rho_0, phi_0, tau_1, b_1, rho_1, phi_1]^T
%
%       OUTPUTS: PHI        [] 8x8 matrix that corresponds to the PHI 
%                              matrix involved in the Fisher Information 
%                              of the vector of parameters epsilon in the 
%                              case of an additive white Gaussian circular 
%                              centered noise
%
%    REFERENCES: [1] 2020 [Lubeigt et al] Joint Delay-Doppler Estimation 
%                    Performance in a Dual Source Context
%                    DOI:10.3390/rs12233894
%                [2] 2022 [Lubeigt et al] Clean-to-Composite Bound Ratio: 
%                    A New Multipath Criterion for GNSS Signal Design and 
%                    Analysis
%                    DOI: 10.1109/TAES.2022.3172023
%--------------------------------------------------------------------------

%% I - Initialisation

NW = 3;
NQ = [2, 3];

wc = 2*pi*fc; % reduced pulsation


tau_0   = epsilon(1);
b_0     = epsilon(2);

tau_1   = epsilon(5);
b_1     = epsilon(6);

Delta_tau   = (tau_1 - tau_0)/fs; % reduced path separation
Delta_b     = b_1 - b_0;
    
%% II -  Evaluation of Matrices W, W_Delta and Q involved

% 1- Q matrices
Q_0 = [1j*wc*b_0,      0, -1;
               0, -1j*wc,  0];

Q_1 = [1j*wc*b_1,      0, -1;
               0, -1j*wc,  0];

Q = zeros(2*NQ);

Q(1:NQ(1),1:NQ(2)) = Q_0;
Q(NQ(1)+1:end,NQ(2)+1:end) = Q_1;


% 2- W matrices

% first term of Phi (DD1)
W = generate_W_Delta(code, fs, fc, 0, 0);

W_Delta = generate_W_Delta(code, fs, fc, Delta_tau, Delta_b);
W_Delta = W_Delta*exp(1j*wc*b_1*Delta_tau);

DD1 = zeros(2*NW);

DD1(1:NW,1:NW) = W;
DD1(NW+1:end,NW+1:end) = W;

DD1(1:NW,NW+1:end) = W_Delta';
DD1(NW+1:end,1:NW) = W_Delta;

% second term of Phi (DD2)
w1 = W(1,1);
W_Delta_11_exp = W_Delta(1,1);

w = W(:,1);
w_Delta_1dot = W_Delta(1,:);
w_Delta_dot1 = W_Delta(:,1);

DD2 = zeros(2*NW);

DD2(1:NW,1:NW) = w*w' + w_Delta_1dot'*w_Delta_1dot ...
                  - 2*real(w*w_Delta_1dot*(W_Delta_11_exp/w1)');
DD2(NW+1:end,NW+1:end) = w*w' + w_Delta_dot1*w_Delta_dot1' ...
                          - 2*real(w_Delta_dot1*w'*(W_Delta_11_exp/w1)');

W_temp = w*w_Delta_1dot + w_Delta_dot1*w' - w*w'*(W_Delta_11_exp/w1) ...
          - w_Delta_dot1*w_Delta_1dot*(W_Delta_11_exp/w1)'; 

DD2(1:NW,NW+1:end) = W_temp';
DD2(NW+1:end,1:NW) = W_temp;

DD2 = (w1/(w1^2 - abs(W_Delta_11_exp)^2)) * DD2;

%% III - Computation of the matrix PHI

PHI = Q*(DD1 - DD2)*Q';

% Note: here one still needs to take the real part and multiply by 
% 2*Fs/sigma^2 (see [2])

end % function generate_FIM