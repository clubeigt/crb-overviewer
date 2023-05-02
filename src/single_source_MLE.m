function [i_tau_est, Fd_est, rho_est, phi_est, sigma2_est] = single_source_MLE(data, replica, data_doppler, Fd_axis, OSF)
%--------------------------------------------------------------------------
%        AUTHOR: Corentin Lubeigt 
%                Lorenzo Ortega
%       CREATED: 05/18/2020
%   
%   DESCRIPTION: This function implements the maximum likelihood estimator
%   (MLE) for the estimation of the time-delay, Doppler frequency, 
%   amplitude and phase of a signal embedded in an additive white Gaussian 
%   noise with unknown variance. In this single source scenario, the MLE 
%   consists of finding the set [i_tau_est,Fd_est] that maximizes the
%   correlation function between the received noisy signal and a local
%   clean replica.
%
%        INPUTS: data         []   Nx1 noisy signal whose parameters are to
%                                  be estimated
%                replica      []   Nx1 clean version of the main signal (no
%                                  delay, no Doppler)
%                data_doppler []   NxNDoppler matrix of complex exponential
%                                  functions at different frequencies. The
%                                  frequencies correspond to the ones
%                                  present in the vector Fd_axis
%                Fd_axis      [Hz] 1xNDoppler vector of tested Doppler
%                                  frequencies
%                OSF          []   oversampling factor
%                
%       OUTPUTS: i_tau_est  []    estimated index of the time-delay
%                Fd_est     [Hz]  esitmated Doppler frequency
%                rho_est    []    estimated amplitude
%                phi_est    [rad] estimated phase
%                sigma2_est []    estimated noise variance
%     
%
%    REFERENCES: [1] 1993 [Ottersten et al] Exact and Large Sample Maximum 
%                    Likelihood Techniques for Parameter Estimation and 
%                    Detection in Array Processing, chapter extracted from 
%                    Radar Array Processing,
%                    DOI: 10.1007/978-3-642-77347-1_4
%--------------------------------------------------------------------------

%% I - Initialisation

i_tau_offset = floor(length(replica)/2); % needed due to the use of fftshift
normCorr	= sqrt(sum(abs(replica).^2)); % norm of the correlation function

% Fourier Transforms of the replica for each possible Doppler in Fd_axis
FT_data = fft(data);
FT_doppler = fft(replica.*data_doppler);
FTc_doppler = conj(FT_doppler);
  
%% II - MLE 

% 1- normalized cross-correlation function
corrFunc = fftshift(ifft(FT_data.*FTc_doppler,[],1),1)./(abs(normCorr^2));

% 2- parameters estimation
% delay/doppler (tau/Fd)
[max_main,i_tau_max] = max(abs(corrFunc).^2,[],2);
[~,i_tau_est]  = max(max_main);

i_Fd_est = i_tau_max(i_tau_est); 
Fd_est = Fd_axis(i_Fd_est); 

% oversampling of the cross-correlation to obtain a fine estimated of the
% time-delay
if (OSF ~=1)
    corrFuncSamp = adapt_signal(corrFunc(:,i_Fd_est),OSF);
    [~, i_tau_est] = max(abs(corrFuncSamp));  

    % amplitude (rho)
    rho_est = norm(corrFuncSamp(i_tau_est));
    % phase (phi)
    phi_est = angle(corrFuncSamp(i_tau_est));
else
    % amplitude (rho)
    rho_est = norm(corrFunc(i_tau_est,i_Fd_est));
    % phase (phi)
    phi_est = angle(corrFunc(i_tau_est,i_Fd_est));
end

%% III- Output

% 1- remove the offset
i_tau_est = i_tau_est - i_tau_offset*OSF - 1;

% 2- estimate noise variance
sigma2_est = var(data - rho_est*exp(1j*phi_est)*circshift(replica.*(data_doppler(:,i_Fd_est)),i_tau_est));

end % end of single_source_MLE function