function [i_tau_est, Fd_est, rho_est, phi_est,sigma2_est] = single_source_MLE(data, signal_clean, data_doppler, Fd_axis, OSF)
%--------------------------------------------------------------------------
%        AUTHOR: Corentin Lubeigt 
%                Lorenzo Ortega
%       CREATED: 05/18/2020
%   
%   DESCRIPTION: 
%
%        INPUTS: data      []      Fourier transform of the noisy composed 
%                                  signal whose parameters we try to estimate
%                signal_clean []   clean version of the main signal (no delay,
%                                  no phase, no doppler) oversampled
%                data_doppler []
%                Fd_axis      [Hz]
%                OSF          []
%                
%       OUTPUTS: i_tau_est  []    estimation of the index of the time delay 
%                                 of the main signal
%                Fd_est     [Hz]  estimation of the doppler drift of the main
%                                 signal
%                rho_est    []    estimation of the amplitude of the main
%                                 signal
%                phi_est    [rad] estimation of the phase of the main signal
%                sigma2_est []    estimation of the phase of the main signal
%     
%
%    REFERENCES: 
%--------------------------------------------------------------------------

%% I - Initialisation
i_tau_range = floor(length(signal_clean)/2);
% if OSF==1
dataSamp = signal_clean;
% else
%     [dataSamp, ~] = adapt_signal(signal_clean,OSF);
% end

% Fourier Transforms of the data_clean set for each possible doppler in
% b_axis
FT_data = fft(data);
FT_data_doppler_clean = fft(dataSamp.*data_doppler);
FTc_data_doppler_clean = conj(FT_data_doppler_clean);

normCorr	= sqrt(sum(abs(dataSamp).^2));
   
%% II - Loop

FT_data_current = FT_data;

% a- Estimation of the main signal parameters
corrFunc_data_current = fftshift(ifft(FT_data_current.*FTc_data_doppler_clean,[],1),1)./(abs(normCorr^2));

% a1) delay/doppler (tau/Fd)
[max_main,i_tau_max] = max(abs(corrFunc_data_current).^2,[],2);
[~,i_tau_est]  = max(max_main);

i_Fd_est = i_tau_max(i_tau_est); 
Fd_est = Fd_axis(i_Fd_est); 

% delay fine
if (OSF ~=1)
    corrFuncSamp = adapt_signal(corrFunc_data_current(:,i_Fd_est),OSF);
    [~, i_tau_est] = max(abs(corrFuncSamp));  

    % a2) amplitude (alpha)
    rho_est = norm(corrFuncSamp(i_tau_est));
    % a3) phase (phi)
    phi_est = angle(corrFuncSamp(i_tau_est));
else
    % a2) amplitude (alpha)
    rho_est = norm(corrFunc_data_current(i_tau_est,i_Fd_est));
    % a3) phase (phi)
    phi_est = angle(corrFunc_data_current(i_tau_est,i_Fd_est));
end


%% III- Output

% a- remove the origin offset
i_tau_est = i_tau_est - i_tau_range*OSF - 1;

% b- estimate sigma2
% sigma2_est = var(data - alpha_est*exp(1j*phi_est)*circshift(dataSamp.*(data_doppler(:,i_Fd_est)),i_tau_est));
sigma2_est = var(data);

end % end of single_source_MLE function