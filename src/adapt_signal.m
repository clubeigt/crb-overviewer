function [dataSamp, timeSamp, N_fft] = adapt_signal(input_signal, OSF, time_axis)
    
N_signal = length(input_signal);
signal_padding = [zeros(N_signal,1);input_signal;zeros(N_signal,1)]; % Zero Padding

pw2	= floor(log2(length(signal_padding)));
if (log2(length(signal_padding)) ~= pw2)
    signal_padding = [signal_padding;zeros((2^(pw2+1))-numel(signal_padding),1)];
end

N_signal = length(signal_padding);
% Oversampling factor of frequency representation of input signals
% to compute the spline samples basis (power of 2)
fftOverSampFac	  = 4;
% Oversampling factor (in comparison with signal bandwidth) to compute the spline samples basis (should be a power of 2)
splineOverSampFac = 8;

N_fft        = N_signal*fftOverSampFac;
N_ifft       = splineOverSampFac*N_fft;
N_spline     = splineOverSampFac*N_signal;
realClass    = {'double'};
isRealDataIn = sum(abs(imag(signal_padding))) == 0; % is there any imaginary part in this signal

% Memory preallocation
iFreq       = zeros(N_fft,1);
dataZp      = complex(zeros(N_fft,1,realClass{1}));         
splineBasis = complex(zeros(N_ifft,1,realClass{1}));		
    
iFreq = [1:fix(N_fft/2)+1 , N_ifft - fliplr(0:fix(N_fft/2)-2)].';
% Fourier Transform
dataZp              = fft(signal_padding*splineOverSampFac,N_fft);
splineBasis(iFreq)  = dataZp;

splineBasis	= ifft(splineBasis);
splineBasis	= splineBasis(1:N_spline);

iDataOut = (1:N_signal*OSF).';
iSampOut = (iDataOut-1)*((1/OSF)*splineOverSampFac); %  ->  interp1(0:splineSize-1,...)

if (isRealDataIn)
    dataSamp = interp1(0:N_spline-1,real(splineBasis),iSampOut,'spline');
else
    dataSamp = interp1(0:N_spline-1,splineBasis,iSampOut,'spline');
end

%% zero padding removal
dataSamp = dataSamp(OSF*length(input_signal):OSF*(length(input_signal) + length(input_signal)) - 1);

N_fft = 2^(nextpow2(length(dataSamp)));

%% time vector
if(nargin==3)
    timeSamp = linspace(0,max(time_axis), length(dataSamp))';
else
    timeSamp = 0;
end

end