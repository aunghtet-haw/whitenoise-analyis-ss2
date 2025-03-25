clear all;
% general plot settings
set(groot,'defaultLineLineWidth',2.0);         % make lines thicker 
set(groot,'defaulttextinterpreter','latex');   % allow LaTeX syntax for labels, title, ...

% band-limited white noise
N = 1000;                % number of samples 
T_s = 0.01;              % sampling period
f_s = 1/T_s;             % sampling frequency
t = 0*T_s:T_s:(N-1)*T_s; % time axis for plotting
sigma = 10;              % standard deviation of sample values
mu = 0;                  % constant offset 
rng('shuffle');          % start MATLAB random number generator
gwn = sigma*randn(1,N);  % Gaussian white noise sequence
x = mu + gwn;            % input signal: deterministic offset plus stochastic white noise

figure(1); clf; % open / clear plot window
plot(t,x,'b');  % plot input signal (time domain)
xlabel('time $n\cdot T_s$ in s'); ylabel('input samples $x[n]$');
title(strcat('Input signal: offset $\mu=\,$', num2str(mu), ...
    ' plus Gaussian white noise with $\sigma=\,$', num2str(sigma)));
subtitle(strcat('sampling period $T_s=\,$', num2str(T_s),['' ...
    ', $N=\,$'], num2str(N),' samples'));
grid on;

figure(2); clf;     % open / clear plot window
histogram(x,20,'FaceColor','b'); 
xlabel('amplitude value $x$'); ylabel('number of values');
mean_x = mean(x);   % average value (offset)
sigma_x = std(x,1); % standard deviation
title('Distribution of input values:')
subtitle(strcat('$E(X)=\,$', num2str(mean_x,3), ...
    ', $\sigma(X)=\,$', num2str(sigma_x,5)));

figure(3); clf;           % open / clear plot window
[r_xx,lags] = xcorr(x,x); % autocorrelation estimate
r_xx = r_xx/N;            % 'biased' normalization by length of signal
r_xx_t = 0*lags;
plot(lags*T_s,r_xx,'b');
xlabel('lag time $k\cdot T_s$ in s'); ylabel('$r_{xx}[k]$');
title('Input signal autocorrelation')
grid on; % xlim([-50,50]);

%Power Spectrum 3 variants
figure(4); clf; % open / clear plot window

r_xx = circshift([r_xx 0],N+1); % zero-pad for even length and shift for first element to represent zero lag
Nr = length(r_xx);              % length is now Nr=2*N
S_xx_a = fft(r_xx)/(Nr);        % DFT to yield two-sided power spectrum, should be real and non-negative
I = max(imag(S_xx_a));          % for check that imaginary part is very small (caused by numerical noise)
S_xx_a = real(S_xx_a);          % remove nuisance imaginary part                                  
S_xx_a_one_sided = S_xx_a(1:Nr/2+1); % elements corresponding to non-negative frequencies 0...f_s/2
S_xx_a_one_sided(2:Nr/2) = 2*S_xx_a_one_sided(2:Nr/2); % double to account for negative frequency phasors
f_a = (0:Nr/2)*f_s/Nr;          % one-sided frequency axis corresponding to the one-sided spectrum

subplot(3,1,1)
plot(f_a,S_xx_a_one_sided,'b');
xlabel('frequency $f$ in Hz'); ylabel('DFT$\big(r_{xx}[k]\big)$');
title('\textbf{Input signal power spectrum}','fontsize',12)
subtitle('a) from DFT of autocorrelation')
grid on;

% S_xx_b_one_sided should hold the one-sided power spectrum of x and
% f_b the corresponding one-sided frequency axis
% your code goes here ... look at ss2lab3task2b (prof approach of calculating)
% Determine the length of the signal
N = length(x); % Total number of samples in the signal
% Compute the FFT of the signal
S_xx_b = fft(x); % Perform the Fast Fourier Transform (FFT) of the signal `x`
% Compute the two-sided power spectrum
P2 = abs(S_xx_b / N).^2; % Normalize FFT by dividing by N and square the magnitude to get power
% Compute the one-sided power spectrum
P1 = P2(1:floor(N/2) + 1); % Extract the first half of the spectrum (one-sided spectrum)
% Adjust power values to account for the omitted half of the spectrum
P1(2:end-1) = 2 * P1(2:end-1); % Double all values except for DC 0Hz (index 1) and Nyquist (last index if N is even)
% Define the frequency axis for the one-sided spectrum
f_b = (0:floor(N/2)) * (f_s / N); % Compute frequency values corresponding to each FFT bin
% Assign the one-sided power spectrum to a variable
S_xx_b_one_sided = P1; % Store the one-sided power spectrum for plotting or comparison
subplot(3,1,2)
plot(f_b,S_xx_b_one_sided,'b');
xlabel('frequency $f$ in Hz'); ylabel('$\left|\frac{2\cdot X[n]}{\sqrt{2}\cdot N}\right|^2$');
title('b) from DFT of signal')
grid on;

% option 'power' in pwelch() yields the power spectrum 
[S_xx_c,f] = pwelch(x,rectwin(N),0,N,f_s,'power');
subplot(3,1,3)
plot(f,S_xx_c,'b-')
xlabel('frequency $f$ in Hz'); ylabel('PS[n]');
title("c) from MATLAB's pwelch")
grid on;

% Calculate the mean power of the signal x in the following ways:
% your code goes here ...
% Calculate the mean power of the signal x in the following ways:

% Mean power from time-domain samples
P1 = mean(x.^2); % Power is the mean square value of the signal in the time domain

% Mean power from the autocorrelation function
P2 = r_xx(1); % The autocorrelation function's value at lag 0 corresponds to mean power since we circshift

% Mean power from two-sided power spectrum (a)
P3 = sum(S_xx_a); % Integrate two-sided spectrum, (multiply by frequency resolution only if we use psdensity)

% Mean power from one-sided power spectrum (a)
P4 = sum(S_xx_a_one_sided); % Integrate one-sided spectrum (already adjusted for energy)

% Mean power from one-sided power spectrum (b)
P5 = sum(S_xx_b_one_sided); % Integrate the one-sided spectrum from FFT of signal
% Mean power from one-sided power spectrum (c) (using pwelch)
P6 = sum(S_xx_c); % Integrate pwelch spectrum using the frequency resolution

% Display results
fprintf('Mean power calculated in various ways:\n');
fprintf('P1 (time domain): %.5f\n', P1);
fprintf('P2 (autocorrelation): %.5f\n', P2);
fprintf('P3 (two-sided spectrum a): %.5f\n', P3);
fprintf('P4 (one-sided spectrum a): %.5f\n', P4);
fprintf('P5 (one-sided spectrum b): %.5f\n', P5);
fprintf('P6 (one-sided spectrum c, pwelch): %.5f\n', P6);

[psd_x_sr,f_sr] = pwelch(x,rectwin(N),0,N,f_s,'psd'); % sr: single rectangular window
% compare with power spectrum from c) divided by spectral resolution:
psd_x_c = S_xx_c/(f_s/N);
e_psd = max(abs(psd_x_sr'-psd_x_c)); % the difference is on the order of 1e-15 (numerical noise)

% averaged psd with reduced spectral resolution using default pwelch call
N_win = 2^8; % reduced DFT length
[psd_x,f] = pwelch(x,[],[],N_win,f_s,'psd'); % length of one-sided power density spectrum is N_win/2+1

% Calculate mean power by integrating over the psd (Riemann sum)
% your code goes here ...
% Integrate over the averaged PSD, scaled by frequency resolution
P7 = sum(psd_x) * (f(2) - f(1)); % should be very close to P1...P6, difference caused by residue window effects

% Display P7 result
fprintf('P7 (from averaged PSD): %.5f\n', P7);
 

% compare raw and averaged power density spectrum
figure(5); clf; % open / clear plot window
plot(f_sr,psd_x_sr,'c',f,psd_x,'b-')
xlabel('frequency $f$ in Hz'); ylabel('power density in W/Hz');
title("Input signal power density spectrum from MATLAB's \tt{pwelch()}")
legend('one rectangular window','default: Hamming window + averaging')
grid on;

% compare raw and averaged power density spectrum
figure(5); clf; % open / clear plot window
plot(f_sr,psd_x_sr,'c',f,psd_x,'b-')
xlabel('frequency $f$ in Hz'); ylabel('power density in W/Hz');
title("Input signal power density spectrum from MATLAB's \tt{pwelch()}")
legend('one rectangular window','default: Hamming window + averaging')
grid on;

figure(6); clf; % open / clear plot window
semilogx(f,10*log10(psd_x),'b'); 
xlabel('frequency $f$ in Hz'); ylabel('$10\cdot \log_{10}\big(S_{yy}(f)\big)$ in dB/Hz');
title('Output signal power density spectrum')
grid on;

%RLC
% definition of the system parameters
R  = 1;           % resistance
C  = 0.01;        % capacitance
L  = 0.100;       % inductance
tau = R*C;        % constant 1
wr = 1/sqrt(L*C); % constant 2

% declaration of the continuous-time system
A2 = [0 -wr^2; 1 -tau*wr^2];
b2 = [1; 0];
c2 = [0 wr^2];
d2 = 0;
sys = ss(A2,b2,c2,d2);

figure(7);clf; % open / clear plot window

% frequency response of RLC resonant circuit as known 
% from Signals and Systems 1 lecture
s = j*2*pi*f; H = 1./(1+s*tau+s.^2/wr^2);
subplot(2,1,1);
semilogx(f,20*log10(abs(H)), 'm-')
ylabel({'amplitude response';'$|H(f)|$ in dB'});
title('\textbf{Theoretical frequency response of $RLC$ resonant filter}','fontsize',12);
subtitle(strcat('$R=\,$', num2str(R), ', $L=\,$', num2str(L), ', $C=\,$', num2str(C)));
grid on;
subplot(2,1,2)
semilogx(f,rad2deg(angle(H)), 'm-')
xlabel('frequency $f$ in Hz'); ylabel({'phase response';['$\varphi\big(H(f)\big)$ in $^\circ$']})
grid on;

% alternatively: use MATLAB bode command
figure(8); clf;
bode(sys);

figure(9); clf; % open / clear plot window
sys_d = c2d(sys,T_s);
%impulse(sys, sys_d);
step(sys,sys_d);
legend('continuous','discrete')

y = lsim(sys_d,x,t); % calculate output of sys_d for given input
figure(10); clf;      % open / clear plot window
plot(t,y,'rx');
xlabel('time $n\cdot T_s$ in s'); ylabel('output samples $y[n]$');
title(strcat('Output signal: $RLC$ filtered with $R=\,$', num2str(R), ', $L=\,$', num2str(L), ', $C=\,$', num2str(C)));
subtitle(strcat('constants $\tau=RC=\,$', num2str(tau,3), ', $\omega_r=\frac{1}{\sqrt{LC}}=\,$', num2str(wr,3)));

figure(11); clf; % open / clear plot window
histogram(y,'FaceColor','r'); 
xlabel('amplitude value $y$'); ylabel('number of values');
mean_y = mean(y);   % average value (offset)
sigma_y = std(y,1); % standard deviation
sigma_y_theo=sigma_x*sqrt(T_s/(2*tau)); % theoretical standard deviation
title("Distribution of output values")
subtitle(strcat('$E(Y)=\,$', num2str(mean_y,3), ...
    ' (theory:\,',num2str(mean_x,3),'), $\sigma(Y)=\,$', num2str(sigma_y,3), ...
    ' (theory:\,',num2str(sigma_y_theo,3), ')'));
grid on;

figure(12); clf;        % open / clear plot window
[r_yy,lags] = xcorr(y); % autocorrelation estimate
r_yy = r_yy/N;          % 'biased' normalization by length of signal
plot(lags*T_s,r_yy,'r');
xlabel('lag time $k\cdot T_s$ in s'); ylabel('$r_{yy}[k]$');
title('Output signal autocorrelation')
grid on; % xlim([-50,50]);

figure(13); clf; % open / clear plot window
[psd_y,f] = pwelch(y,[],[],N_win,f_s,'psd');
semilogx(f,10*log10(psd_y),'r'); 
xlabel('frequency $f$ in Hz'); ylabel('$10\cdot \log_{10}\big(S_{yy}(f)\big)$ in dB/Hz units');
title('Output signal power density spectrum')
grid on;

% Calculate the amplitude response estimate from power density spectra
amplitude_response_estimated = sqrt(psd_y ./ psd_x);

% Plot the amplitude response
figure(14); clf; % open / clear plot window
semilogx(f, 20 * log10(amplitude_response_estimated), 'k', 'DisplayName', 'experimental'); % black for experimental
hold on;
semilogx(f, 20 * log10(abs(H)), 'm-', 'DisplayName', 'theoretical'); % magenta for theoretical
xlabel('frequency $f$ in Hz', 'Interpreter', 'latex');
ylabel('$20 \cdot \log_{10} \left(\sqrt{\frac{S_{yy}(f)}{S_{xx}(f)}} \right)$', 'Interpreter', 'latex');
title('Amplitude response from auto power density spectrum', 'Interpreter', 'latex');
legend('show', 'Location', 'best');
grid on;

% Compute cross power spectral density using cpsd()
[pxy, f] = cpsd(y, x, [], [], [], f_s); % Cross-power spectral density between input and output
[psd_x, f] = pwelch(x, [], [], N_win, f_s); % Power spectral density of input
[psd_y, f] = pwelch(y, [], [], N_win, f_s); % Power spectral density of output

% Calculate amplitude response (magnitude of transfer function)
amplitude_response = abs(pxy) ./ psd_x; % H(f) magnitude

% Calculate phase response (angle of transfer function)
phase_response = unwrap(angle(pxy)); % H(f) phase in radians

% Plot the amplitude response
figure(15); clf; % open/clear plot window
subplot(2, 1, 1);
semilogx(f, 20 * log10(amplitude_response), 'k', 'DisplayName', 'Experimental');
hold on;
semilogx(f, 20 * log10(abs(H)), 'm-', 'DisplayName', 'Theoretical');
xlabel('Frequency $f$ [Hz]', 'Interpreter', 'latex');
ylabel('Amplitude Response $|H(f)|$ [dB]', 'Interpreter', 'latex');
title('Amplitude Response (Experimental vs. Theoretical)', 'Interpreter', 'latex');
legend('show', 'Location', 'best');
grid on;

% Plot the phase response
subplot(2, 1, 2);
semilogx(f, rad2deg(phase_response), 'k', 'DisplayName', 'Experimental');
hold on;
semilogx(f, rad2deg(angle(H)), 'm-', 'DisplayName', 'Theoretical');
xlabel('Frequency $f$ [Hz]', 'Interpreter', 'latex');
ylabel('Phase Response $\phi(H(f))$ [degrees]', 'Interpreter', 'latex');
title('Phase Response (Experimental vs. Theoretical)', 'Interpreter', 'latex');
legend('show', 'Location', 'best');
grid on;
