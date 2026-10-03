% Find the RMS value of the sinusoidal signal using both analytical (Equation 2.14) and digital
% (Equation 2.13) approaches.

% Solution, Digital
% Generate a 1-cycle sine wave. Assume a sample interval of Ts = 0.005 s and make N = 500 points
% long (both arbitrary). So solving Equation 1.4 for the total time: T = NTs = 0.005(500) = 2.5 s.
% To generate a single cycle given these parameters, the frequency of the sine wave should be
% f = 1/ T = 1/2.5 = 0.4 Hz. Set the amplitude of the sine wave, A = 1.0.

N = 500;               % Number of points for waveform
Ts = .005;             % Sample interval = 5 msec
t = (1:N)*Ts;          % Generate time vector (t = N Ts)
f = 1/(Ts*N);          % Sine wave freq. for 1 cycle
A = 1;                 % Sine wave amplitude
x = A*sin(2*pi*f*t);   % Generate sine wave
RMS = sqrt(mean(x.^2)) % Take the RMS value and output.


% Calculating the mean, variance, or standard deviation of a signal using MATLAB is straightforward.
% Assuming a signal vector x, these measurements are determined using one of the code
% lines below.
xm = mean(x); % Evaluate mean of x
xvar = var(x); % Variance of x normalizing by N-1
xstd = std(x); % Evaluate the standard deviation of x,