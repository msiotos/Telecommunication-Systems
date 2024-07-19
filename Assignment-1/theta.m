clear all;  
close all;  
% Parameters
T = 12;  % Let T=12
Ts = 0.01;  
t = -25:Ts:25; 
phi = zeros(1, length(t));
phi(abs(t) <= T/2) = 1/sqrt(T); 

% Reverse the time vector and phi for convolution
t_rev = -t(end:-1:1);
phi_rev = phi(end:-1:1);          
tconv = t(1) + t_rev(1) : Ts : t(end) + t_rev(end);

% Convolution to find the autocorrelation function
Rff = conv(phi, phi_rev) * Ts;
figure;
plot(tconv, Rff);
xlabel('\tau');
ylabel('Rff');
title('Autocorrelation function of \phi for \Theta 1');

figure();  
Rff = zeros(1, length(tconv));  % Initialize Rff
Rff(-T <= tconv & tconv < -T/2) = -1 - tconv(-T <= tconv & tconv < -T/2) / T;
Rff(-T/2 <= tconv & tconv < 0) = 1 + 3 * tconv(-T/2 <= tconv & tconv < 0) / T;
Rff(0 <= tconv & tconv < T/2) = 1 - 3 * tconv(0 <= tconv & tconv < T/2) / T;
Rff(T/2 <= tconv & tconv < T) = -1 + tconv(T/2 <= tconv & tconv < T) / T;
plot(tconv, Rff);
xlabel('\tau');
ylabel('Rff');
title('Autocorrelation function of \phi for \Theta 3');
