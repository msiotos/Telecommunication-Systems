%A
%a1
T=0.01;
a=[0 0.5 1];
over = 10;
Ts = T/over;
A = 4;

[phi1,t1] = srrc_pulse(T, Ts, A, a(1));
[phi2,t2] = srrc_pulse(T, Ts, A, a(2));
[phi3,t3] = srrc_pulse(T, Ts, A, a(3));

figure ();
plot(t1, phi1,'DisplayName','a = 0');
hold on; 
plot(t2, phi2,'DisplayName','a = 0.5');
plot(t3, phi3,'DisplayName','a = 1');
legend('show');
title('SRRC pulses for different roll-off');
xlabel('t');
ylabel('SRRC pulses'); 

%a2
Fs = 1/Ts;
N = 2048;
F= -Fs/2:Fs/N:Fs/2-Fs/N;


fourier1 = fftshift(fft(phi1,N)*Ts); 
fourier2 = fftshift(fft(phi2,N)*Ts);
fourier3 = fftshift(fft(phi3,N)*Ts);

phasm1 = abs(fourier1).^2;
phasm2 = abs(fourier2).^2;
phasm3 = abs(fourier3).^2;

figure();
plot(F, phasm1,'DisplayName','a = 0');
hold on;
plot(F, phasm2,'DisplayName','a = 0.5');
plot(F, phasm3,'DisplayName','a = 1');
legend('show');

title('Energy spectral density plot');
xlabel('Frequency (Hz)');
ylabel('Energy spectrum'); 
hold off;

figure(3);
semilogy(F, phasm1,'DisplayName','a = 0');
hold on;
semilogy(F, phasm2,'DisplayName','a = 0.5');
semilogy(F, phasm3,'DisplayName','a = 1');
legend('show');

title('Energy spectral density semilogy');
xlabel('Frequency (Hz)');
ylabel('Energy spectrum'); 


%a3
disp('Theory')
BW1=(1+0)/(2*T)
BW2=(1+0.5)/(2*T)
BW3=(1+1)/(2*T)

c1=T/10^3;
for ki= 1:length(F)-1
    c1=[c1 T/10^3];
end
c2=T/10^5;
for ki= 1:length(F)-1
    c2=[c2 T/10^5];
end
semilogy(F,c1);
semilogy(F,c2);
legend('a = 0', 'a = 0.5', 'a = 1','c1','c2');

%B

for k = 0:1:2
    figure;
    plot(t1, phi1, 'DisplayName', '\phi(t)');
    hold on;
    plot(t1 + k*T, phi1, 'DisplayName', '\phi(t - kT)');
    title(['a = 0 and k = ' num2str(k)]);
    hold off;
    legend('show');
    
    figure;
    plot(t2, phi2, 'DisplayName', '\phi(t)');
    hold on;
    plot(t2 + k*T, phi2, 'DisplayName', '\phi(t - kT)');
    title(['a = 0.5 and k = ' num2str(k)]);
    hold off;
    legend('show');
    
    figure;
    plot(t3, phi3, 'DisplayName', '\phi(t)');
    hold on;
    plot(t3 + k*T, phi3, 'DisplayName', '\phi(t - kT)');
    title(['a = 1 and k = ' num2str(k)]);
    hold off;
    legend('show');
end

for k = 0:3
    figure;
    offset1 = [zeros(1, length(0:Ts:k*T)) phi1(1:end - length(0:Ts:k*T))];
    result1 = phi1 .* offset1;
    plot(t1, result1, 'DisplayName', 'a = 0');
    integrated1(k + 1) = sum(result1) * Ts;
    
    %figure;
    offset2 = [zeros(1, length(0:Ts:k*T)) phi2(1:end - length(0:Ts:k*T))];
    result2 = phi2 .* offset2;
    hold on;
    plot(t2, result2, 'DisplayName', 'a = 0.5');
    integrated2(k + 1) = sum(result2) * Ts;
    
    %figure;
    offset3 = [zeros(1, length(0:Ts:k*T)) phi3(1:end - length(0:Ts:k*T))];
    result3 = phi3 .* offset3;
    plot(t3, result3, 'DisplayName', 'a = 1');
    title(['Product of \phi(t) and \phi(t - kT) for k = ' num2str(k)]);
    legend('show');
    hold off;
    integrated3(k + 1) = sum(result3) * Ts;
end

disp('Integration for a = 0: '); disp(integrated1)
disp('Integration for a = 0.5: '); disp(integrated2)
disp('Integration for a = 1: '); disp(integrated3)
