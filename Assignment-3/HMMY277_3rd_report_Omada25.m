%Omada 25

clear all;
close all;

%Data
%Askisi A
N = 200;
A = 1;
A_srrc = 4;
T = 0.01;
over = 10;
Ts = T/over;
a=0.5;

Fs=1/Ts;
Nf=2048;

%1

bit_seq = (sign(randn(4*N,1))+1)/2;

%2,3

b1 = bit_seq(1:2*N);
b2 =  bit_seq(2*N + 1:4*N);

X_I = bits_to_4_PAM(b1, A);
X_Q = bits_to_4_PAM(b2, A);

%4

f = (-Fs/2):(Fs/Nf):(Fs/2)-(Fs/Nf);   

[phi,t1]=srrc_pulse(T, Ts, A_srrc, a);

%create signals
Xi_n = 1/Ts*upsample(X_I,over);
Xi_t = conv(Xi_n,phi)*Ts;

Xq_n = 1/Ts*upsample(X_Q,over);
Xq_t = conv(Xq_n,phi)*Ts;

%time vector
t2 = 0:Ts:N*T-Ts;
ti_conv = linspace(t1(1)+t2(1), t1(end)+t2(end),length(Xi_t));
tq_conv = linspace(t1(1)+t2(1), t1(end)+t2(end),length(Xq_t));

figure();
plot(ti_conv,Xi_t);
grid on;
title('Output waveform Xi(t)');
xlabel('t(s)');
ylabel('Xi(t)');

figure();
plot(tq_conv,Xq_t);
grid on;
title('Output waveform Xq(t)');
xlabel('t(s)');
ylabel('Xq(t)');

Ti_total = length(ti_conv)*T;
Tq_total = length(tq_conv)*T;

%PXF (from prev. assignment)

PXFi = ((abs(fftshift(fft(Xi_t,Nf))).^2)*Ts)./Ti_total;
PXFq = ((abs(fftshift(fft(Xq_t,Nf))).^2)*Ts)./Tq_total;

figure();
plot(f, PXFi);
title('Plot of periodogram of Xi(t)');
xlabel('F');
ylabel('Px(F)')
grid on;

figure();
plot(f, PXFq);
title('Plot of periodogram of Xq(t)');
xlabel('F');
ylabel('Px(F)')
grid on;

%5

F0 = 200;

Xi_mod_t = 2*Xi_t.*cos(2*pi*F0*ti_conv)*Ts;
Xq_mod_t = -2*Xq_t.*sin(2*pi*F0*tq_conv)*Ts;

figure();
plot(ti_conv,Xi_mod_t);
grid on;
title('Xi(t) multiplied with 2cos(2piFot)');
xlabel('t(s)');
ylabel('Xi(t)');

figure();
plot(tq_conv,Xq_mod_t);
grid on;
title('Xq(t) multiplied with -2sin(2piFot)');
xlabel('t(s)');
ylabel('Xq(t)');

T_total_i=length(Xi_mod_t)*Ts; %
T_total_q=length(Xq_mod_t)*Ts; %

PXFimod = ((abs(fftshift(fft(Xi_mod_t,Nf))).^2)*Ts)./T_total_i;
PXFqmod = ((abs(fftshift(fft(Xq_mod_t,Nf))).^2)*Ts)./T_total_q;

figure();
plot(f, PXFimod);
title('Plot of periodogram of Xi(t)_m_o_d');
xlabel('F');
ylabel('Px(F)')
grid on;

figure();
plot(f, PXFqmod);
title('Plot of periodogram of Xq(t)_m_o_d');
xlabel('F');
ylabel('Px(F)')
grid on;

%6

X_t_mod = Xi_mod_t + Xq_mod_t;

figure();
%ti_conv=tq_conv
plot(ti_conv,X_t_mod);
grid on;
title('Xmod(t) waveform');
xlabel('t(s)');
ylabel('X(t)');


T_total=length(X_t_mod)*Ts; 
PXFmodTotal = ((abs(fftshift(fft(X_t_mod,Nf))).^2)*Ts)./T_total;


figure();
plot(f, PXFmodTotal);
title('Plot of periodogram of Xmod(t)');
xlabel('F');
ylabel('Px(F)')
grid on;

%8

SNR = 20;
var_w = 10*(A^2)/(Ts*10^(SNR/10));
noise = sqrt(var_w)*randn(1,length(X_t_mod));
X_mod_noise = X_t_mod + noise;

%9

Xi_demod = X_mod_noise.*cos(2*pi*F0*ti_conv)*Ts;
Xq_demod = X_mod_noise.*(-1*sin(2*pi*F0*tq_conv))*Ts;

figure();
plot(ti_conv,Xi_demod);
grid on;
title('Xi(t)_d_e_m_o_d_u_l_a_t_e_d multiplied with cos(2piFot)');
xlabel('t(s)');
ylabel('Xi(t)');

figure();
plot(tq_conv,Xq_demod);
grid on;
title('Xq(t)_d_e_m_o_d_u_l_a_t_e_d multiplied with -sin(2piFot)');
xlabel('t(s)');
ylabel('Xq(t)');

T_total_i=length(Xi_demod)*Ts; 
XF_i=fftshift(fft(Xi_demod,Nf))*Ts;  
PxF_i=(abs(XF_i).^2)/T_total_i 

T_total_q=length(Xq_demod)*Ts; 
XF_q=fftshift(fft(Xq_demod,Nf))*Ts; 
PxF_q=(abs(XF_q).^2)/T_total_q


figure();
plot(f, PxF_i);
title('Plot of periodogram of Xi(t)_d_e_m_o_d_u_l_a_t_e_d multiplied with cos(2piFot)');
xlabel('F');
ylabel('Px(F)')
grid on;

figure();
plot(f, PxF_q);
title('Plot of periodogram of Xq(t)_d_e_m_o_d_u_l_a_t_e_d multiplied with -sin(2piFot)');
xlabel('F');
ylabel('Px(F)')
grid on;

%10
Xi_demod = conv(Xi_demod,phi)*Ts;
Xq_demod = conv(Xq_demod,phi)*Ts;

t_conv2 = min(ti_conv)+min(t1):Ts:max(t1)+max(ti_conv)
f = -Fs/2:Fs/Nf:Fs/2-Fs/Nf; 

figure()
plot(t_conv2,Xi_demod)
title('Output waveform of filtered Xi(t)');
xlabel('t(s)');
ylabel('Xi(t)');

figure()
plot(t_conv2,Xq_demod);
grid on;
title('Output waveform of filtered Xq(t)');
xlabel('t(s)');
ylabel('Xq(t)');

T_total_i=length(Xi_demod)*Ts; 
XF_i=fftshift(fft(Xi_demod,Nf))*Ts;
PxF_i=(abs(XF_i).^2)/T_total_i

T_total_q=length(Xq_demod)*Ts; 
XF_q=fftshift(fft(Xq_demod,Nf))*Ts;
PxF_q=(abs(XF_q).^2)/T_total_q 


figure();
plot(f, PxF_i);
title('Periodogram of filtered Xi(t)');
xlabel('F');
ylabel('Px(F)')
grid on;

figure();
plot(f, PxF_q);
title('Periodogram of filtered Xq(t)');
xlabel('F');
ylabel('Px(F)');
grid on;

%11

timesteps_for_sample_i = (2*A*T/Ts)+1 : over : length(Xi_demod)-(2*A*T/Ts);
timesteps_for_sample_q = (2*A*T/Ts)+1 : over : length(Xq_demod)-(2*A*T/Ts);

W_cos_sampled = Xi_demod(timesteps_for_sample_i);
W_sin_sampled = Xq_demod(timesteps_for_sample_q);

for i=1:N
Samples(i,1)=W_cos_sampled(i);
Samples(i,2)=W_sin_sampled(i); 
end 

scatterplot(Samples)
%12

for i=1:N
    X_I_possible(i) = detect_4_PAM(W_cos_sampled(i),A);
    X_Q_possible(i) = detect_4_PAM(W_sin_sampled(i),A); 
end

%13

errors_13 = 0;
actual_symbols = [X_I ; X_Q];
possible_symbols = [X_I_possible ; X_Q_possible];

for i=1:N
    if((actual_symbols(1, i) ~= possible_symbols(1, i)) || (actual_symbols(2, i) ~= possible_symbols(2, i)))
    errors_13 = errors_13 + 1;
    end
end

fprintf('Total errors for Task 13: %d\n', errors_13);

%15

errors_15 = 0;
iBits=PAM_4_to_bits(X_I_possible,A); 
qBits=PAM_4_to_bits(X_Q_possible,A);
bit_est =[iBits qBits];

for i=1:length(bit_est)
   if (bit_est(i) ~= bit_seq(i))
     errors_15 = errors_15 +1;
   end
end

fprintf('Total errors for Task 15: %d\n', errors_15);