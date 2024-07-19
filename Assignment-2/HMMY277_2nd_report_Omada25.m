%Askisi A
%A1
clear all; 
close all;
clc;
    
T = 10 ^ (-3);
over = 10; 
Ts = T/over;
A = 4;
a = 0.5;
Nf=2048;


[phi, t1] = srrc_pulse (T, Ts, A, a);   

Fs=1/Ts;

F = (-Fs/2):(Fs/Nf):(Fs/2)-(Fs/Nf);

Ft = fftshift(fft(phi,Nf))*Ts;

ESD=abs(Ft).^2;

figure();
semilogy(F,ESD);
grid on;
title('Energy Spectral Density');
xlabel('F(Hz)');
ylabel('|\Phi(F)|^2');

%A2
N=100;

b=(sign(randn(N,1))+1)/2;
X_n = bits_to_2PAM(b);
X_d=1/Ts*upsample(X_n,over);
t_Xd = 0:Ts:N*T-Ts;
X_t =conv(X_d,phi)*Ts;
t_conv = min(t1)+min(t_Xd):Ts:max(t1)+max(t_Xd);


figure()
plot(t_conv,X_t)
title('X(t)=\Sigma(Xn*\phi(t-n*T))');
xlabel('t(s)');
ylabel('X(t)');
grid on

%A3
T_total=length(t_conv)*Ts;
XF=fftshift(fft(X_t,Nf)*Ts);  
P_x=(abs(XF).^2)/T_total;

%plot 
figure();
plot(F,P_x)
title('Plot of one periodogram with 2-PAM configuration');
xlabel('F');
ylabel('Px(F)')
grid on;

%semilogy
figure();
semilogy(F,P_x);
title('Semilogy of periodogram with 2-PAM configuration');
xlabel('F');
ylabel('Px(F)')
grid on;

K = 500;
sum=0;

for k=0:K-1
    b=(sign(randn(N,1))+1)/2;
    X_n=bits_to_2PAM(b);
    X_d=1/Ts*upsample(X_n,over);
    X_t=conv(phi,X_d)*Ts;
    XF=fftshift(fft(X_t,Nf)*Ts);
    P_x=abs(XF).^2/T_total;
    sum=sum+P_x;
end

Sx=sum/K;
Sx_th= (var(X_n)/T).*(ESD);


figure();
semilogy(F,Sx,F,Sx_th);
title('Estimated Sx(F) and theoritical Sx(F) semilogy');
xlabel('F(Hz)');
ylabel('Sx(F)');
legend('estimated', 'theoritical');

%A4
b=(sign(randn(N,1))+1)/2;
X_n_4 = bits_to_4PAM(b);
t_Xd_4=(0:Ts:length(X_n_4)*T-Ts);
X_d_4=1/Ts*upsample(X_n_4,over);
X_t_4 =conv(phi,X_d_4)*Ts;
t_conv_4= min(t1)+min(t_Xd_4):Ts:max(t1)+max(t_Xd_4);

figure()
plot(t_conv_4,X_t_4);
title('X(t)=\Sigma(Xn*\phi(t-n*T))');
xlabel('t(s)');
ylabel('X(t)');
grid on

T_total=length(t_conv_4)*Ts; 
XF_4=fftshift(fft(X_t_4,Nf))*Ts;  
P_x_4=(abs(XF_4).^2)/T_total;

sum_4=0;

for k=0:K-1
    
    b=(sign(randn(N,1))+1)/2;
    X_n_4 = bits_to_4PAM(b);
    X_d_4=1/Ts*upsample(X_n_4,over);
    X_t_4 = conv(phi,X_d_4)*Ts;
    XF_4=fftshift(fft(X_t_4,Nf))*Ts;
    
    P_x_4=(abs(XF_4).^2)/T_total;
    sum_4=sum_4+P_x_4;

end

Sx_4=sum_4/K;
Sx_th_4= (var(X_n_4)/T).*(ESD);


%plot 
figure();
plot(F,P_x_4)
title('Plot of one periodogram with 4-PAM configuration');
xlabel('F');
ylabel('Px(F)')
grid on;


%semilogy
figure();
semilogy(F, P_x_4);
title('Semilogy of periodogram with 4-PAM configuration');
xlabel('F');
ylabel('Px(F)')
grid on;
    
figure();
semilogy(F,Sx_4,F,Sx_th_4);
title('Estimated Sx(F) and theoritical Sx(F) semilogy');
xlabel('F(Hz)');
ylabel('Sx(F)');
legend('estimated', 'theoritical');
grid on;

figure()
semilogy(F,Sx,F,Sx_4);
title('P(X)\_2PAM and P(X)\_4PAM')
xlabel('F(Hz)');
ylabel('P(X)\_2PAM , P(X)\_4PAM')
legend('2PAM','4PAM')

figure()
plot(F,Sx,F,Sx_4);
title('P(X)\_2PAM , P(X)\_4PAM')
xlabel('F(Hz)');
ylabel('P(X)\_2PAM , P(X)\_4PAM')
legend('2PAM','4PAM')

%A5
T = 2*T;
over = 2*over;

[phi, t1] = srrc_pulse (T, Ts, A, a);   

Fs=1/Ts;
F = (-Fs/2):(Fs/Nf):(Fs/2)-(Fs/Nf);
Ft = fftshift(fft(phi,Nf))*Ts;
ESD=abs(Ft).^2;
N=100;

b=(sign(randn(N,1))+1)/2;
X_n = bits_to_2PAM(b);
X_d=1/Ts*upsample(X_n,over);
t_Xd = 0:Ts:N*T-Ts;
X_t =conv(X_d,phi)*Ts;
t_conv = min(t1)+min(t_Xd):Ts:max(t1)+max(t_Xd);
T_total=length(t_conv)*Ts;
XF=fftshift(fft(X_t,Nf)*Ts);  
P_x=(abs(XF).^2)/T_total;

%plot 
figure();
plot(F,P_x)
title('Plot of one periodogram with 2-PAM configuration');
xlabel('F');
ylabel('Px(F)')
grid on;

%semilogy
figure();
semilogy(F,P_x);
title('Semilogy of periodogram with 2-PAM configuration');
xlabel('F');
ylabel('Px(F)')
grid on;

K = 500;
sum=0;

for k=0:K-1
    b=(sign(randn(N,1))+1)/2;
    X_n=bits_to_2PAM(b);
    X_d=1/Ts*upsample(X_n,over);
    X_t=conv(phi,X_d)*Ts;
    XF=fftshift(fft(X_t,Nf)*Ts);
    P_x=abs(XF).^2/T_total;
    sum=sum+P_x;
end

Sx=sum/K;
Sx_th= (var(X_n)/T).*(ESD);


figure();
semilogy(F,Sx,F,Sx_th);
title('Estimated Sx(F) and theoritical Sx(F) semilogy');
xlabel('F(Hz)');
ylabel('Sx(F)');
legend('estimated', 'theoritical');
T = T/2;
over = over/2;

%B
F0 = 100;
t = linspace(-A,A,2*F0);

X = randn(1,5);         
phi_b = (2*pi)*rand(1,5);

figure
hold on;
for i = 1 : 5
    Y = X(i)*cos(2*pi*F0*t + phi_b(i));
    txt = ['Y',num2str(i)];
    plot(t,Y,'DisplayName',txt)
end
grid on;
hold off
xlabel("time")
ylabel("Stochastic processes")
legend show;
title("Shared plot")
