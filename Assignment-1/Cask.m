T=0.01;
over=10;
a=0.5;
A=4;
N=100;
Ts=T/over;

%c1
b = (sign(randn(N, 1)) + 1)/2; 
disp('This is a random N-bits sequence')


%c2
%a
%transforming the N bits in to 2-pam
X = bits_to_2PAM(b);

%b
X_delta=1/Ts*upsample(X,over);
t_delta=(0:Ts:N*T-Ts);
figure()
plot(t_delta,X_delta);
title('X_delta pulse')
xlabel('t')
ylabel('X_delta(t)')
 
%c
[phi_apok,t_apok]=srrc_pulse(T,Ts,A,a);
X_t=conv(phi_apok,X_delta)*Ts;
tconv=[t_delta(1)+t_apok(1):Ts:t_delta(end)+t_apok(end)];
figure()
plot(tconv,X_t);
title('convolution of \phi(t) and X_delta')
ylabel('convolution of \phi(t) and X_delta')
xlabel('t')

%d

Z=conv(phi_apok,X_t)*Ts;
tconv2=[tconv(1)+t_apok(1):Ts:tconv(end)+t_apok(end)];
figure()
hold on
plot(tconv2,Z);
stem([0:N-1]*T,X);
title('convolution of \phi(-t) and X_t')
ylabel('convolution of \phi(-t) and X_t')
xlabel('t')
legend('Z(t)','X')
 