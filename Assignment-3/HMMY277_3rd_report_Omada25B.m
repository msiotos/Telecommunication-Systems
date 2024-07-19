%Omada 25
clear all;
close all;

%Askisi B
K = 200;
Esymbols_per_SNR = zeros(1, 9);
Ebits_per_SNR = zeros(1, 9);
s=1;
for SNR = 0:2:16
    for j=1:K
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
      ti_conv = t1(1)+t2(1): Ts : t1(end) + t2(end);
      tq_conv = t1(1)+t2(1): Ts : t1(end) + t2(end);
      Ti_total = length(ti_conv)*T;
      Tq_total = length(tq_conv)*T;
      %5
      F0 = 200;
      Xi_mod_t = 2*Xi_t.*cos(2*pi*F0*ti_conv)*Ts;
      Xq_mod_t = -2*Xq_t.*sin(2*pi*F0*tq_conv)*Ts;
      T_total_i=length(Xi_mod_t)*Ts; %
      T_total_q=length(Xq_mod_t)*Ts; %
      PXFimod = ((abs(fftshift(fft(Xi_mod_t,Nf))).^2)*Ts)./T_total_i;
      PXFqmod = ((abs(fftshift(fft(Xq_mod_t,Nf))).^2)*Ts)./T_total_q;
      %6
      X_t_mod = Xi_mod_t + Xq_mod_t;
      T_total=length(X_t_mod)*Ts;
      PXFmodTotal = ((abs(fftshift(fft(X_t_mod,Nf))).^2)*Ts)./T_total;
      %8
      %SNR = 20;
      var_w = (10*(A^2))/(Ts*10^(SNR/10));
      noise = sqrt(var_w)*randn(1,length(X_t_mod));
      X_mod_noise = X_t_mod + noise;
      %9
      W_cos = X_mod_noise.*cos(2*pi*F0*ti_conv)*Ts;
      W_sin = X_mod_noise.*(-sin(2*pi*F0*tq_conv))*Ts;
      ti_total=length(W_cos)*Ts;
      tq_total=length(W_sin)*Ts;
      PXFW_cos=(abs(fftshift(fft(W_cos,Nf))*Ts).^2)/ti_total;
      PXFW_sin=(abs(fftshift(fft(W_sin,Nf))*Ts).^2)/tq_total;
      %10
      filtered_W_cos = conv(W_cos,phi)*Ts;
      filtered_W_sin = conv(W_sin,phi)*Ts;
      %ti_conv=tq_conv
      t_conv2 = min(ti_conv)+min(t1) : Ts : max(t1)+max(ti_conv);
      ti_total2 = length(filtered_W_cos)*Ts;
      tq_total2 = length(filtered_W_sin)*Ts;
      PXFfiltered_W_cos=(abs(fftshift(fft(filtered_W_cos,Nf))*Ts).^2)/ti_total2;
      PXFfiltered_W_sin=(abs(fftshift(fft(filtered_W_sin,Nf))*Ts).^2)/tq_total2;
      %11
      timesteps_for_sample_i = (2*A*T/Ts)+1 : over : length(filtered_W_cos)-(2*A*T/Ts);
      timesteps_for_sample_q = (2*A*T/Ts)+1 : over : length(filtered_W_sin)-(2*A*T/Ts);
      W_cos_sampled = filtered_W_cos(timesteps_for_sample_i);
      W_sin_sampled = filtered_W_sin(timesteps_for_sample_q);
      for i=1:N
       Samples(i,1)=W_cos_sampled(i);
       Samples(i,2)=W_sin_sampled(i); 
      end 
      %scatterplot(Samples)
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
    end
    Esymbols_per_SNR(s) = errors_13;
    Ebits_per_SNR(s) = errors_15;
    s = s+1;
end

function y = Q(x)
    y = 0.5 * erfc(x / sqrt(2));
end

th = 1;
for SNR=0:2:16 
 
    Esymbols_per_SNR_theor(th) = 3*Q(sqrt(0.2.*(10.^(SNR/10)))); 
    Ebits_per_SNR_theor(th) = Esymbols_per_SNR_theor(th)./4; 
    th=th+1;     
end

%B2
SNR=[0:2:16]; 
figure(); 
semilogy(SNR, Esymbols_per_SNR); 
hold on; 
semilogy(SNR, Esymbols_per_SNR_theor); 
hold off; 
grid on;
legend ('Predicted','Theoretical') 
title('Symbols Error Rate');

%B3
figure(); 
semilogy(SNR, Ebits_per_SNR); 
hold on; 
semilogy(SNR, Ebits_per_SNR_theor); 
hold off; 
grid on;
legend ('Predicted','Theoretical') 
title('Bits Error Rate');

