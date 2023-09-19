%signal and system project 2
clear
clc
close all

% 求取信道增益 channel gain
M = 4;
xp = randn(1,32); 
H = channel_est(xp,M);
figure
subplot(2,1,1)
stem(0:length(H)-1,real(H))
title("Channel gain-real component")
subplot(2,1,2)
stem(0:length(H)-1,imag(H))
title("Channel gain-imaginary component")

%重新生成随机信号进行仿真
data = randn(1,32);
figure
stem(0:length(data)-1,data);
title("Data to be sent")
datasig_cp = trans_end(data,M);
wireless_sig = transmitter(datasig_cp);
signal_picked_up = wireless(wireless_sig);
sig_demodulated = receiver(signal_picked_up);
recovered_sig = recv_end(sig_demodulated(1:length(datasig_cp)),M,H);
figure
subplot(3,1,1)
stem(0:length(recovered_sig)-1,recovered_sig)
hold on
stem(0:length(data)-1,data)
legend("Recovered data", "Data")
subplot(3,1,2)
stem(0:length(data)-1,(recovered_sig-data))
title("Absolute Error ")
subplot(3,1,3)
stem(0:length(data)-1,(recovered_sig-data)./data)
title("Relative Error ")
figure
stem(0:length(H)-1,ifft(H));
title("Verifying the calculated h[n]")
legend("inverse FFT of H")


% block 1
function sig_with_cp = trans_end(sig,M)
%to find the signal with cp
%here set M as the length of cyclic prefix
%sig is the input signal
    x = ifft(sig).*length(sig);
    figure
    subplot(2,1,1)
    stem(0:length(x)-1,real(x));
    title("Inverse fft of data---real part")
    
    subplot(2,1,2)
    stem(0:length(x)-1,imag(x));
    title("Inverse fft of data---odd part")
    
    figure
    sig_with_cp = [x(end-M+1:end),x];
    subplot(2,1,1)
    stem(0:length(sig_with_cp)-1,real(sig_with_cp));
    title("Inverse fft of data(with CP)---real part")
    subplot(2,1,2)
    stem(0:length(sig_with_cp)-1,imag(sig_with_cp));
    title("Inverse fft of data(with CP)---odd part")
end
 
% block2
function recovered_sig = recv_end(recsig,M,H)
    de_cp_sig = recsig(M+1:end);
    figure
    subplot(2,1,1)
    stem(0:length(de_cp_sig)-1,real(de_cp_sig))
    title("Received signal without CP---real part")
    subplot(2,1,2)
    stem(0:length(de_cp_sig)-1,imag(de_cp_sig))
    title("Received signal without CP---imaginary part")
    
    sig = fft(de_cp_sig)./length(de_cp_sig);
    recovered_sig = real(sig(1:length(H))./H);
    
    figure
    subplot(2,1,1)
    stem(0:length(recovered_sig)-1,real(recovered_sig))
    title("Recoverd signal ---real part")
    subplot(2,1,2)
    stem(0:length(recovered_sig)-1,imag(recovered_sig))
    title("Recovered signal---imaginary part")
end
 
% block3
function wireless_sig = transmitter(sig_with_cp)
    % T = 1 us;
    T = 0.000001;
    %carrier frequency 
    wc = 100000000;
    wc_rad= wc.*2.*pi;
    %because sampling frequency should be larger than 2*wc
    %we can set sampling frequency as 10*wc
    fs = wc.*10;
    x = sig_with_cp;
    N = length(x);
    %time step
    timestep= 1/fs;
    t0 = 0:timestep:(T-timestep);
    %number of points in one T
    np = length(t0);
    t = 0:timestep:(N.*T-timestep);
    
    xp = zeros(1,N*length(t0));
    for k = 0:1:N-1
        xp(k*np+1) = x(k+1);
    end
    %
    xc = zeros(1,N*np);
    for k = 0:1:N-1
        xc((k*np+1):(k*np+np)) = x(k+1);
    end
    % real and imaginary parts
    xr = real(xc);
    xi = imag(xc);
    
    %generate sinsig
    sinsig = sin(wc_rad.*t);
    cossig = cos(wc_rad.*t);
    %modulation
    xcos = xr.*cossig;
    xsin = xi.*sinsig;
    
    %the final output
    wireless_sig= xcos+xsin;
    
    %plot figures
    figure
    %xc
    subplot(2,1,1)
    plot(real(xc))
    hold on
    stem(real(xp))
    hold off
    legend("xc","xp")
    title("Real part of xp and xc")
    
    subplot(2,1,2)
    plot(imag(xc))
    hold on
    stem(imag(xp))
    hold off
    legend("xc","xp")
    title("Imaginary part xp and xc")
    
    figure
    plot(t,wireless_sig)
    xlabel("time/s")
    title("Modulated signal")
    
    figure
    fft_wireless_sig = fft(wireless_sig);
    subplot(2,1,1)
    plot(-length(fft_wireless_sig)./2:length(fft_wireless_sig)./2-1,fftshift(real(fft_wireless_sig)))
    title("Real part of fft of modulated signal")
    subplot(2,1,2)
    plot(-length(fft_wireless_sig)./2:length(fft_wireless_sig)./2-1,fftshift(imag(fft_wireless_sig)))
    title("Imaginary part of fft of modulated signal")
end
 
% block 4
function sig_demodulated = receiver(recv_sig)
% basic signal parameters
    wc=100000000;
    wc_rad = wc.*2.*pi;
    fs = 10*wc;
    sig = recv_sig;
    timestep = 1/fs;
 
%demodulation
np = length(sig);%number of points
t = 0:timestep:(np.*timestep-timestep);
 
sigcos = cos(wc_rad.*t).*sig;
sigsin = sin(wc_rad.*t).*sig;
 
sig_real=fft(sigcos);
sig_imag=fft(sigsin);
figure
subplot(2,1,1)
plot(-length(sigcos)./2:length(sigcos)./2-1,fftshift(abs(sig_real)))
title("fft of the received signal---real part")
subplot(2,1,2)
plot(-length(sigcos)./2:length(sigcos)./2-1,fftshift(abs(sig_imag)))
title("fft of the received signal---imaginary part")
 
figure
subplot(2,2,1)
plot(-length(sigcos)./2:length(sigcos)./2-1,fftshift(abs(sig_real)))
title("fft of the received signal---real part")
subplot(2,2,2)
plot(-length(sigcos)./2:length(sigcos)./2-1,fftshift(abs(sig_imag)))
title("fft of the received signal---imaginary part")
  
subplot(2,2,3)
sig_real(length(sigcos)./8:length(sigcos)*7./8)=0;
plot(-length(sigcos)./2:length(sigcos)./2-1,fftshift(abs(sig_real)))
title("fft of the filtered signal---real part")
subplot(2,2,4)
sig_imag(length(sigsin)./8:length(sigsin)*7./8)=0;
plot(-length(sigsin)./2:length(sigsin)./2-1,fftshift(abs(sig_imag)))
title("fft of the filtered signal---imaginary part")
 
sigcos = 2.*ifft(sig_real);
sigsin = 2.*ifft(sig_imag);
  
% final synthesis
sig_dem = sigcos+sigsin.*1j;
%ADC
%sig_dem = [sig_dem,0]
figure
subplot(2,1,1)
plot(0:length(sig_dem)-1,real(sig_dem))
xlabel("time/s")
title("Demodulated signal---real part")
subplot(2,1,2)
plot(0:length(sig_dem)-1,imag(sig_dem))
xlabel("time/s")
title("Demodulated signal---imaginary part")
N = length(sig_dem)/(1000);
yd = zeros(1,N);
for k = 1:N
yd(k) = mean(sig_dem(((k-1).*1000+100):(k.*1000-100)));  % here 1000 = fs.*T matlab gives warning because fs.*T is calulated as 1000.00(float type)
end %求取一个采样周期内的平均值 作为ADC

figure
sig_demodulated = yd; % 与原信号相比，存在很小的复数部分
subplot(2,1,1)
stem(0:length(sig_demodulated)-1,real(sig_demodulated))
title("Signal demodulated---real part")
subplot(2,1,2)
stem(0:length(sig_demodulated)-1,imag(sig_demodulated))
title("Signal demodulated---imaginary part")
 
end
 

function H = channel_est(commonsig,M)
%channel gain
    commonsig_cp = trans_end(commonsig,M);
    wireless_sig = transmitter(commonsig_cp);
    signal_picked_up = wireless(wireless_sig);
 
    sig_demodulated = receiver(signal_picked_up);
    H = ones(1,32);
    recsig = sig_demodulated(1:length(commonsig_cp));
    de_cp_sig = recsig(M+1:end);
    sig = fft(de_cp_sig)./length(de_cp_sig);
    recovered_sig = sig(1:length(H));
    % recovered_sig = recv_end(sig_demodulated(1:length(commonsig_cp)),M,H);
    figure
    stem(0:length(commonsig)-1,recovered_sig(1:length(commonsig)))
    hold on
    H = recovered_sig./commonsig;% ÔöÒæÏòÁ¿
    
end
 
 
function sig_picked_up = wireless(wireless_sig)
% 这个函数计发射信号和无线系统的h(t)进行卷积
fs = 1000000000;
T = 0.000001;
timestep = 1/fs;
t0 = 0:timestep:4.*T;
h = [0.5,zeros(1,1499),0.4,zeros(1,999),0.35,zeros(1,499),0.3,zeros(1,1000)];
sig_picked_up = conv(wireless_sig,h);
 
N = (length(wireless_sig)+length(h))./1000;
t = 0:timestep:(N).*T-2.*1/fs;

figure()
subplot(2,1,1)
plot(t0,h)
xlabel("time/s")
title("h(t) of actual wireless channel")
subplot(2,1,2)
plot(t,sig_picked_up)
xlabel("time/s")
title("In-air modulated signal")
end
