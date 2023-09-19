function [f] = Nchannel(N,fcut,fs,s,T)%fcut为低通滤波器的截止频率,s为speech signal,fs为采样频率
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
sumsig=[zeros(1,length(s))];%这步很重要
fmin=200;
fmax=7000;
dmin=log10(fmin/165.4+1)/0.06;
dmax=log10(fmax/165.4+1)/0.06;
d=(dmax-dmin)/N;
for i=1:N+1
    D(i)=d*(i-1)+dmin;
    f(i)=165.4*(10^(0.06*D(i))-1);%f(i)的第一个值始终为200Hz,f(i)的最后一个值始终为7000Hz
end
for i=1:N
   [b,a]=butter(4,[f(i) f(i+1)]/(fs/2));
    y=filter(b,a,s);
    Y=abs(y);
    [c,d]=butter(4,fcut/(fs/2));
    env=filter(c,d,Y);%通过低通滤波器得到的波
    f1=(f(i)+f(i+1))/2;
    n=1:length(s);%正弦波的长度一定要注意
    dt=n*(1/fs);
    sinesig=sin(2*pi*f1*dt);%产生相对应的正弦波信号
    x=env.*sinesig;%矩阵相乘一定要注意
    sumsig=sumsig+x;
end

sumsig=sumsig/norm(sumsig)*norm(s);%能量归一化
sound(sumsig,fs);

sig = repmat(sumsig,1,10);
[Pxx,w] = periodogram(sig,[],512,fs);

fftsum=fft(sumsig);%用fft近似求sum的CTFT
fftsum=fftshift(fftsum);%把零频率移至最中间
range=[-length(fftsum)/2:length(fftsum)/2-1];

figure;subplot(4,1,1),plot(sumsig);
title({['synthesized signal','  N=',num2str(N),'  fcut=',num2str(fcut)]});%num2str将数字转换为字符，并解决格式上的问题
xlabel('t');
subplot(4,1,2);plot(Pxx);
title('PSD');xlabel('frequency/Hz');
subplot(4,1,3);semilogy(w, Pxx);
title('PSD/dB');xlabel('frequency/Hz');
subplot(4,1,4),plot(range,abs(fftsum));
title({['FFT of synthesizeed signal','  N=',num2str(N),'   fcut=',num2str(fcut)]});%num2str将数字转换为字符，并解决格式上的问题
xlabel('frequency');

filename= ['project1 Task',num2str(T),' N=',num2str(N),'.wav'];
audiowrite(filename,sumsig,fs);
end