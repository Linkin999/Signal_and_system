# Signals and Systems -- Project 1

> 张旭东 12011923		张宇潇 12011918		张一弛 11813226		杨恒宇 12012709
>

## Introduction

#### Tone-vocoder

本项目和下面的代码提供了一个多频带包络线索，可以合成语音信号，它可以使用带通滤波器将音频样本分割成不同的部分，然后使用校正低通滤波器，最后使用频移重建信号。

首先，输入信号将通过不同类型的带通滤波器从原始信号中分离出不同的频率。其次，信号的每一部分都将进行全波整流和具有完全相同频率的低通滤波器，以获得包络信号。最后，将得到的信号移位并求和以构造输出信号。

#### Code

```matlab
% function of Tone-vocoder

function [f] = Nchannel(N,fcut,fs,s,T)%fcut为低通滤波器的截止频率,s为speech signal,fs为采样频率,T为tasks数（保存用）
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
ylabel('PSD');xlabel('frequency/Hz');
subplot(4,1,3);semilogy(w, Pxx);
ylabel('PSD/dB');xlabel('frequency/Hz');
subplot(4,1,4),plot(range,abs(fftsum));
title({['FFT of synthesizeed signal','  N=',num2str(N),'   fcut=',num2str(fcut)]});%num2str将数字转换为字符，并解决格式上的问题
xlabel('frequency');

% 保存图片
figurename = ['figure/T',num2str(T),'N',num2str(N),'F',num2str(fcut),'.png'];
saveas(gcf,figurename);
close;

% 保存音频
filename= ['https://raw.githubusercontent.com/2235161562/SS_wav/pre1/wav/project1 Task',num2str(T),' N=',num2str(N),' fcut=',num2str(fcut),'.wav'];
audiowrite(filename,sumsig,fs);
end
```

## Task 1

在本任务中，我们需要将低通滤波器的截止频率设置为50Hz。为音调声码器定义不同的N，让N等于1,2,4,6,8等等，并让信号通过音调声码器。然后将继续信号保存为音频波，并通过收听分析输出音频。

#### 在N=1时

![T1N1](figure/T1N1F50.png)

#### 在N=2时

![T1N2](figure/T1N2F50.png)

#### 在N=4时

![T1N4](figure/T1N4F50.png)

#### 在N=6时

![T1N6](figure/T1N6F50.png)

#### 在N=8时

![T1N8](figure/T1N8F50.png)

#### 在N=32时

![T1N32](figure/T1N32F50.png)

#### 综合分析

​	从输出音频文件中，我们发现当N等于1到4时，输出音频是无法区分任何信息，声音非常尖锐。如果我们把N提高到6和8时，我们可以区分音频文件中的单词，尽管声音听起来仍然有些不自然，会有一些机器的尖锐。如果我们继续将N提高，声音的效果会更加的好。

​	同时从FFT图中我们可以得出结论，光谱与我们分析预期结果一致，每个输出信号都有N个峰值，对应于分布在不同频率的N个频带数。

#### 拓展尝试

​	验证分析，如果将N提高到了20，得到的声音合成感明显，但是声音的区分度会更好一点。

（见project1 Task1 N=32 fcut=50.wav）

​	在此如果将N提高到了80，得到的声音能产生语调的分别，声音非常明显。

见project1 Task1 N=80 fcut=50.wav

​	但当N达到一个较高的数值的时候，合成的声音会出现失真。

​	当N达到110，合成的声音会再次失真，变成一声比较快速的嘟嘟声。

见project1 Task1 N=110 fcut=50.wav

​	在N=100时，该失真不会出现。

见project1 Task1 N=100 fcut=50.wav

​	当N达到200时，基本无法识别信号

见project1 Task1 N=200 fcut=50.wav

​	信号的重构（包络）过程似乎不是线性的，这可能导致了在n>100上信号的失真与错位。因此，本实验的可懂度重建音频N低于或等于100，而这可能的原因也许和频率映射到耳蜗位置的功能有关。

#### Code

```matlab
% task 1
% fcut = 50
N = [1 2 4 6 8 16 32 64 70 80 90 100 110 128 200];
for n=1:length(N)
Nchannel(N(n),50,fs,x,1);
end
```

## Task2

将频带数N设置为4。将LPF截止频率定义为20 Hz、50 Hz、100 Hz和400 Hz，并让信号通过音调声码器。将继续信号保存为音频波，并通过收听分析输出音频。

#### 在 f<sub>cut</sub> = 20Hz时

![T2N2F20](figure/T2N2F20.png)

#### 在 f<sub>cut</sub> = 50Hz时

![T2N2F50](figure/T2N2F50.png)

#### 在 f<sub>cut</sub> = 100Hz时

![T2N2F100](figure/T2N2F100.png)

#### 在 f<sub>cut</sub> = 400Hz时

![T2N2F400](figure/T2N2F400.png)

#### 综合分析

​	从图中，我们发现，随着低通频率的频率变得更高，输出信号在峰值之间包含更多的散射信号，这是由于通过低通滤波器的信号增加所致。从频谱上，我们还可以通过比较频谱的低频部分来验证这一发现。随着低通滤波器频率的增加，输出音频文件也变得更加自然。

#### 拓展尝试

观察得到，当N合适时，随着截止频率fcut增加时，信号回声减小，质量得到改善，但会逐渐趋于饱和。而当截止频率很低时，信号噪声较大，难以识别。但当N的值很小时，变化更加明显

对比以下音频：

project1 Task2 N=2 fcut=10.wav

project1 Task2 N=2 fcut=100.wav

project1 Task2 N=2 fcut=200.wav

project1 Task2 N=2 fcut=400.wav

project1 Task2 N=2 fcut=1000.wav

当截至频率增加时，我们听到的声音越来越低沉，也逐渐清晰，之后趋于稳定（300Hz左右）。当截止频率过大时，信号质量反而出现异样，但总体质量仍相对较好。

#### 与 Task1 的比较

​	从Task1和Task2中，我们发现有两种不同的方法可以覆盖尽可能多的频率，同时在原始信号通过音调声码器后重建原始信号：增加频带数和增加包络LPF的低通频率。第一种方法在采样时增加频带以实现更紧密的覆盖，而第二种方法增加LPF（包络提取）的通过频率的带宽以使用更少的频带实现类似的结果。

#### Code
```matlab
% task 2
N = [2 4 90];
fcut = [20 50 100 400];
for n=1:length(N)
    for f=1:length(fcut)
    Nchannel(N(n),fcut(f),fs,x,2);
    end
end
```

## Task 3

在这个任务中我们要生成信噪比为 -5db  的噪音信号，然后以50Hz 的低通滤波器提取包络。

由于信噪比很低，此时声音信号中大部分是噪音。我们使用的是提供的语音 `C_01_01.wav`，语音的内容是“寄挂号信到北京需要多少钱”，原信号的时域和频域图像如下：

![T3Origin](figure/T3Origin.png)

#### 在N = 2时：

![T3N2F50](figure/T3N2F50.png)

从图中我们可以看出此时生成的信号仍然大部分是噪音信号，实际听感也基本上是噪音，无法分辨语音内容，说明N过小。从频域来看，构成信号的频率成分较少，与实际语音信号差距较大。

#### 在N= 4时

![T3N4F50](figure/T3N4F50.png)

从图中我们可以看到，相较于 N= 2的情形有一定的改进，即频域的频率组成更多（但是还远不够），但是时域信号来看仍然与语音信号差距较大，噪音仍然很强。实际听感仍无法分辨语音内容。

#### 在N= 6时

![T3N6F50](figure/T3N6F50.png)

从图中我们可以看到，相较于 N= 4的情形有一定的改进，即频域的频率组成更多（但是还远不够），但是时域信号来看仍然与语音信号差距较大，噪音仍然很强。实际听感仍无法分辨语音内容。

#### 在 N= 8时

![T3N8F50](figure/T3N8F50.png)

仍然噪音很大无法分辨语音。

#### 在 N＝16时

![T3N16F50](figure/T3N16F50.png)

相比之下，有一定改善，能隐约听到信息的内容，但是噪音很大，对于初次听的人很难分辨。

#### 综合分析

随着 N 的增加，意味着合成的时候使用了更多的频率的正弦信号叠加，使用更多的正弦信号可以逐渐逼近原来的声音信号，即合成语音的效果越好，越容易分辨。

#### 与 Task1 的比较

Task3 与 Task1 的区别主要是 Task3  要求使用的输入语音信号信噪比为 -5db ，即噪音信号能量很高的信号，在这样的输入信号下，当 Ｎ 的取值和 Task1一样的时候合成的语音效果远远不如Task1。

实测当 N 的取值在100左右的时候可以比较清楚地分辨语音信号的内容（但是会有比较高频的底噪），效果不如 Task1 。 

#### Code
```matlab
% −5dB
N=length(x);
sig=repmat(x,1,10);
[Pxx1,w1]=periodogram(sig,[],512,fs);
e=fir2(3000,w1/(fs/2),sqrt(Pxx1/max(Pxx1))); noise=1-2*rand(1,length(e)+N);
SSN=filter(e,1,noise);
SSN=SSN((length(e)+1):end);
normratio=10^(-5/20);
SSN=SSN*norm(x)/norm(SSN)/normratio; signal=x+SSN;

% task 3
N = [2 4 6 8 16];
for n=1:length(N)
Nchannel(N(n),50,fs,signal,3);
end
```

## Task 4

在这一个任务中，我们同样需要的是信噪比为 -5db的信号，固定 Ｎ=6，改变低通滤波器的截止频率。

#### 在 f<sub>cut</sub> = 20Hz时

![T4N6F20](figure/T4N6F20.png)

可以从图中发现时域信号仍然呈现明显噪音的特点，实际听感完全无法分辨语音内容。

#### 在 f<sub>cut</sub> = 50Hz时

![T4N6F50](figure/T4N6F50.png)

可以从图中发现时域信号仍然呈现明显噪音的特点，实际听感完全无法分辨语音内容。

#### 在 f<sub>cut</sub> = 100Hz时

![T4N6F100](figure/T4N6F100.png)

也同样几乎都是噪音，能听到几乎无法分辨的一点语音信息（语音信号的某些低频成分被保留下来），但是无法听到能正常辨识的语音信号。

#### 在 f<sub>cut</sub> = 400Hz时

![T4N6F400](figure/T4N6F400.png)

也同样几乎都是噪音，能听到几乎无法分辨的一点语音信息（语音信号的某些低频成分被保留下来），但是无法听到能正常辨识的语音信号。

#### 综合分析

随着低通滤波器的截止频率提高，意味着在最后提取信号的包络的时候会有更多高频成分被提取，即提取了更多合成之前各个频率段信号的细节，但是问题在于信号的信噪比非常小，而且给定的 Ｎ=6 在这样的信噪比下是很小的（从Task3中可以看出），这使得合成信号的频率组成非常单一，最后的声音听起来主要是比较高频的底噪（受到最后乘的正弦信号的频率影响），提取包络也无法得到理想的语音信息了。

#### Code
```matlab
% task 4
fcut = [20 50 100 400];
for n=1:length(fcut)
Nchannel(6,fcut(n),fs,signal,4);
end
```

## 讨论和感悟

从这个项目中，我们了解了Tone-vocoder的原理以及发声单元的数量对再现音频质量的影响。我们还了解了很多关于LaTex、Adobe Audition、Matlab和其他帮助我们完成这个项目的软件。