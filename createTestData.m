function createTestData()
%% BPSK调制
Fs=5e6; %采样率
N=10000000; %样本数
CA_op=7;
t=(0:N-1)/Fs; 
fi=1e3;  %载波频率（多普勒频移）
thetai=0;%初始相位
Noise=0.3;%噪声功率
load('CAlist.mat');
CA=CAs(CA_op,:); %CA码
CAWidth=Fs/(1.023e6+fi);
s=double(CA(t2n(t,Fs,0,CAWidth)));

cdata=s.*exp(1j*(2*pi*fi*t+thetai))+randn(1,N)*Noise;
save('testData.mat','cdata');
end

function n=t2n(t,Fs,Phi,CAWidth)
    n=mod(ceil((t*Fs+1+Phi)/CAWidth),1023);
    n(n==0)=1023;
end