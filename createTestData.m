function createTestData()
%% BPSK����
Fs=5e6; %������
N=10000000; %������
CA_op=7;
t=(0:N-1)/Fs; 
fi=1e3;  %�ز�Ƶ�ʣ�������Ƶ�ƣ�
thetai=0;%��ʼ��λ
Noise=0.3;%��������
load('CAlist.mat');
CA=CAs(CA_op,:); %CA��
CAWidth=Fs/(1.023e6+fi);
s=double(CA(t2n(t,Fs,0,CAWidth)));

cdata=s.*exp(1j*(2*pi*fi*t+thetai))+randn(1,N)*Noise;
save('testData.mat','cdata');
end

function n=t2n(t,Fs,Phi,CAWidth)
    n=mod(ceil((t*Fs+1+Phi)/CAWidth),1023);
    n(n==0)=1023;
end