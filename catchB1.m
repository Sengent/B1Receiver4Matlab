function [fd,Cd,rate,dd]=catchB1(g,cdata,d1)
Fs=5e6;%������
L=length(cdata);
t=(0:L-1)/Fs;
%g=ca(a,b);

T=floor(t*2.046e6)+1;
cam=g(T);
CA=fft(double(cam));

% +-10kHz
%F=round(1e4*L/Fs);
F=20;
d=zeros(F+1,L);

%
for m=0:F
    %fd=(m-floor(F/2)+1)*1e4/F;
    data=cdata.*exp(-1j*2*pi*250*m*t);
    DATA=fft(data);
    
    %if(m>0)
    %    DATA=DATA([m+1:end,1:m]);
    %elseif(m<0)
    %    DATA=DATA([end+m+1:end,1:end+m]);
    %end
    %DATA=DATA;
    temp=ifft(conj(CA).*DATA);
    d(m+F+1,:)=abs(temp).^2;
    if(m==0)continue;end
    data=cdata.*exp(1j*2*pi*500*m*t);
    DATA=fft(data);
    %DATA=[DATA(1),conj(DATA(end:-1:2502)),conj(DATA(2501:-1:2))];
    %if(m>0)
    %    DATA=DATA([m+1:end,1:m]);
    %elseif(m<0)
    %    DATA=DATA([end+m+1:end,1:end+m]);
    %end
    %DATA=DATA;
    temp=ifft(conj(CA).*DATA);
    d(F+1-m,:)=abs(temp).^2;
end

dd=d+d1;
figure(100);mesh(dd);
d=dd(:,1:5000);
PEAK1=max(max(d));

[x,y]=find(d==PEAK1);
%Ĩ��first peak
for m=y-5:y+5
    for n=x-2:x+2
        mm=m;
        nn=n;
        if(m>L)
            mm=m-L;
        elseif(m<1)
            mm=m+L;
        end
        if(n>2*F+1||n<1)
            continue;
        end
        d(nn,mm)=0;
    end
end

PEAK2=max(max(d));
rate=PEAK1/PEAK2;
Cd=mod(y,5000);
fd=(x-F-1)*250;
%data=cdata.*exp(-1j*2*pi*fd*t);
%DATA=fft(data);
%temp=ifft(conj(CA).*DATA);
%plot(abs(temp).^2);
