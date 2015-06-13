function [fd,Cd,rate,dd]=acquisit(g,cdata,d1)
global Fs CB_F L;
%Fs=5e6;%采样率
%Fs=5e6;%采样率
LENGTH=length(cdata);
t=(0:LENGTH-1)/Fs;
%g=ca(a,b);

T=floor(t*CB_F)+1;
cam=g(T);
CA=fft(double(cam));

% +-5kHz
%F=round(1e4*L/Fs);
F=20;
d=zeros(F+1,LENGTH);
step=250;
%
for m=0:F
    %fd=(m-floor(F/2)+1)*1e4/F;
    data=cdata.*exp(-1j*2*pi*step*m*t);
    DATA=fft(data);
    
    temp=ifft(conj(CA).*DATA);
    d(m+F+1,:)=abs(temp).^2;
    if(m==0)continue;end
    data=cdata.*exp(1j*2*pi*step*m*t);
    DATA=fft(data);
    
    temp=ifft(conj(CA).*DATA);
    d(F+1-m,:)=abs(temp).^2;
end

dd=d+d1;
figure(100);mesh(dd);
d=dd(:,1:L);
PEAK1=max(max(d));

[x,y]=find(d==PEAK1);
%抹掉first peak
for m=y-5:y+5
    for n=x-2:x+2
        mm=m;
        nn=n;
        if(m>LENGTH)
            mm=m-LENGTH;
        elseif(m<1)
            mm=m+LENGTH;
        end
        if(n>2*F+1||n<1)
            continue;
        end
        d(nn,mm)=0;
    end
end

PEAK2=max(max(d));
rate=PEAK1/PEAK2;
Cd=mod(y,L);
fd=(x-F-1)*step;

