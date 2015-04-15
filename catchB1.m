function [fd,Cd,rate,dd]=catchB1(g,cdata,d1)
Fs=5e6;%²ÉÑùÂÊ
L=length(cdata);
t=(1:L)/Fs;
%g=ca(a,b);

T=mod(ceil(t*2.046e6)-1,2046)+1;
cam=g(T);
CA=fft(double(cam));

% +-10kHz
%F=round(1e4*L/Fs);
F=40;
d=zeros(F+1,L);

%
for m=-F:F
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
end

dd=d+d1;
figure(100);mesh(dd);
d=dd;
PEAK1=max(max(d));

[x,y]=find(d==PEAK1);
%Ä¨µôfirst peak
for m=y-5:y+5
    for n=x-2:x+2
        mm=m;
        nn=n;
        if(m>L)
            mm=m-L;
        elseif(m<1)
            mm=m+L;
        end
        if(n>2*F+1)
            nn=n-2*F+1;
        elseif(n<1)
            nn=n+2*F+1;
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
