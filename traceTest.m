clear;
load('data.mat');
load('CBlist.mat');
CB=CBs(7,:);
Fs=5e6;
length=2000;


Phi=0;
for m=1:1000    
    [Fi,Phi,rate]=catchB1(CB,cdata(1:5000));
    if(rate>1.5)
        break;
    end
    cdata=cdata(5001:end);
end
cdata=cdata(Phi+1:end);
[Fi,Phi,rate]=catchB1(CB,cdata(1:5000));

CAWidth=Fs/(2.046e6);
Thetai=0;
CM=[0,0];
CTcoh=10;


P=zeros(1,length-1);
OMG=zeros(1,length-1);
Ud=zeros(1,length-1);
PHI=zeros(1,length-1);

JUMP=zeros(1,length-1);
rp=0;
ud1=0;
omge1=0;

UD=zeros(1,4);

Fi=Fi+389;

for m=1:length-1
    data=cdata(m*5000-4999:m*5000);
    [d,CTcoh,CM,rp,omge,ud,jump]=trace(data,CB,Fi,Thetai,Phi,CTcoh,CM,rp,omge1,ud1);
    P(m)=rp;
    Ud(m)=ud;
    
    %omge=0;
    %{
    n=mod(m,4);
    UD(n+1)=ud;
    
    if(n==3)
        tmp=(abs(UD(2)-UD(3))<pi/2*1e3);
        tmp=tmp+(abs(UD(4)-UD(3))<pi/2*1e3)*2;
        tmp=tmp+(abs(UD(4)-UD(2))<pi/2*1e3)*4;
        if(tmp==0)
            omge=0;
        elseif(tmp==1)
            omge=(UD(2)+UD(3))/2;
        elseif(tmp==2)
            omge=(UD(3)+UD(4))/2;
        elseif(tmp==4)
            omge=(UD(2)+UD(4))/2;
        elseif(tmp==7)
            omge=(UD(2)+UD(3)+UD(4))/3;
        end
    end
    %}
    OMG(m)=omge;
    Phi=Phi+d;
    PHI(m)=d;
    JUMP(m)=jump;
        Fi=Fi+omge;
    
    Thetai=2*pi*Fi*5000/Fs+Thetai;
    ud1=ud;
    omge1=omge;
    
end
figure(2)
plot(P(500:end),'.');
figure(3);
hold off
plot(Ud(2:end),'r');
hold on;
plot(OMG(2:end));
figure(4);plot(PHI);
figure(5)
JUMPS=zeros(1,20);
n=1;
for m=1:20:length-20
    JUMPS=JUMPS+JUMP(m:m+19);
end
bar(JUMPS);
