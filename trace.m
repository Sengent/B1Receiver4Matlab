function [d,CTcoh,CM,rp,omge,ud,jump]=trace(data,CA,Fi,Thetai,Phi,CTcoh,CM,rp,omge1,ud1)
    Fs=5e6;
    t=(0:4999)/Fs;
   
    CAWidth=Fs/(2.046e6);
    
    %混频
    u=data.*exp(-1j*(2*pi*Fi*t+Thetai));
    
    %本地C/A码
    CA_p=CA(t2n(t,Fs,Phi,CAWidth));
    CA_l=CA(t2n(t-0.5*CAWidth/Fs,Fs,Phi,CAWidth));
    CA_e=CA(t2n(t+0.5*CAWidth/Fs,Fs,Phi,CAWidth));
    
    %{
    % 
    figure(1)
    hold off;    plot(real(u(2600:2850))); hold on;plot(CA_p(2600:2850),'g')
    %}
    
    %相关
    p=sum(double(CA_p).*u);
    l=sum(double(CA_l).*u);
    e=sum(double(CA_e).*u);
    %hold on;    plot(p);
    %幅值
    P=abs(p);
    Lt=abs(l);Et=abs(e);
    CTcoh=CTcoh-1;
    %TODO 积分清除
    if(CTcoh==0)
        E=CM(1)+Et;L=CM(2)+Lt;
        CM(:)=0;
        CTcoh(:)=20;
        %鉴别器
        d=(E-L)/(E+L)/2;
    else
        CM(1)=CM(1)+Et;
        CM(2)=CM(2)+Lt;
        %(Et-Lt)/(Et+Lt)
        d=0;
    end    
    
    %锁频环 鉴频器
    omge=0;
    n=0;
    ud=0;
    jump=abs(angle(p)-angle(rp))>pi;
    if(rp~=0)  
        %{
        Ip=real(p);
        Qp=imag(p);
        Ip1=real(rp);
        Qp1=imag(rp);

        Pdot=Ip1*Ip+Qp1*Qp;
        Pcross=Ip1*Qp-Qp1*Ip;
        %}
        tmp=p*conj(rp);
        tmp=tmp/abs(tmp);
        Pdot=real(tmp);
        Pcross=imag(tmp);

        ud=Pcross*sign(Pdot);
        %ud=atan2(Pcross,Pdot);
        %ud=Pcross*sign(Pdot);
        %{
        if(ud>pi/2)
            ud=ud-pi;
        elseif(ud<-pi/2)
            ud=ud+pi;
        end
        %}
        ud=ud*1e3/2/pi;
        
        %跳跃点检测
        %{
        if(abs(ud-ud1)>0.75*pi*1e3)
            %ud=ud1;
            omge=0;
            drop=1;
            return
        end
        %}
        
        %环路滤波器
        K=20;
        omgn=7.547;
        a2=1.414;
        
        b0=a2*omgn+(1e-3)/2*omgn^2;
        b1=-a2*omgn+(1e-3)/2*omgn^2;
               
        omge=omge1+1/K*(b0*ud+b1*ud1);
        
        %omge=ud/80;
    end

    rp=p;
end

function n=t2n(t,Fs,Phi,CAWidth)
    n=mod(ceil((t*Fs+1+Phi)/CAWidth),2046);
    n(n==0)=2046;
end