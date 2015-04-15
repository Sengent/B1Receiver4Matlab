clear
load('data.mat');
load('CBlist.mat');
CB=CBs(7,:);
NH=[0, 0, 0, 0, 0, 1, 0, 0, 1, 1,0, 1, 0, 1, 0, 0, 1, 1, 1, 0]*2-1;

Fs=5e6;%Hz
test_time=4000;%ms
fll_time=1000;
freq_0=0;%中频
CB_width=Fs/(2.046e6);
bi=763;

% DLL 环路滤波器参数
B_L_dll=0.1;
omg_N_dll=B_L_dll*4;
T=1e-3;
T_coh=40;
% FLL 环路滤波器参数
B_L_fll=10;
omg_N_fll=B_L_fll/0.53;
a2=1.414;
% PLL 环路滤波器参数
B_L_pll=10;
omg_N_pll=B_L_pll/0.7845;
a3=1.1;
b3=2.4;

% catch
[freq_i,code_phase,rate]=catchB1(CB,cdata(1:5e3));

m=code_phase;
code_phase_0=0;
code_phase=0;
theta_i=0;
theta_0=0;
t=(0:4999)/Fs;

p=zeros(1,test_time);
e=zeros(1,test_time);
l=zeros(1,test_time);
code_phase_e=zeros(1,floor(test_time/T_coh)+1);
w_dll=zeros(1,floor(test_time/T_coh)+1);
code_phase_e_dll=zeros(1,floor(test_time/T_coh)+1);
freq_e=zeros(1,fll_time);
freq_fll=zeros(1,fll_time);
w_1_fll=zeros(1,fll_time);
w_2_fll=zeros(1,fll_time);
freq_fll(2)=freq_i;
w_2_fll(2)=freq_i;
%code_lock_state=zeros(1,test_time);
%freq_state=zeros(1,fll_time);
%phase_lock_state=zeros(1,test_time);
fll_state=true;
phase_e=zeros(1,test_time-fll_time);
phase_pll=zeros(1,test_time-fll_time);
w_1_pll=zeros(1,test_time-fll_time);
w_2_pll=zeros(1,test_time-fll_time);
%w_1_pll(1)=

n=2;
nn=2;
phase=0;
while m<test_time*5e3-5e3
    data=cdata(m+1:m+5e3);
    %CB_width=Fs/(1.023e6*(1+1/bi));
    % 混频
    u=data.*exp(-1j*(2*pi*freq_i*t+theta_i));
    %% Code Tracking Loop
    % 复制本地CA码
    n_c=mod(floor((t*Fs+code_phase)/CB_width+1),1023);
    n_c(n_c==0)=1023;
    cb_p=CB(n_c);
    cb_e=cb_p([3:end,1,2]);
    cb_l=cb_p([end-1,end,1:end-2]);
    % 相关器
    p(n)=sum(double(cb_p).*u);
    l(n)=sum(double(cb_l).*u);
    e(n)=sum(double(cb_e).*u);
    % Code Discriminator
    if(mod(n,T_coh)==0)
        Lt=sum(abs(l(n-T_coh+1:n)));    Et=sum(abs(e(n-T_coh+1:n)));
        %Pt=abs(p(n));
        nnn=floor(n/T_coh)+1;
        code_phase_e(nnn)=(Et-Lt)/(Et+Lt)/2;
        % Code Loop Filter
        w_dll(nnn)=w_dll(nnn-1)+2*omg_N_dll^2*T*code_phase_e(nnn);
        code_phase_e_dll(nnn)=w_dll(nnn)+2*omg_N_dll*code_phase_e(nnn);
    end
    
    % Code Lock Detector    
    if(n<=2)
        
    elseif(fll_state)
       %% Frequency Lock Loop
        tmp=p(n)*conj(p(n-1));
        tmp=tmp/abs(tmp);
        dot=real(tmp);
        cross=imag(tmp);
        % FLL Discriminator
        freq_e(n)=cross*sign(dot)*1e3/4;
        %freq_e(n)=-atan2(dot,cross)*1e3/2/pi;
        % FLL filter
        w_1_fll(n)=w_1_fll(n-1)+omg_N_fll^2*T*freq_e(n);
        w_2_fll(n)=w_2_fll(n-1)+T*((w_1_fll(n)+w_1_fll(n-1))/2+a2*omg_N_fll*freq_e(n));
        freq_fll(n)=(w_2_fll(n)+w_2_fll(n-1))/2;
        %w_1_fll(n)=w_1_fll(n-1)+omg_N_fll^2*T*freq_e(n);
        %freq_e_fll(n)=freq_e_fll(n-1)+w_1_fll(n)+2*cacy*omg_N_fll*freq_e(n);
        freq_i=freq_0+freq_fll(n);  
        if(n>=fll_time)
            fll_state=false;
            %w_2_pll(2)=angle(p(n))+5e3/Fs*2*pi*freq_i;
            theta_0=angle(p(n));
            theat_i=theta_0;
            nn=2;
        end
    else
       %% Phase Lock Loop
        I_ps=real(p(n));
        Q_ps=imag(p(n));
        p_angle=atan2(Q_ps,I_ps);
        if p_angle<-pi/2;
            p_angle=p_angle+pi;
        elseif p_angle>pi/2
            p_angle=p_angle-pi;
        end
        phase_e(nn)=p_angle;
        % PLL filter
        w_1_pll(nn)=w_1_pll(nn-1)+omg_N_pll^3*T*phase_e(nn);
        w_2_pll(nn)=w_2_pll(nn-1)+((w_1_pll(nn)+w_1_pll(nn-1))/2+a3*omg_N_pll^2*phase_e(nn))*T;
        phase_pll(nn)=(w_2_pll(nn)+w_2_pll(nn-1))/2;
        %theta_i=theta_0+phase_pll(n);
        % 锁定检测
        phase=phase_pll(nn);
        nn=nn+1;
    end
   
    
    code_phase=code_phase-freq_i/bi/1000;
    %code_phase=code_phase_0+code_phase_e_dll(n);
    if(mod(n,T_coh)==0)
        code_phase=code_phase+code_phase_e_dll(floor(n/T_coh)+1);
    end
    %code_phase=code_phase_0;
    
    %区间调整
    if code_phase>0.8
        code_phase=code_phase-1;
        m=m-1;
    elseif code_phase<-0.8
        code_phase=code_phase+1;
        m=m+1;
    end
    
    m=m+5e3;
    theta_0=theta_0+5e3/Fs*2*pi*freq_i;
    theta_i=theta_0+phase;
    %theta_i=mod(theta_i,pi*2);
    
    n=n+1;    
end
figure(1)
hold off;plot(code_phase_e/Fs*3e8)
hold on;plot(code_phase_e_dll/Fs*3e8,'r')

figure(2)
%hold off;plot(freq_e);
hold on;plot(freq_fll,'r');
figure(3)
hold off;plot(phase_e);
 phase=mod(phase_pll,pi);
 phase(phase>pi/2)=phase(phase>pi/2)-pi;
 hold on;plot(phase,'r');
%hold on;plot(phase_pll,'g');
