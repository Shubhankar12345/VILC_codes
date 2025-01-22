%% Initlisation 
clear all;close all;clc
%% Tragectory Genration
dt=1/100;

t=0:1/100:20+2*dt;      t2=0:1/100:20;
tau1=0*t'; tau2=0*t';
TAU1=[t' tau1]; TAU2=[t' tau2];
l1 = 1 ; l2 = 1 ;

x=0*t;
t1=t-max(t)/2;
y=((tanh(t1)-min(tanh(t1)))/2)*0.2 + 0.2;

XD=[t2' x(1,1:end-2)']; YD=[t2' y(1,1:end-2)'];

[q1,q2]=tracgen(x,y);

dq1=diff(q1')/dt; dq2=diff(q2')/dt;
ddq1=diff(dq1')/dt; ddq2=diff(dq2')/dt;
q1=q1(1,1:end-2); q1=q1'; q2=q2(1,1:end-2); q2=q2';

dq1=dq1(1:end-1,1);
dq2=dq2(1:end-1,1);

ddq1=ddq1';
ddq2=ddq2';

q1d=[t2' q1];
q2d=[t2' q2];
    
th1i=q1(1);
th2i=q2(1);
%% Iteration Control
ICN=50;
BC=ICN;
PC=BC+ICN;
WC=PC+ICN;
RPC=WC+ICN;
ITRN=RPC;
%% Rates
RR_ILC=0.7;  LR_ILC=0.1;

RR_THM=0.5;  LR_THM=0.1;

RR_FF=0.7;   LR_FF=0.1;

RR_KA=0.95;  LR_KA=0.6;

RR_DA=0.95;  LR_DA=0.6;

RR_TOM=0.9;  LR_TOM=0.1;
%% Variable Intialisation
tau1=0*t2';
tau2=0*t2';
u1=0;m1=0;u_ilc=0;
u2=0;m2=0;u_ff=0;
q1m=q1*0;q2m=q2*0;
thetaM=0;SimOUT.TH1.Data=0*q1;SimOUT.dTH1.Data=0*dq1;tauMF=0;
tauT=[q1*0 q2*0];SimOUT.TH2.Data=0*q1;SimOUT.dTH2.Data=0*dq1;
TE=0;
% Stiffness Variables 
Kp=[2 0;0 2];Dp=[1 0;0 1];
Ka=Kp*0; Da=Dp*0;
K=Kp+Ka;
D=Dp+Da;
F=0;

figure(1)
SimOUT=sim('Impdc_CNTRL_v1.slx');
for i=1:1:ITRN
    disp(i)
    th1i=q1(1) +  0*rand/100; % Change in initial condition
    th2i=q2(1) +  0*rand/100; % Change in initial condition
    % Input to the system
    TAU1=[t2' tau1];
    TAU2=[t2' tau2];    
    if(i <= BC)         
        F=0; %Perturbation Flag off
        %Result Plotter
        figure(1)
        subplot(141)
        xfk = l1*cos(SimOUT.TH1.Data) + l2*cos(SimOUT.TH1.Data + SimOUT.TH2.Data) ;
        yfk = l1*sin(SimOUT.TH1.Data) + l2*sin(SimOUT.TH1.Data + SimOUT.TH1.Data) ;
        axis([-0.03 0.03 0.15 0.45]);
%         plot(t2,SimOUT.TH1.Data)
        plot(xfk,yfk) 
        hold on
%         subplot(422)
%         plot(t2,SimOUT.TH2.Data)
%         hold on
        Bth1=SimOUT.TH1.Data;
        Bth2=SimOUT.TH2.Data;
        if(i == 1)
            figure(2)
            subplot(441)
            [Xko,Yko]=ellipsoids(K,0,0);
            plot(Xko,Yko,'b','Linewidth',2);axis square;axis([-15 15 -15 15]);
            subplot(449)
            [XDo,YDo]=ellipsoids(D,0,0);
            plot(XDo,YDo,'b','Linewidth',2);axis square;axis([-15 15 -15 15]);
        end
        hold on
        if(i == BC)
            figure(2)
            subplot(4,4,5)
            [Xko,Yko]=ellipsoids(K,0,0);
            plot(Xko,Yko,'r','Linewidth',2);axis square;axis([-15 15 -15 15]);
            subplot(4,4,13)
            [XDo,YDo]=ellipsoids(D,0,0);
            plot(XDo,YDo,'r','Linewidth',2);axis square;axis([-15 15 -15 15]);
        end
    end
    if(i > BC && i<= PC) %(i > ITRN/4 && i<= 2*ITRN/4)
        F=1; %Perturbation Flag ON        
        % Result Plotter
        xfk = l1*cos(SimOUT.TH1.Data) + l2*cos(SimOUT.TH1.Data + SimOUT.TH2.Data) ;
        yfk = l1*sin(SimOUT.TH1.Data) + l2*sin(SimOUT.TH1.Data + SimOUT.TH1.Data) ;
        figure(1)
        subplot(142)
        plot(xfk,yfk) 
        axis([-0.03 0.03 0.15 0.45])
        hold on
%         subplot(424)
%         plot(t2,SimOUT.TH2.Data)
%         hold on
        Pth1=SimOUT.TH1.Data;
        Pth2=SimOUT.TH2.Data;
        if(i == BC+1)
            figure(2)
            subplot(442)
            [Xko,Yko]=ellipsoids(K,0,0);
            plot(Xko,Yko,'b','Linewidth',2);axis square;axis([-15 15 -15 15]);
            subplot(4,4,10)
            [XDo,YDo]=ellipsoids(D,0,0);
            plot(XDo,YDo,'b','Linewidth',2);axis square;axis([-15 15 -15 15]);
        end
        hold on
        if(i == PC)
            figure(2)
            subplot(4,4,6)
            [Xko,Yko]=ellipsoids(K,0,0);
            plot(Xko,Yko,'r','Linewidth',2);axis square;axis([-15 15 -15 15]);
            subplot(4,4,14)
            [Xdo,Ydo]=ellipsoids(D,0,0);
            plot(Xdo,Ydo,'r','Linewidth',2);axis square;axis([-15 15 -15 15]);
        end
    end
    if(i > PC && i<= WC ) %(i > 2*ITRN/4 && i<= 3*ITRN/4 )
     
        F=0; %Perturbation flag OFF

        % Result Plotter
        figure(1)
        subplot(143)
        xfk = l1*cos(SimOUT.TH1.Data) + l2*cos(SimOUT.TH1.Data + SimOUT.TH2.Data) ;
        yfk = l1*sin(SimOUT.TH1.Data) + l2*sin(SimOUT.TH1.Data + SimOUT.TH1.Data) ;
        axis([-0.03 0.03 0.15 0.45]) ;
        plot(xfk,yfk) ;
        hold on
%         subplot(426)
%         plot(t2,SimOUT.TH2.Data)
%         hold on
        Wth1=SimOUT.TH1.Data;
        Wth2=SimOUT.TH2.Data;
        if(i == PC+1)
            figure(2)
            subplot(443)
            [Xko,Yko]=ellipsoids(K,0,0);
            plot(Xko,Yko,'b','Linewidth',2);axis square;axis([-15 15 -15 15]);
            subplot(4,4,11)
            [Xdo,Ydo]=ellipsoids(D,0,0);
            plot(Xdo,Ydo,'b','Linewidth',2);axis square;axis([-15 15 -15 15]);
        end
        hold on
        if(i == WC)
            figure(2)
            subplot(4,4,7)
            [Xko,Yko]=ellipsoids(K,0,0);
            plot(Xko,Yko,'r','Linewidth',2);axis square;axis([-15 15 -15 15]);
            subplot(4,4,15)
            [Xdo,Ydo]=ellipsoids(D,0,0);
            plot(Xdo,Ydo,'r','Linewidth',2);axis square;axis([-15 15 -15 15]);
        end
    end
    if( i> WC)    
        F=1; %Perturbation flag ON
        figure(1)
        xfk = l1*cos(SimOUT.TH1.Data) + l2*cos(SimOUT.TH1.Data + SimOUT.TH2.Data) ;
        yfk = l1*sin(SimOUT.TH1.Data) + l2*sin(SimOUT.TH1.Data + SimOUT.TH1.Data) ;
        subplot(144)
        plot(xfk,yfk) ;
        axis([-0.03 0.03 0.15 0.45]) ;
        hold on
%         subplot(428)
%         plot(t2,SimOUT.TH2.Data)
%         hold on
        RPth1=SimOUT.TH1.Data;
        RPth2=SimOUT.TH2.Data;
        if(i == WC+1)
            figure(2)
            subplot(444)
            [Xko,Yko]=ellipsoids(K,0,0);
            plot(Xko,Yko,'b','Linewidth',2);axis square;axis([-15 15 -15 15]);
            subplot(4,4,12)
            [Xdo,Ydo]=ellipsoids(D,0,0);
            plot(Xdo,Ydo,'b','Linewidth',2);axis square;axis([-15 15 -15 15]);

        end
        hold on
        if(i == RPC)
            figure(2)
            subplot(4,4,8)
            [Xko,Yko]=ellipsoids(K,0,0);
            plot(Xko,Yko,'r','Linewidth',2);axis square;axis([-15 15 -15 15]);
            subplot(4,4,16)
            [Xdo,Ydo]=ellipsoids(D,0,0);
            plot(Xdo,Ydo,'r','Linewidth',2);axis square;axis([-15 15 -15 15]);
        end
    end
    SimOUT=sim('Impdc_CNTRL_v1.slx');
    %% Stiffness Update 

    
%     nDelta=-1*pinv([SimOUT.TH1.Data-q1 SimOUT.TH2.Data-q2 SimOUT.dTH1.Data-dq1 SimOUT.dTH2.Data-dq2])*(tauT-tauMF);

    DeltaK=-1*pinv([SimOUT.TH1.Data-q1 SimOUT.TH2.Data-q2])*(tauT-tauMF);
    DeltaD=-1*pinv([SimOUT.dTH1.Data-dq1 SimOUT.dTH2.Data-dq2])*(tauT-tauMF);

    Ka=RR_KA*Ka + LR_KA*tanh(DeltaK);
    Da=RR_DA*Da + LR_DA*tanh(DeltaD);
    K=Ka+Kp;
    D=Da+Dp;
   
    %% Model Update 
    tauOP=[tau1';tau2'];
    thetaC = [SimOUT.TH1.Data SimOUT.TH2.Data SimOUT.dTH1.Data SimOUT.dTH2.Data SimOUT.ddTH1.Data SimOUT.ddTH2.Data];
    tau_ic=[SimOUT.tau1.Data SimOUT.tau2.Data];
    tauMF=RR_TOM*tauMF + LR_TOM*((tauT) -tauMF);
    %% Errors & UPDATE
    %Task error
    TE=[(q1 - SimOUT.TH1.Data) (q2 - SimOUT.TH2.Data)];
    %SPE Error
    thetaM=RR_THM*thetaM + LR_THM*(thetaC-thetaM);
    % MOE
    PP=fwdynamics(tauOP, [SimOUT.TH1.Data(1);SimOUT.TH2.Data(1)], [SimOUT.dTH1.Data(1);SimOUT.dTH2.Data(1)], [SimOUT.ddTH1.Data(1);SimOUT.ddTH2.Data(1)], dt);
    MOE=(fwdynamics(tauOP, [SimOUT.TH1.Data(1);SimOUT.TH2.Data(1)], [SimOUT.dTH1.Data(1);SimOUT.dTH2.Data(1)], [SimOUT.ddTH1.Data(1);SimOUT.ddTH2.Data(1)], dt) - thetaM);
    % Update Equations
    u_ilc=RR_ILC*u_ilc + LR_ILC*TE;
    %Feedforward
    MOE_thd=MOE(:,1:2);
    MOE_dthd=MOE(:,3:4);
    MOE_ddthd=MOE(:,5:6);

    u_ff=RR_FF*u_ff + (LR_FF*ivdynamics(MOE_thd',MOE_dthd',MOE_ddthd')');  
    tauT=u_ilc+u_ff;    
    %%tauMF=tauT;
    tau1=tauT(:,1);
    tau2=tauT(:,2);
    E(i,1)=norm(q1-SimOUT.TH1.Data)+norm(q2-SimOUT.TH2.Data);
    E(i,2)=norm(u_ilc + u_ff);
    E(i,3)=trace(K);
    E(i,4)=trace(D);
    clc
end
figure(2)
subplot(441)
grid minor
xlabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
subplot(442)
grid minor
xlabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
subplot(443)
grid minor
xlabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
subplot(444)
grid minor
xlabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
subplot(445)
grid minor
xlabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
subplot(446)
grid minor
xlabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
subplot(447)
grid minor
xlabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nm}{rad}$','fontsize',18,'interpreter','latex')
subplot(448)
grid minor
xlabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
subplot(449)
grid minor
xlabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
subplot(4,4,10)
grid minor
xlabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
subplot(4,4,11)
grid minor
xlabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
subplot(4,4,12)
grid minor
xlabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
subplot(4,4,13)
grid minor
xlabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
subplot(4,4,14)
grid minor
xlabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
subplot(4,4,15)
grid minor
xlabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
subplot(4,4,16)
grid minor
xlabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
ylabel('$\frac{Nms}{rad}$','fontsize',18,'interpreter','latex')
orient(gcf,'landscape')
print('-painters','-fillpage','-dpdf','Result_Opti_Imped_Stiff_Space')

figure(1)
subplot(141)
ylim([0 0.3])
grid minor
plot(l1*cos(q1)+l2*cos(q1+q2),l1*sin(q1)+l2*sin(q1+q2),'--r','linewidth',2); hold on; plot(l1*cos(Bth1)+l2*cos(Bth1+Bth2),l1*sin(Bth1)+l2*sin(Bth1+Bth2),'--g','linewidth',2);
axis([-0.03 0.03 0.15 0.45]); 
xlabel('Time (s)','fontsize',18,'interpreter','latex')
ylabel('rad','fontsize',18,'interpreter','latex')
title('Baseline','fontsize',18,'interpreter','latex')
subplot(142)
ylim([2.7 3])
grid minor
plot(l1*cos(q1)+l2*cos(q1+q2),l1*sin(q1)+l2*sin(q1+q2),'--r','linewidth',2);hold on;plot(l1*cos(Pth1)+l2*cos(Pth1+Pth2),l1*sin(Pth1)+l2*sin(Pth1+Pth2),'--g','linewidth',2);
axis([-0.03 0.03 0.15 0.45]); 
xlabel('Time (s)','fontsize',18,'interpreter','latex')
ylabel('rad','fontsize',18,'interpreter','latex')
title('Perturbation','fontsize',18,'interpreter','latex')
subplot(143)
ylim([0 0.3])
grid minor
plot(l1*cos(q1)+l2*cos(q1+q2),l1*sin(q1)+l2*sin(q1+q2),'--r','linewidth',2);hold on;plot(l1*cos(Wth1)+l2*cos(Wth1+Wth2),l1*sin(Wth1)+l2*sin(Wth1+Wth2),'--g','linewidth',2);
axis([-0.03 0.03 0.15 0.45]); 
xlabel('Time (s)','fontsize',18,'interpreter','latex')
ylabel('rad','fontsize',18,'interpreter','latex')
title('Washout','fontsize',18,'interpreter','latex')
subplot(144)
ylim([2.7 3])
grid minor
plot(l1*cos(q1)+l2*cos(q1+q2),l1*sin(q1)+l2*sin(q1+q2),'--r','linewidth',2);hold on;plot(l1*cos(RPth1)+l2*cos(RPth1+RPth2),l1*sin(RPth1)+l2*sin(RPth1+RPth2),'--g','linewidth',2);
axis([-0.03 0.03 0.15 0.45]); 
xlabel('Time (s)','fontsize',18,'interpreter','latex')
ylabel('rad','fontsize',18,'interpreter','latex')
title('Reperturbation','fontsize',18,'interpreter','latex')
orient(gcf,'landscape')
print('-painters','-fillpage','-dpdf','Result_Opti_Imped_ON_Joint_SPACE')
figure(3)
subplot(311)
stem(E(:,1))
grid minor
xlabel('Iterations','interpreter','latex','fontsize',18)
ylabel('Norm of Error','interpreter','latex','fontsize',18)
title('Error Plot of RR SCM','interpreter','latex','fontsize',18)

subplot(312)
stem(E(:,3))
grid minor
xlabel('Iterations','interpreter','latex','fontsize',18)
ylabel('Norm of K','interpreter','latex','fontsize',18)
title('Error Plot of RR SCM','interpreter','latex','fontsize',18)

subplot(313)
stem(E(:,4))
grid minor
xlabel('Iterations','interpreter','latex','fontsize',18)
ylabel('Norm of D','interpreter','latex','fontsize',18)
title('Error Plot of RR SCM','interpreter','latex','fontsize',18)

orient(gcf,'landscape')
print('-painters','-fillpage','-dpdf','Result_Opti_Imped_ON_Error_Space')

% figure(5)
% subplot(141)
% grid minor
% xlabel('$\frac{N}{m}$','fontsize',18,'interpreter','latex')
% ylabel('$\frac{N}{m}$','fontsize',18,'interpreter','latex')
% subplot(142)
% grid minor
% xlabel('$\frac{N}{m}$','fontsize',18,'interpreter','latex')
% ylabel('$\frac{N}{m}$','fontsize',18,'interpreter','latex')
% subplot(143)
% grid minor
% xlabel('$\frac{N}{m}$','fontsize',18,'interpreter','latex')
% ylabel('$\frac{N}{m}$','fontsize',18,'interpreter','latex')
% subplot(144)
% grid minor
% xlabel('$\frac{N}{m}$','fontsize',18,'interpreter','latex')
% ylabel('$\frac{N}{m}$','fontsize',18,'interpreter','latex')
% orient(gcf,'landscape')
% print('-painters','-fillpage','-dpdf','Result_Opti_Imped_ON_Cartesian_Space')

% subplot(212)
% stem(E(:,2),'r')
% hold on
% stem(E(:,3),'b')
% stem(E(:,4),'k')
% grid minor
% xlabel('Iterations','interpreter','latex','fontsize',18)
% ylabel('Norm of Error','interpreter','latex','fontsize',18)
% title('Error Plot of RR SCM','interpreter','latex','fontsize',18)
% legend('ILC','IMC','Total','Location','northwest','interpreter','latex','fontsize',18)

function [q1,q2]=tracgen(x,y)
l1=1;
l2=1;

q2=acos((x.^2 + y.^2 - l1.^2 - l2.^2)/(2*l1*l2));
q1=atan(y./x) - atan(l2*sin(q2)./(l1 + (l2*cos(q2))));
q2(end)=q2(end-1);

for i=1:1:length(q1)
    if(isnan(q1(i))==1)
        q1(i)=0;
    end
    
    if(isnan(q2(i))==1)
        q2(i)=0;
    end
end
end

function [Xo,Yo]=ellipsoids(J,x0,y0)
[evt,evl]=eig(J*J');
if(evl(2,2)>evl(1,1))
    a=sqrt(evl(2,2));
    b=sqrt(evl(1,1));
    ax=1;
    ay=0;
    bx=evt(1,2);
    by=evt(2,2);
    th=atan2( ax*by-ay*bx, ax*bx+ay*by );
end
if(evl(2,2)<=evl(1,1))
    a=sqrt(evl(1,1));
    b=sqrt(evl(2,2));
    ax=1;
    ay=0;
    bx=evt(1,2);
    by=evt(2,2);
    th=atan2( ax*by-ay*bx, ax*bx+ay*by );
end
thf= th + pi/2;
te=-pi:0.01:pi;
x56=a*cos(te);
y56=b*sin(te);
R5=[cos(th) -sin(th) ; sin(th) cos(th)];
P5= R5*[x56;y56];
Xo= P5(1,:)+x0;
Yo= P5(2,:)+y0;
R5=[ cos(thf) -sin(thf) ; sin(thf) cos(thf)];
P5= R5*[x56;y56];
end

function theta_hat = fwdynamics(Tau, thi, dthi, ddthi, dt)
    m1 = 1 ;
    m2 = 1 ; g = 0 ;
    l1 = 1 ; lc1 = l1/2 ;
    l2 = 1 ; lc2 = l2/2 ;
    I1 = 1/12 ; I2 = 1/12 ;    
    b1 = 0.5 ; b2 = 0.5 ;
    th_b = [thi] ;
    dth_b = [dthi] ;
    ddth_b = [ddthi] ;
    for i = 1:length(Tau)-1
        M = [m1*lc1^2 + m2*(l1^2+lc2^2+2*l1*lc2^2+2*l1*lc2*cos(th_b(2,i))) + I1 + I2, m2*(lc2^2 + l1*lc2*cos(th_b(2,i))) + I2;
            m2*(lc2^2 + l1*lc2*cos(th_b(2,i))) + I2, m2*lc2^2 + I2] ;
        C = [-m2*l1*lc2*sin(th_b(2,i))*dth_b(2,i), -m2*l1*lc2*sin(th_b(2,i))*dth_b(2,i) - m2*l1*lc2*sin(th_b(2,i))*dth_b(1,i);
            m2*l1*lc2*sin(th_b(2,i))*dth_b(1,i), 0] ;
        G = [(m1*lc1 + m2*l1)*g*cos(th_b(1,i)) + m2*lc2*g*cos(th_b(1,i) + th_b(2,i));
            m2*lc2*g*cos(th_b(1,i)+th_b(2,i))] ;
        D = [b1*dth_b(1,i);b2*dth_b(2,i)] ;
        ddth_b(:,i+1) = inv(M)*(Tau(:,i) - C*dth_b(:,i) - G - D) ;
        dth_b(:,i+1) = dth_b(:,i) + ddth_b(:,i+1).*dt ;
        th_b(:,i+1) = th_b(:,i) + dth_b(:,i+1).*dt ;
    end    
    theta_hat = [th_b(1,:)' th_b(2,:)' dth_b(1,:)' dth_b(2,:)' ddth_b(1,:)' ddth_b(2,:)'] ;
end


function  Tau_iv = ivdynamics(thd,dthd,ddthd)
    m1 = 1 ;
    m2 = 1 ; g = 0 ;
    l1 = 1 ; lc1 = l1/2 ;
    l2 = 1 ; lc2 = l2/2 ;
    I1 = 1/12 ; I2 = 1/12 ;    
    b1 = 0.5 ; b2 = 0.5 ;
%     pi1 = m1*lc1^2 + m2*(l1^2 + lc2^2) + I1 + I2 ;
%     pi2 = m2*lc2*l1 ;
%     pi3 = m2*l1*lc2 ;
%     pi4 = m1*lc1 + m2*l1 ;
%     pi5 = m2*l2 ;
%     pii = [pi1;pi2;pi3;pi4;pi5] ;
    for i = 1:length(thd(1,:))
%         y11 = ddthd(1,i) ;
%         y12 = cos(thd(2,i))*(2*ddthd(1,i) + ddthd(2,i)) + sin(thd(2,i))*(dthd(1,i)^2 - 2*dthd(1,i)*dthd(2,i)) ;
%         y13 = ddthd(2,i) ;
%         y14 = g*cos(thd(1,i)) ;
%         y15 = g*cos(thd(1,i) + thd(2,i)) ;
%         y21 = 0 ;
%         y22 = cos(thd(2,i))*ddthd(1,i) + sin(thd(2,i))*dthd(1,i)^2 ;
%         y23 = ddthd(2,i) ;
%         y24 = 0 ;
%         y25 = g*cos(thd(1,i) + thd(2,i)) ;
%         Y = [y11,y12,y13,y14,y15;
%             y21,y22,y23,y24,y25] ;
%         D = [b1*dthd(1,i);b2*dthd(2,i)] ;
%         Tau_iv(:,i) = Y*pii + D ;
        M = [m1*lc1^2 + m2*(l1^2+lc2^2+2*l1*lc2^2+2*l1*lc2*cos(thd(2,i))) + I1 + I2, m2*(lc2^2 + l1*lc2*cos(thd(2,i))) + I2;
            m2*(lc2^2 + l1*lc2*cos(thd(2,i))) + I2, m2*lc2^2 + I2] ;
        C = [-m2*l1*lc2*sin(thd(2,i))*dthd(2,i), -m2*l1*lc2*sin(thd(2,i))*dthd(2,i) - m2*l1*lc2*sin(thd(2,i))*dthd(1,i);
            m2*l1*lc2*sin(thd(2,i))*dthd(1,i), 0] ;
        G = [(m1*lc1 + m2*l1)*g*cos(thd(1,i)) + m2*lc2*g*cos(thd(1,i) + thd(2,i));
            m2*lc2*g*cos(thd(1,i)+thd(2,i))] ;
        D = [b1*dthd(1,i);b2*dthd(2,i)] ;  
        Tau_iv(:,i) = M*[ddthd(1,i);ddthd(2,i)] + C*[dthd(1,i);dthd(2,i)] + G + D ;
    end
end
