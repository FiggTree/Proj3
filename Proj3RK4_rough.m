clear
Qval(1)=700000;
CAval(1)=.25;
Tval(1)=350;
Tsval(1)=350;
timeval(1)=0;
a=timeval(1);
b=100;
h=.25;
N=(b-a)./h;

tH=10;
tTs=10;
F=10000;
Qspec=700000;
Qspec2=900000;

for i=1:N
    timeval(i+1)=timeval(i)+h;
    K1=dQ(Qval(i),Qspec,tH);
    K2=dQ(Qval(i)+K1.*h./2,Qspec,tH);
    K3=dQ(Qval(i)+K2.*h./2,Qspec,tH);
    K4=dQ(Qval(i)+K3.*h,Qspec,tH);
    Qval(i+1)=Qval(i) + (K1+2.*K2+2.*K3+K4).*h./6;
    K1=dCA(CAval(i),Tval(i),F);
    K2=dCA(CAval(i)+K1.*h./2,Tval(i),F);
    K3=dCA(CAval(i)+K2.*h./2,Tval(i),F);
    K4=dCA(CAval(i)+K3.*h,Tval(i),F);
    CAval(i+1)=CAval(i) + (K1+2.*K2+2.*K3+K4).*h./6;
    K1=dT(Tval(i),CAval(i),Qval(i),F);
    K2=dT(Tval(i)+K1.*h./2,CAval(i),Qval(i),F);
    K3=dT(Tval(i)+K2.*h./2,CAval(i),Qval(i),F);
    K4=dT(Tval(i)+K3.*h,CAval(i),Qval(i),F);
    Tval(i+1)=Tval(i) + (K1+2.*K2+2.*K3+K4).*h./6;
    K1=dTs(Tsval(i),Tval(i),tTs);
    K2=dTs(Tsval(i)+K1.*h./2,Tval(i),tTs);
    K3=dTs(Tsval(i)+K2.*h./2,Tval(i),tTs);
    K4=dTs(Tsval(i)+K3.*h,Tval(i),tTs);
    Tsval(i+1)=Tsval(i) + (K1+2.*K2+2.*K3+K4).*h./6;
end

%vary matrix F+50%, F-50%, tH+50%, tH-50%, tTs+50%, tTs-50%
QvalVary(1:6,1)=700000;
CAvalVary(1:6,1)=.25;
TvalVary(1:6,1)=350;
TsvalVary(1:6,1)=350;
FVary=[F*1.5;F*.5;F;F;F;F];
tHVary=[tH;tH;tH*1.5;tH*.5;tH;tH];
tTsVary=[tTs;tTs;tTs;tTs;tTs*1.5;tTs*.5];
for j=1:6
    for i=1:N
        K1=dQ(QvalVary(j,i),Qspec,tHVary(j));
        K2=dQ(QvalVary(j,i)+K1.*h./2,Qspec,tHVary(j));
        K3=dQ(QvalVary(j,i)+K2.*h./2,Qspec,tHVary(j));
        K4=dQ(QvalVary(j,i)+K3.*h,Qspec,tHVary(j));
        QvalVary(j,i+1)=QvalVary(j,i)+(K1+2.*K2+2.*K3+K4).*h./6;
        K1=dCA(CAvalVary(j,i),TvalVary(j,i),FVary(j));
        K2=dCA(CAvalVary(j,i)+K1.*h./2,TvalVary(j,i),FVary(j));
        K3=dCA(CAvalVary(j,i)+K2.*h./2,TvalVary(j,i),FVary(j));
        K4=dCA(CAvalVary(j,i)+K3.*h,TvalVary(j,i),FVary(j));
        CAvalVary(j,i+1)=CAvalVary(j,i)+(K1+2.*K2+2.*K3+K4).*h./6;
        K1=dT(TvalVary(j,i),CAvalVary(j,i),QvalVary(j,i),FVary(j));
        K2=dT(TvalVary(j,i)+K1.*h./2,CAvalVary(j,i),QvalVary(j,i),FVary(j));
        K3=dT(TvalVary(j,i)+K2.*h./2,CAvalVary(j,i),QvalVary(j,i),FVary(j));
        K4=dT(TvalVary(j,i)+K3.*h,CAvalVary(j,i),QvalVary(j,i),FVary(j));
        TvalVary(j,i+1)=TvalVary(j,i)+(K1+2.*K2+2.*K3+K4).*h./6;
        K1=dTs(TsvalVary(j,i),TvalVary(j,i),tTsVary(j));
        K2=dTs(TsvalVary(j,i)+K1.*h./2,TvalVary(j,i),tTsVary(j));
        K3=dTs(TsvalVary(j,i)+K2.*h./2,TvalVary(j,i),tTsVary(j));
        K4=dTs(TsvalVary(j,i)+K3.*h,TvalVary(j,i),tTsVary(j));
        TsvalVary(j,i+1)=TsvalVary(j,i)+(K1+2.*K2+2.*K3+K4).*h./6;
    end
end

% 4 figures holding the different differential components constant across
% different variances. 
figure
hold on
plot(timeval,QvalVary(1,:),'r') % +50% Feed
plot(timeval,QvalVary(2,:),'b') % -50% Feed
plot(timeval,QvalVary(3,:),'r.') % +50% tH
plot(timeval,QvalVary(4,:),'b.') % -50% tH
plot(timeval,QvalVary(5,:),'r--') % +50% tTs
plot(timeval,QvalVary(6,:),'b--') % -50% tTs
plot(timeval,Qval,'k') % no variance
legend('plus feed','minus feed','plus tH','minus tH','plus tTs','minus tTs','no variance');
title('Q heat added to reactor per second')
ylabel('Heat per second (cal/s)')
xlabel('Time (s)')
%no change in the amount of heat added per second becuase Q is not varied,
%and Q=Qspec the whole time making the change of Q per second per second
%zero always creating a constant rate of addition. 
figure
hold on
plot(timeval,CAvalVary(1,:),'r'); % +50% Feed
plot(timeval,CAvalVary(2,:),'b'); % -50% Feed
plot(timeval,CAvalVary(3,:),'r.'); % +50% tH
plot(timeval,CAvalVary(4,:),'b.'); % -50% tH
plot(timeval,CAvalVary(5,:),'r--'); % +50% tTs
plot(timeval,CAvalVary(6,:),'b--'); % -50% tTs
plot(timeval,CAval,'k') % no variance
legend('plus feed','minus feed','plus tH','minus tH','plus tTs','minus tTs','no variance');
ylim([.1 1.1]);
title('CA concentration of A in reactor')
ylabel('Concentration (gmol/L)')
xlabel('Time (s)')
%The rate at which the concentration changes would rely only on a changing
%feed rate and temperature. A variance of tH or tTs has no effect on the rate of change of
%concentration since temperature relies on feed rate and a constat Q. 
%A variance in Feed does have an effect on concentration.
%Increasing the feed rate amt increses the final equilibrium concentration
%and visaversa. This is due to the differntial equation equilibrating when
%the derivative = 0. This happens when F.*(CA0-CA))./dense ==
%(Vr.*k0.*CA.*exp(-ER./T). The second varies CA and T(F), while the first
%part varies CA and F. When F is large, the first part concentration
%difference is amplified, requiring a larger CA in the second part to
%cancel out and equal 0. When F is small, the first part is scaled up less
%and the second part requires a smaller CA to cancel out and equal 0.
figure
hold on
plot(timeval,TvalVary(1,:),'r'); % +50% Feed
plot(timeval,TvalVary(2,:),'b'); % -50% Feed
plot(timeval,TvalVary(3,:),'r.'); % +50% tH
plot(timeval,TvalVary(4,:),'b.'); % -50% tH
plot(timeval,TvalVary(5,:),'r--'); % +50% tTs
plot(timeval,TvalVary(6,:),'b--'); % -50% tTs
plot(timeval,Tval,'k') % no variance
legend('plus feed','minus feed','plus tH','minus tH','plus tTs','minus tTs','no variance');
ylim([341 358]);
title('T temperature in reactor')
ylabel('Temperature (K)')
xlabel('Time (s)')
%The rate of change of temperature relies on the Feed rate, concentration
%and Q, but it was determined earlier that Q was constant throughout this
%process so it only relies on the varying F and CA. 
figure
hold on
plot(timeval,TsvalVary(1,:),'r'); % +50% Feed
plot(timeval,TsvalVary(2,:),'b'); % -50% Feed
plot(timeval,TsvalVary(3,:),'r.'); % +50% tH
plot(timeval,TsvalVary(4,:),'b.'); % -50% tH
plot(timeval,TsvalVary(5,:),'r--'); % +50% tTs
plot(timeval,TsvalVary(6,:),'b--'); % -50% tTs
plot(timeval,Tsval,'k') % no variance
legend('plus feed','minus feed','plus tH','minus tH','plus tTs','minus tTs','no variance');
ylim([341 358]);
title('Ts temperature on sensor')
ylabel('Temperature (K)')
xlabel('Time (s)')


% 7 figures holding the different variances constant across different
% differential components.
ylabels1={'Q','CA','T','Ts'};
figure
novarytbl=[Qval',CAval',Tval',Tsval'];
stackedplot(timeval,novarytbl,'Title','No Variance','DisplayLabels',ylabels1);

figure
plusFvarytbl=[QvalVary(1,:)',CAvalVary(1,:)',TvalVary(1,:)',TsvalVary(1,:)'];
stackedplot(timeval,plusFvarytbl,'Title','50% more Feed','DisplayLabels',ylabels1);

figure
minusFvarytbl=[QvalVary(2,:)',CAvalVary(2,:)',TvalVary(2,:)',TsvalVary(2,:)'];
stackedplot(timeval,minusFvarytbl,'Title','50% less Feed','DisplayLabels',ylabels1);

figure
plustHvarytbl=[QvalVary(3,:)',CAvalVary(3,:)',TvalVary(3,:)',TsvalVary(3,:)'];
stackedplot(timeval,plustHvarytbl,'Title','50% more tH','DisplayLabels',ylabels1);

figure
minustHvarytbl=[QvalVary(4,:)',CAvalVary(4,:)',TvalVary(4,:)',TsvalVary(4,:)'];
stackedplot(timeval,minustHvarytbl,'Title','50% less tH','DisplayLabels',ylabels1);

figure
plustTsvarytbl=[QvalVary(5,:)',CAvalVary(5,:)',TvalVary(5,:)',TsvalVary(5,:)'];
stackedplot(timeval,plustTsvarytbl,'Title','50% more tTs','DisplayLabels',ylabels1);

figure
minustTsvarytbl=[QvalVary(6,:)',CAvalVary(6,:)',TvalVary(6,:)',TsvalVary(6,:)'];
stackedplot(timeval,minustTsvarytbl,'Title','50% less tTs','DisplayLabels',ylabels1);
%
%
Qval2(1)=700000;
CAval2(1)=.25;
Tval2(1)=350;
Tsval2(1)=350;
timeval2(1)=0;
a2=timeval2(1);
b2=200;
h2=.25;
N2=(b2-a2)./h2;

for i=1:N2
    timeval2(i+1)=timeval2(i)+h2;
    if timeval2(i)<100
        K1=dQ(Qval2(i),Qspec,tH);
        K2=dQ(Qval2(i)+K1.*h2./2,Qspec,tH);
        K3=dQ(Qval2(i)+K2.*h2./2,Qspec,tH);
        K4=dQ(Qval2(i)+K3.*h2,Qspec,tH);
        Qval2(i+1)=Qval2(i) + (K1+2.*K2+2.*K3+K4).*h2./6;
    else
        K1=dQ(Qval2(i),Qspec2,tH);
        K2=dQ(Qval2(i)+K1.*h2./2,Qspec2,tH);
        K3=dQ(Qval2(i)+K2.*h2./2,Qspec2,tH);
        K4=dQ(Qval2(i)+K3.*h2,Qspec2,tH);
        Qval2(i+1)=Qval2(i) + (K1+2.*K2+2.*K3+K4).*h2./6;
    end
    K1=dCA(CAval2(i),Tval2(i),F);
    K2=dCA(CAval2(i)+K1.*h2./2,Tval2(i),F);
    K3=dCA(CAval2(i)+K2.*h2./2,Tval2(i),F);
    K4=dCA(CAval2(i)+K3.*h2,Tval2(i),F);
    CAval2(i+1)=CAval2(i) + (K1+2.*K2+2.*K3+K4).*h2./6;
    K1=dT(Tval2(i),CAval2(i),Qval2(i),F);
    K2=dT(Tval2(i)+K1.*h2./2,CAval2(i),Qval2(i),F);
    K3=dT(Tval2(i)+K2.*h2./2,CAval2(i),Qval2(i),F);
    K4=dT(Tval2(i)+K3.*h2,CAval2(i),Qval2(i),F);
    Tval2(i+1)=Tval2(i) + (K1+2.*K2+2.*K3+K4).*h2./6;
    K1=dTs(Tsval2(i),Tval2(i),tTs);
    K2=dTs(Tsval2(i)+K1.*h2./2,Tval2(i),tTs);
    K3=dTs(Tsval2(i)+K2.*h2./2,Tval2(i),tTs);
    K4=dTs(Tsval2(i)+K3.*h2,Tval2(i),tTs);
    Tsval2(i+1)=Tsval2(i) + (K1+2.*K2+2.*K3+K4).*h2./6;
end

%comparison of Qspec changing at 100s and standard graph. should overlap
figure
hold on
plot(timeval,Qval,'k')
plot(timeval2,Qval2,'r')
title('Q with Qspec change at 100s')

figure
hold on
plot(timeval,CAval,'k')
plot(timeval2,CAval2,'r')
title('CA with Qspec change at 100s')

figure
hold on
plot(timeval,Tval,'k')
plot(timeval2,Tval2,'r')
title('T with Qspec change at 100s')

figure
hold on
plot(timeval,Tsval,'k')
plot(timeval2,Tsval2,'r')
title('Ts with Qspec change at 100s')

%extra graph, shows how Ts trails T
figure
hold on
plot(timeval,Tval,'k')
plot(timeval,Tsval,'r')
legend('T','tTs')
title('extra graph: temperature sensor delay')
