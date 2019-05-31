function dT = dT(T,CA,Q,F)
CP=1;
T0=350;
Vr=100;
Hrxn=160000;
k0=1.97.*10.^24;
ER=20000;
dense=1000;
dT=((F.*CP.*(T0-T))-(Vr.*Hrxn.*CA.*k0.*exp(-ER./T))+Q)./(Vr.*dense.*CP);
end
