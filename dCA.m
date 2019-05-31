function dCA = dCA(CA,T,F)
Vr=100;
CA0=1;
dense=1000;
k0=1.97.*10.^24;
ER=20000;
dCA=(((F.*(CA0-CA))./dense)-(Vr.*k0.*CA.*exp(-ER./T)))./Vr;
end
