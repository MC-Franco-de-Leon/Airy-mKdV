function ReadD
close all
thetap
H1p
H2p
H3p
TCtH6=theta;
HCt6(:,1)=H1;
HCt6(:,2)=H2;
HCt6(:,3)=H3;
cd ..
cd dir5
clear vars theta
clear vars H1
clear vasr H2
clear vasr H3
thetap
H1p
H2p
H3p
TCtH5=theta;
HCt5(:,1)=H1;
HCt5(:,2)=H2;
HCt5(:,3)=H3;
cd ..
cd dir4
clear vars theta
clear vars H1
clear vasr H2
clear vasr H3
thetap
H1p
H2p
H3p
TCtH4=theta;
HCt4(:,1)=H1;
HCt4(:,2)=H2;
HCt4(:,3)=H3;
cd ..
cd dir3
clear vars theta
clear vars H1
clear vasr H2
clear vasr H3
thetap
H1p
H2p
H3p
TCtH3=theta;
HCt3(:,1)=H1;
HCt3(:,2)=H2;
HCt3(:,3)=H3;
cd ..
cd dir2
clear vars theta
clear vars H1
clear vasr H2
clear vasr H3
thetap
H1p
H2p
H3p
TCtH2=theta;
HCt2(:,1)=H1;
HCt2(:,2)=H2;
HCt2(:,3)=H3;
cd ..
cd dir1
clear vars theta
clear vars H1
clear vasr H2
clear vasr H3
thetap
H1p
H2p
H3p
TCtH1=theta;
HCt1(:,1)=H1;
HCt1(:,2)=H2;
HCt1(:,3)=H3;
tp
clear vars theta
clear vars H1
clear vasr H2
clear vasr H3
cd ..
cd dir6

visualizeHH(t,HCt1,HCt2,HCt3,HCt4,HCt5,HCt6); 
Mt=TimeConv(TCtH1,TCtH2,TCtH3,TCtH4,TCtH5,TCtH6,t);
clear vars aux
aux(:,1)=Mt(:,1);
aux(:,2:5)=Mt(:,7:10);
Ordert=OTable(aux)
Errort=ETable(Mt(:,1:6))

M1=max(abs(HCt1(:,1)));
M2=max(abs(HCt2(:,1)));
M3=max(abs(HCt3(:,1)));
M4=max(abs(HCt4(:,1)));
M5=max(abs(HCt5(:,1)));
M6=max(abs(HCt6(:,1)));
M1=[M1 M2 M3 M4 M5 M6]
M=max(M1)

cd ..
myfile =fopen('Tables.txt', 'w');
fprintf(myfile,'%c',Ordert);
fprintf(myfile,'%c',Errort);
fprintf(myfile,'Max. val. of H1 is: %d ',M);


