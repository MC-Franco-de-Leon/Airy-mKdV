function ReadDS
clear all 
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
cd Pru5
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
cd Pru4
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
cd Pru3
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
cd Pru2
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
cd Pru1
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
cd Pru6
visualizeHH(t,HCt1,HCt2,HCt3,HCt4,HCt5,HCt6); 
Mt=SpaceConv(TCtH1,TCtH2,TCtH3,TCtH4,TCtH5,TCtH6,t);
clear vars aux
aux(:,1)=Mt(:,1);
aux(:,2:5)=Mt(:,7:10);
OTable(aux)
ETable(Mt(:,1:6))
