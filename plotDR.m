function plotDR
close all;
clear all;

R0=1;
m=3;

d0=[.05 .1 .15 .2 .25 ];

c=(m^3-m)/R0^3;

D=d0;
R=R0*[1 1 1 1 1];
tp;
slp;
Nframes=length(t);

clear vars Delta
clear vars r0
Deltap;
r0p;
Md(:,5)=Delta;
Mr(:,5)=r0;
cd ..

cd D4
clear vars Delta
clear vars r0
Deltap;
r0p;
Md(:,4)=Delta;
Mr(:,4)=r0;
cd ..

cd D3
clear vars Delta
clear vars r0
Deltap;
r0p;
Md(:,3)=Delta;
Mr(:,3)=r0;
cd ..

cd D2
clear vars Delta
clear vars r0
Deltap;
r0p;
Md(:,2)=Delta;
Mr(:,2)=r0;
cd ..

cd D1
clear vars Delta
clear vars r0
Deltap;
r0p;
Md(:,1)=Delta;
Mr(:,1)=r0;
cd ..
cd D5








for j=1:5 %different delta0
    for i=1:Nframes % for each frame in time
        eD(i,j)=abs(D(j)-Md(i,j));
        eR(i,j)=abs(R(j)-Mr(i,j));
    end
end

cmap = hsv(Nframes);  %# Creates a Nframes-by-3 set of colors from the HSV colormap


figure(1)
plot(log(d0),log(d0.^2),'k*-');
for i=1:Nframes
    poly=polyfit(log(D),log(eD(i,:)),1);
    Slopd(i)=poly(1);
    hold on
    plot(log(d0),log(eD(i,:).^2),'--s','Color',cmap(i,:));
%    title(sprintf('Slope = %d ',Slopd(i)));
    pause(.20)    
    axis equal;
%    hold off
end
    legend('x^2','t=.01024','t=.02048','t=.03072','t=.04096','t=.05120','t=.06144','t=.07168','t=.08192','t=.09216','t=.1024')
    xlabel('Initial perturbation \delta_0')
    ylabel('|\delta_L-\delta_N|^2')
Md=mean(Slopd)
hold off

    
%figure(2)
plot(log(d0),log(d0.^2),'k*-')
for i=1:Nframes
    poly=polyfit(log(D),log(eR(i,:)),1);
    Slopr(i)=poly(1);
%    hold on
    plot(log(d0),log(eR(i,:).^2),'--s','Color',cmap(i,:));
%    title(sprintf('Slope = %d ',Slopr(i)));
    pause(.20)
    axis equal;
end
    %xlabel(sprintf('Time = %d', t(i)))
    legend('x^2','t=.01024','t=.02048','t=.03072','t=.04096','t=.05120','t=.06144','t=.07168','t=.08192','t=.09216','t=.1024')

    xlabel('Initial perturbation \delta_0')
    ylabel('|R_L-R_N|^2')
    Mr=mean(Slopr)
    hold off



end