function plotDRframes
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
deloradp;
Mdr(:,5)=delorad;
Md(:,5)=Delta;
Mr(:,5)=r0;
cd ..

cd D4
clear vars Delta
clear vars r0
clear vars delorad
Deltap;
r0p;
deloradp;
Mdr(:,4)=delorad;
Md(:,4)=Delta;
Mr(:,4)=r0;
cd ..

cd D3
clear vars Delta
clear vars r0
clear vars delorad
Deltap;
r0p;
deloradp;
Mdr(:,3)=delorad;
Md(:,3)=Delta;
Mr(:,3)=r0;
cd ..

cd D2
clear vars Delta
clear vars r0
clear vars delorad
Deltap;
r0p;
deloradp;
Mdr(:,2)=delorad;
Md(:,2)=Delta;
Mr(:,2)=r0;
cd ..

cd D1
clear vars Delta
clear vars r0
clear vars delorad
Deltap;
r0p;
deloradp;
Mdr(:,1)=delorad;
Md(:,1)=Delta;
Mr(:,1)=r0;
cd ..
cd D5








for j=1:5 %different delta0
    for i=1:Nframes % for each frame in time
        eD(i,j)=abs(D(j)-Md(i,j));
        eR(i,j)=abs(R(j)-Mr(i,j));
        edR(i,j)=abs(D(j)/R(j)-Mdr(i,j));

    end
end

cmap = hsv(Nframes);  %# Creates a Nframes-by-3 set of colors from the HSV colormap


auxfig1=figure(1)
auxfig2=figure(2)

for i=1:Nframes
    figure(1)
    poly=polyfit(log(D),log(eD(i,:)),1);
    Slopd(i)=poly(1);
    figure(1)
    plot(d0.^2,d0.^2,'k*-');
    hold on 
    plot(d0.^2,eD(i,:),'--s','Color',cmap(i,:));
    axis equal;
    legend('x',sprintf('Time = %d', t(i)))
    xlabel('Initial perturbation \delta_0^2')
    ylabel('|\delta_L-\delta_N|')
    M1(i)=getframe(auxfig1); % leaving gcf out crops the frame in the movie.

    hold off
    
end

    
Md=mean(Slopd)
%hold off

    
for i=1:Nframes
    figure(2)
    poly=polyfit(log(D),log(eR(i,:)),1);
    Slopr(i)=poly(1);
    plot(d0.^2,d0.^2,'k*-')
    hold on
    plot(d0.^2,eR(i,:),'--s','Color',cmap(i,:));
    axis equal;
    legend('x',sprintf('Time = %d', t(i)))
    xlabel('Initial perturbation \delta_0^2')
    ylabel('|R_L-R_ N|')
    M2(i)=getframe(auxfig2); % leaving gcf out crops the frame in the movie.

    hold off    
    
end

    Mr=mean(Slopr)
    hold off

movie2avi(M1,'MovieDelta.avi');
movie2avi(M2,'MovieRadius.avi');

    
%     
% for i=1:Nframes
% 
%     figure(16)
% 
%     plot(d0.^2,d0.^2,'k*-')
%     hold on
%     plot(d0.^2,edR(i,:),'--s','Color',cmap(i,:));
%     axis equal;
%     legend('x',sprintf('Time = %d', t(i)))
%     xlabel('Initial perturbation \delta_0^2')
%     ylabel('|D_L/R_L-D_N/R_N|')
%     hold off     
% end

end