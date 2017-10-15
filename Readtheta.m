function [Th]=Readtheta(theta)
for i=1:20
    Th(:,i)=theta(:,i*10);
end

clear vars theta
%clear vars t
%clear vars H1
%clear vars H2
%clear vars H3

