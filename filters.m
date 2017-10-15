function filters

N=101;

px=linspace(0,1,N);



lambda=.5;
close all


for i=1:N
    avar=i/N;
    if (0<= avar && avar<= lambda )
    ceniceros(i)=1;
    elseif (lambda<avar && avar<=1) 
    ax=1+((-avar+lambda)/(1-lambda));
    ceniceros(i)=(35-84*ax+70*ax^2-20*ax^3)*ax^4;
    end
end


for i=1:N
    avar=i/N;
    if (0<= avar && avar<= lambda )
    LM(i)=1;
    elseif (lambda<avar && avar<=1) 
             aux3=(1-((avar-lambda)/(1-lambda))^2);
             aquotient=exp(-1/aux3);
             
             LM(i)=exp(1)*aquotient;   
    end
end


lambda=.5
             
    for i=1:N
    avar=i/N;
    if (0<= avar && avar<= lambda )
    M(i)=1;
    elseif (lambda<avar && avar<=1) 
             aux3=1/(1-((avar-lambda)/(1-lambda))^1)^4
             aquotient=exp(-aux3)
             auxrho=exp(1)*aquotient
             
             M(i)=exp(1)*aquotient;   
    end
end
figure(1)
plot(px,ceniceros,'r');
hold on
plot(px,LM,'b')
hold on
plot(px,M,'m*')

