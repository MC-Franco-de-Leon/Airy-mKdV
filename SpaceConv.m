function Ma=SpaceConv(Th,Th2,Th4,Th8,Th16,Th32,t)
[a,b]=size(Th);

h=1/a;
n=a;
format shortE
M=double(zeros(9,b-1));

for j=1:b-1

for i=1:n
    error1(i)=Th(i,j+1)-Th2(2*i-1,j+1);
    error2(i)=Th2(2*i-1,j+1)-Th4(4*i-3,j+1);
    error3(i)=Th4(4*i-3,j+1)-Th8(8*i-7,j+1);
    error4(i)=Th8(8*i-7,j+1)-Th16(16*i-15,j+1);
    error5(i)=Th16(16*i-15,j+1)-Th32(32*i-31,j+1);
end

norm(error1,inf);
norm(error2,inf);
norm(error3,inf);
norm(error4,inf);
norm(error5,inf);

l21=0;
l22=0;
l23=0;
l24=0;
l25=0;

for i=1:n
    l21=l21+h*(error1(i))^2;
    l22=l22+h*(error2(i))^2;
    l23=l23+h*(error3(i))^2;
    l24=l24+h*(error4(i))^2;
    l25=l25+h*(error5(i))^2;
end
   


M(1,j)=sqrt(l21);
M(2,j)=sqrt(l22);
M(3,j)=sqrt(l23);
M(4,j)=sqrt(l24);
M(5,j)=sqrt(l25);


M(6,j)= log2(M(1,j)/M(2,j));
M(7,j)= log2(M(2,j)/M(3,j));
M(8,j)= log2(M(3,j)/M(4,j));
M(9,j)= log2(M(4,j)/M(5,j));

end

Ma(:,1)=t(2:end);
Ma(:,2)=M(1,:);
Ma(:,3)=M(2,:);
Ma(:,4)=M(3,:);
Ma(:,5)=M(4,:);
Ma(:,6)=M(5,:);

Ma(:,7)=M(6,:);
Ma(:,8)=M(7,:);
Ma(:,9)=M(8,:);
Ma(:,10)=M(9,:);
