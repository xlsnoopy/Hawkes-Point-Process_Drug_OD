function [u]=variablebw2D(x,y,bw2d,p)
N2=size(x);
N=max(N2);
bw=zeros(1,N);
d=zeros(1,N);
u=zeros(size(x));

for i=1:N
%     for j=1:N
%         d(j)=(x(i)-x(j))^2+(y(i)-y(j))^2;
%         d(j)=d(j)^.5;
%     end
%     ds=sort(d);
%     bw(i)=ds(k);
    bw(i) = bw2d;
end

KN=sum(p);
parfor i=1:N    
   for j=[1:(i-1),(i+1):N]
       u(i)=u(i)+p(j)*(1/(2*pi*bw(j)^2))*exp((-(x(i)-x(j))^2-(y(i)-y(j))^2)/(2*bw(j)^2))/(KN);
   end
end

end