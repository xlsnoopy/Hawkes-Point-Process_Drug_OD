function [v]=variablebw1D(t,bw1d,p)
N2=size(t);
N=max(N2);
bw=zeros(1,N);
d=zeros(1,N);
v=zeros(size(t));

for i=1:N
%     for j=1:N
%         d(j)=(t(i)-t(j))^2;
%         d(j)=d(j)^.5;
%     end
%     ds=sort(d);
%     bw(i)=ds(k);
    bw(i) = bw1d;
end

KN=sum(p);
parfor i=1:N
   for j= [1:(i-1), (i+1):N]
       v(i)=v(i)+p(j)*(1/(2*pi*bw(j)^2))*exp((-(t(i)-t(j))^2)/(2*bw(j)^2))/(KN);
   end
end

end