%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to estimate branching probabilities

function p=updatep_m1(t,x,y,p,K0,w,mu,sigma,u,v)

N=size(t,1);
priors=zeros(size(p,1),1);
for i=1:size(priors,1)
    priors(i)=sum(p(i,:));
end
parfor i=1:N
    prow = p(i,:);
    for j=1:(i-1)
        % probability i triggered by j is proportional to triggering
        % kernel evaluated at inter-point times and distances
        prow(j)=K0*w*exp(-w*(t(i)-t(j)))*...
            exp((-(x(i)-x(j))^2-(y(i)-y(j))^2)/(2*sigma^2))/(2*pi*sigma^2)*priors(j);
          prow(i) =  (mu)* u(i) * v(i);

    end
    p(i,:) = prow;
end

end
