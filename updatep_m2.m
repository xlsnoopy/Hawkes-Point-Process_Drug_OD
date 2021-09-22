%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to estimate branching probabilities

function [p,lam]=updatep_m2(t,x,y,p,K0,w,mu,sigma,u,v)

N=size(t,1);

for i=1:N
    for j=1:(i-1)
        % probability i triggered by j is proportional to triggering
        % kernel evaluated at inter-point times and distances
        p(i,j)=K0*w*exp(-w*(t(i)-t(j)))*...
            exp((-(x(i)-x(j))^2-(y(i)-y(j))^2)/(2*sigma^2))/(2*pi*sigma^2);
          p(i,i) =  (mu)* u(i) * v(i);
    end
    % save intensity at each event for analysis 
    lam(i)=sum(p(i,1:i));
    %normalize probabilities to sum to 1
    p(i,1:i)=p(i,1:i)/sum(p(i,1:i));
end
end
