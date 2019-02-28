%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to estimate branching probabilities

function p=updatep_m1(t,x,y,p,K0,w,mu1, mu2,sigma,u,v)%v_d,v_w,v_m,v_y)

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
          prow(i) =  (mu1+mu2)* u(i) * v(i);%v_d(i) * v_w(i) * v_m(i) * v_y(i);

        % save intensity at each event for analysis 
        % lam(i)=sum(p(i,1:i));
    end
    p(i,:) = prow;
    %normalize probabilities to sum to 1
    %p(i,1:i)=p(i,1:i)/sum(p(i,1:i));
end


% for i=1:N
%     for j=1:(i-1)
%         % probability i triggered by j is proportional to triggering
%         % kernel evaluated at inter-point times and distances
%         p(i,j)=K0*w*exp(-w*(t(i)-t(j)))*...
%             exp((-(x(i)-x(j))^2-(y(i)-y(j))^2)/(2*sigma^2))/(2*pi*sigma^2)*...
%             priors(j);
%           p(i,i) =  mu* u(i) * v(i);%v_d(i) * v_w(i) * v_m(i) * v_y(i);
% 
%         % save intensity at each event for analysis 
%         % lam(i)=sum(p(i,1:i));
%     end
%     %normalize probabilities to sum to 1
%     %p(i,1:i)=p(i,1:i)/sum(p(i,1:i));
% end

end
