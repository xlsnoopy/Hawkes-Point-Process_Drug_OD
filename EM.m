function [K0, w, mu, sigma, u, v, old_P] = EM(Data, K0, w, mu, sigma, niter, bw2d, bw1d)
N = size(Data,1);
old_P = diag(ones(N,1));
%% EM-H iterations
for iter = 1:niter
    % estimate group for unknown data
    
    % Background rates
    u = variablebw2D(Data(:,1), Data(:,2), bw2d, diag(old_P)); %0.1
    v = variablebw1D(Data(:,3), bw1d, diag(old_P)); %100
%         bgP = diag(old_P{i});
%         for j = 1:N
%             v_h{i}(j) = sum(bgP(Data(:,14) == Data(j,14)))/sum(bgP);
%             v_d{i}(j) = sum(bgP(Data(:,15) == Data(j,15)))/sum(bgP);
%             v_m{i}(j) = sum(bgP(Data(:,16) == Data(j,16)))/sum(bgP);
%             v_y{i}(j) = sum(bgP(Data(:,17) == Data(j,17)))/sum(bgP);
%         end
    
        % E-Step: calculate probabilities
     new_P = updatep_m2(Data(:,3), Data(:,1), Data(:,2), old_P, ...
        K0,w,mu, sigma,u, v);% v_h{i}, v_d{i}, v_m{i}, v_y{i});
        


        % M-Step: update params
    [K0,w,mu,sigma]=updatepar(Data, new_P);
%     [K0,w,mu,sigma]
%     u
%     v
    old_P = new_P;
    
    % H-step: resample hour of day for group-known data
%     lam = zeros(1,24);
%     for i = 1:4
%         bgP = diag(old_P{i});
%         for hhat = 1:24
%             density(i,hhat) = sum(bgP(Data(:,14) == hhat))/sum(bgP);
%         end
%     end
%     for j = 1:N
%         if Data(j,5) == 2
%             Data(j,14) = sum(rand >= cumsum([0, density(Data(j,4),:)]));
%         end
%     end
disp([K0, w, mu, sigma]);
end
end





