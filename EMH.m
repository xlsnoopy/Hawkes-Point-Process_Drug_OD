function [K0, w, mu1, mu2, sigma, u, v, new_P] = EMH(Data, K0, w, mu1, mu2, sigma, niter, bw2d, bw1d)
N = size(Data,1);
new_P = [];
for i = 1:4
    A = zeros(N,N);
    A(logical(eye(N))) = Data(:,5+i);
    A(1,1) = 0.25;
    new_P{i} = A;
end
%% EM-H iterations
for iter = 1:niter
    % estimate group for unknown data
    
    % Background rates
    u = []; v = []; v_h = [];  v_d = []; v_m = []; v_y = [];
    for i = 1:4        
        u{i} = variablebw2D(Data(:,1), Data(:,2), bw2d, diag(new_P{i})); %0.1
        v{i} = variablebw1D(Data(:,3), bw1d, diag(new_P{i})); %100
%         for j = 1:N
%             v_h{i}(j) = sum(bgP(Data(:,14) == Data(j,14)))/sum(bgP);
%             v_d{i}(j) = sum(bgP(Data(:,15) == Data(j,15)))/sum(bgP);
%             v_m{i}(j) = sum(bgP(Data(:,16) == Data(j,16)))/sum(bgP);
%             v_y{i}(j) = sum(bgP(Data(:,17) == Data(j,17)))/sum(bgP);
%         end
    end
    
    for i = 1:4
        % E-Step: calculate probabilities
        new_P{i} = updatep_m1(Data(:,3), Data(:,1), Data(:,2), new_P{i}, ...
            K0(i),w(i),mu1(i), mu2(i),sigma(i),u{i}, v{i});% v_h{i}, v_d{i}, v_m{i}, v_y{i});
    end
   
    for j = 1:size(new_P{i},1)
         for i = 1:4
            if Data(j,5) == 2
                if Data(j,4) ~= i
                    new_P{i}(j,1:j) = zeros(size(new_P{i}(j,1:j)));
                end
            end
         end
         temp_sum = 0; 
         for i = 1:4
            temp_sum=temp_sum+sum(new_P{i}(j,1:j));
         end
         for i=1:4
            new_P{i}(j,1:j) = new_P{i}(j,1:j)/temp_sum;
        end
    end

    for i =1:4
        % M-Step: update params
        [K0(i),w(i),mu1(i), mu2(i), sigma(i)]=updatepar(Data, new_P{i});
    end
    
    
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
disp(iter);
disp(mu1+mu2);
end
end