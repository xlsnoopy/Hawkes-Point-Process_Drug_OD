function [K0, w, mu, sigma, u, v, new_P] = EMH(Data, K0, w, mu, sigma, niter, bw2d, bw1d)
% this program computes the model parameters "K0", "w", "mu", "sigma" using
% EM iterations
% it also computes the spatial and temporal background rates, called "u"
% and "v", using kernel density estimation with bandwidths, 
% called "bw2d" and "bw1d"

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
        u{i} = variablebw2D(Data(:,1), Data(:,2), bw2d, diag(new_P{i})); 
        v{i} = variablebw1D(Data(:,3), bw1d, diag(new_P{i})); 
    end
    
    for i = 1:4
        % E-Step: calculate probabilities
        new_P{i} = updatep_m1(Data(:,3), Data(:,1), Data(:,2), new_P{i}, ...
            K0(i),w(i),mu(i),sigma(i),u{i}, v{i});
    end
   
    % for each labelled data, reset their probabilities to 1 if its from
    % i-th group, 0 otherwise
    % compute sum of probabilities of each event of all groups
    % and normalize to be summed up to 1
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
        [K0(i),w(i),mu(i), sigma(i)]=updatepar(Data, new_P{i});
    end
end
end