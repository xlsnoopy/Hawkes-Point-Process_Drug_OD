%% 50 experiments
% this section of the program runs 50 repeated simulations on our proposed
% model, and compares the parameters at convergence with their true values

iters = 50;
K0list = zeros(iters,4);
wlist = zeros(iters,4);
mulist = zeros(iters,4);
siglist = zeros(iters,4);
for rep = 1:iters
Geodata = [];

% simulate 4 point processes with 4 different settings of parameters
[times1, x1, y1, p] = Hawkes_Simulation(0.03, 0.9, 0.1, 2500, 0.01, ....
    [1,2,3,4]);
topic = ones(size(times1,1),1);
background = [ones(p,1); zeros(size(times1,1)-p,1)];
Geodata = [Geodata; x1,y1,times1,topic, background];
[times2, x2, y2,p] = Hawkes_Simulation(0.01,0.8, 0.5, 2500, 0.001, ...
    [4,3,2,1]);
topic = ones(size(times2,1),1)*2;
background = [ones(p,1); zeros(size(times2,1)-p,1)];
Geodata = [Geodata; [x2, y2, times2, topic, background]];
[times3, x3, y3,p] = Hawkes_Simulation(0.02,0.6, 1, 2500, 0.02, ...
    [4,4,1,1]);
topic = ones(size(times3,1),1)*3;
background = [ones(p,1); zeros(size(times3,1)-p,1)];
Geodata = [Geodata; [x3, y3, times3, topic, background]];
[times4, x4, y4,p] = Hawkes_Simulation(0.05,0.75, 0.3, 2500, 0.003, ...
    [1,4,1,4]);
topic = ones(size(times4,1),1)*4;
background = [ones(p,1); zeros(size(times4,1)-p,1)];
Geodata = [Geodata; [x4, y4, times4, topic,background]];

% Normalize latitude, longitude and time
Geodata(:,1) = (Geodata(:,1)-min(Geodata(:,1)))/(max(Geodata(:,1))-min(Geodata(:,1))); % Lat
Geodata(:,2) = (Geodata(:,2)-min(Geodata(:,2)))/(max(Geodata(:,2))-min(Geodata(:,2))); % Long
Geodata(:,3) = (Geodata(:,3)-min(Geodata(:,3)));%/(max(Data(:,3))-min(Data(:,3))); % Time

% Randomly assign 30% of data to labelled (EMS), the rest unlabelled are called OP
perc = 0.3;
N = size(Geodata,1);
label = ones(N,1)*2;
X1 = datasample(1:N,floor(perc*N),'Replace', false);
label(X1) = 1; % label EMS as 1, OP as 2
N_ems = floor(perc*N);
N_op = N-N_ems;

Data = [Geodata(:,1:end-1), label];

Op = Data(Data(:,5) == 2,:);
probs_op=zeros(N_op,4);
% Initialize probs_op with ith column = 1 where its topic = i
for i=1:N_op
    probs_op(i,Op(i,4))=1;
end
Opioid = [Op, probs_op];
EMS = Data(Data(:,5) == 1,:);
probs_ems=ones(N_ems,4);
% Initialize probs_ems with .25 for each column
for i = 1:N_ems
    probs_ems(i,:) = ones(4,1)*.25;
end
EMS = [EMS, probs_ems];
Data = [EMS; Opioid];
Data = sortrows(Data, 3);

% Initialize parameters
K0 = ones(4,1)*.5; w = ones(4,1)*.1; mu = ones(4,1)*50; sigma = ones(4,1)*.01;
% Perform EM update - converges within 50 iterations
[K0, w, mu, sigma, u, v, new_P] = EMH(Data, K0, w, mu, sigma, 50, 0.1, 100);

% disp([K0, w, mu, sigma]);
K0list(rep,:) = K0;
wlist(rep,:) = w;
mulist(rep,:) = mu;
siglist(rep,:) = sigma;
end 

% Plot the parameters  
subplot(2,2,1)
histogram(K0list(:,1),'FaceColor',[0.4,0.6,0.7])
hline1 = line([.9 .9], get(gca, 'ylim'));
hline2 = line([mean(K0list(:,1)),mean(K0list(:,1))], get(gca,'ylim'));
hline1.Color = 'r'; hlin2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
title('Group 1');
subplot(2,2,2)
histogram(K0list(:,2),'FaceColor',[0.4,0.6,0.7])
hline1 = line([.8 .8], get(gca, 'ylim'));
hline2 = line([mean(K0list(:,2)),mean(K0list(:,2))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
title('Group 2');
subplot(2,2,3)
histogram(K0list(:,3),'FaceColor',[0.4,0.6,0.7])
hline1 = line([.6 .6], get(gca, 'ylim'));
hline2 = line([mean(K0list(:,3)),mean(K0list(:,3))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
title('Group 3');
subplot(2,2,4)
histogram(K0list(:,4),'FaceColor',[0.4,0.6,0.7])
hline1 = line([.75 .75], get(gca, 'ylim'));
hline2 = line([mean(K0list(:,4)),mean(K0list(:,4))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
title('Group 4');
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
sgtitle('Histogram of K0 on 50 Realizations')
set(gcf, 'PaperPosition', [0 0 3.5 3.5]); 
set(gcf, 'PaperSize', [3.5 3.5]); 
exportgraphics(gcf,'K0.pdf','Resolution',300);


subplot(2,2,1)
histogram(wlist(:,1), 'FaceColor',[0.4,0.6,0.7])
hline1 = line([.1 .1], get(gca, 'ylim'));
hline2 = line([mean(wlist(:,1)),mean(wlist(:,1))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
title('Group 1');
subplot(2,2,2)
histogram(wlist(:,2), 'FaceColor',[0.4,0.6,0.7])
hline1 = line([.5 .5], get(gca, 'ylim'));
hline2 = line([mean(wlist(:,2)),mean(wlist(:,2))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
title('Group 2');
subplot(2,2,3)
histogram(wlist(:,3), 'FaceColor',[0.4,0.6,0.7])
hline1 = line([1 1], get(gca, 'ylim'));
hline2 = line([mean(wlist(:,3)),mean(wlist(:,3))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
title('Group 3');
subplot(2,2,4)
histogram(wlist(:,4), 'FaceColor',[0.4,0.6,0.7])
hline1 = line([.3 .3], get(gca, 'ylim'));
hline2 = line([mean(wlist(:,4)),mean(wlist(:,4))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);title('Group 4');
sgtitle('Histogram of w on 50 Realizations')
set(gcf, 'PaperPosition', [0 0 3.5 3.5]); 
set(gcf, 'PaperSize', [3.5 3.5]); 
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf,'w.pdf','Resolution',300);

subplot(2,2,1)
h1 = histogram(siglist(:,1), 'FaceColor',[0.4,0.6,0.7]);
hline1 = line([.01 .01], get(gca, 'ylim'));
hline2 = line([mean(siglist(:,1)),mean(siglist(:,1))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
ax = gca;
ax.XRuler.Exponent = 0;
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
title('Group 1');
subplot(2,2,2)
histogram(siglist(:,2), 'FaceColor',[0.4,0.6,0.7])
hline1 = line([.001 .001], get(gca, 'ylim'));
hline2 = line([mean(siglist(:,2)),mean(siglist(:,2))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
ax = gca;
ax.XRuler.Exponent = 0;
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
title('Group 2');
subplot(2,2,3)
histogram(siglist(:,3), 'FaceColor',[0.4,0.6,0.7])
hline1 = line([.02 .02], get(gca, 'ylim'));
hline2 = line([mean(siglist(:,3)),mean(siglist(:,3))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
ax = gca;
ax.XRuler.Exponent = 0;
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
title('Group 3');
subplot(2,2,4)
histogram(siglist(:,4), 'FaceColor',[0.4,0.6,0.7])
hline1 = line([.003 .003], get(gca, 'ylim'));
hline2 = line([mean(siglist(:,4)),mean(siglist(:,4))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
ax = gca;
ax.XRuler.Exponent = 0;
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);title('Group 4');
sgtitle('Histogram of sigma on 50 Realizations')
set(gcf, 'PaperPosition', [0 0 3.5 3.5]); 
set(gcf, 'PaperSize', [3.5 3.5]); 
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf,'sigma.pdf','Resolution',300);

subplot(2,2,1)
histogram(mulist(:,1), 'FaceColor',[0.4,0.6,0.7])
hline1 = line([67 67], get(gca, 'ylim'));
hline2 = line([mean(mulist(:,1)),mean(mulist(:,1))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
title('Group 1');
subplot(2,2,2)
histogram(mulist(:,2), 'FaceColor',[0.4,0.6,0.7])
hline1 = line([28 28], get(gca, 'ylim'));
hline2 = line([mean(mulist(:,2)),mean(mulist(:,2))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
title('Group 2');
subplot(2,2,3)
histogram(mulist(:,3), 'FaceColor',[0.4,0.6,0.7])
hline1 = line([55 55], get(gca, 'ylim'));
hline2 = line([mean(mulist(:,3)),mean(mulist(:,3))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);
title('Group 3');
subplot(2,2,4)
histogram(mulist(:,4), 'FaceColor',[0.4,0.6,0.7])
hline1 = line([132 132], get(gca, 'ylim'));
hline2 = line([mean(mulist(:,4)),mean(mulist(:,4))], get(gca,'ylim'));
hline1.Color = 'r'; hline2.Color = 'b';
hline1.LineStyle = '-.'; hline2.LineStyle = '-';
set(hline1,'LineWidth',1.2);
set(hline2,'LineWidth', 1.2);title('Group 4');
sgtitle('Histogram of mu on 50 Realizations')
set(gcf, 'PaperPosition', [0 0 3.5 3.5]); 
set(gcf, 'PaperSize', [3.5 3.5]); 
set(gcf, 'PaperPositionMode', 'auto');
exportgraphics(gcf,'mu.pdf','Resolution',300);


%% Changing proportion of labeled and unlabeled data
% this section of program computes the log likelihood using our proposed
% model, and compares it to the baseline models
% with increasing size of unlabelled data and decreasing size of labelled
% data, our model has a higher log likelihood

% Simulate point processes
Geodata = [];
[times1, x1, y1, p] = Hawkes_Simulation(0.03, 0.9, 0.1, 2500, 0.01, ....
    [1,2,3,4]);
topic = ones(size(times1,1),1);
background = [ones(p,1); zeros(size(times1,1)-p,1)];
Geodata = [Geodata; x1,y1,times1,topic, background];
[times2, x2, y2,p] = Hawkes_Simulation(0.01,0.8, 0.5, 2500, 0.001, ...
    [4,3,2,1]);
topic = ones(size(times2,1),1)*2;
background = [ones(p,1); zeros(size(times2,1)-p,1)];
Geodata = [Geodata; [x2, y2, times2, topic, background]];
[times3, x3, y3,p] = Hawkes_Simulation(0.02,0.6, 1, 2500, 0.02, ...
    [4,4,1,1]);
topic = ones(size(times3,1),1)*3;
background = [ones(p,1); zeros(size(times3,1)-p,1)];
Geodata = [Geodata; [x3, y3, times3, topic, background]];
[times4, x4, y4,p] = Hawkes_Simulation(0.05,0.75, 0.3, 2500, 0.003, ...
    [1,4,1,4]);
topic = ones(size(times4,1),1)*4;
background = [ones(p,1); zeros(size(times4,1)-p,1)];
Geodata = [Geodata; [x4, y4, times4, topic,background]];

Geodata(:,1) = (Geodata(:,1)-min(Geodata(:,1)))/(max(Geodata(:,1))-min(Geodata(:,1))); % Lat
Geodata(:,2) = (Geodata(:,2)-min(Geodata(:,2)))/(max(Geodata(:,2))-min(Geodata(:,2))); % Long
Geodata(:,3) = (Geodata(:,3)-min(Geodata(:,3)));%/(max(Data(:,3))-min(Data(:,3))); % Time

iters = 1:9;
ll1 = zeros(length(iters),2); % all data
ll2 = zeros(length(iters),1); % unlabeled only
ll3 = zeros(length(iters),1); % labeled only

% compute log likelihood of each dataset, with changing proportions of
% unlabelled and labelled data
for rep = 1:length(iters)
perc = iters(rep);
N = size(Geodata,1);
label = ones(N,1)*2;
% randomly draw 10% to 90% of data for EMS, rest are OP
X1 = datasample(1:N,floor(perc*N/10),'Replace', false); 
label(X1) = 1;
N_ems = floor(perc*N/10);
N_op = N-N_ems;
Data = [Geodata(:,1:end-1), label];

Op = Data(Data(:,5) == 2,:);
probs_op=zeros(N_op,4);
for i=1:N_op
    probs_op(i,Op(i,4))=1;
end
Opioid = [Op, probs_op];
EMS = Data(Data(:,5) == 1,:);
probs_ems=ones(N_ems,4);
p_op = ones(4,1)*.25;
for i = 1:N_ems
    probs_ems(i,:) = p_op;
end
EMS = [EMS, probs_ems];
Data = [EMS; Opioid];
Data = sortrows(Data, 3);

% proposed model
K0 = ones(4,1)*.5; w = ones(4,1)*.1; mu = ones(4,1)*100; sigma = ones(4,1)*.01;
[K0, w, mu, sigma, u, v, new_P] = EMH(Data, K0, w, mu, sigma, 50, 0.1, 100);
N = N_ems + N_op;

% verify that converged probabilities recover ground truth
group = zeros(N,1);
for i = 1:N
    if Data(i,5)==1
        % draw a group number for each event from OP based on its probabilities
        group(i) = sum(rand>=cumsum([0, sum(new_P{1}(i,:)), sum(new_P{2}(i,:)), sum(new_P{3}(i,:)), sum(new_P{4}(i,:))]));
    else
        % assign its group number for each event from EMS
        group(i) = Data(i,4);
    end
end

% calculate log likelihood
ll_op = zeros(4,1);
ll_ems = zeros(4,1);
for i = 1:N
    k = group(i);
    if Data(i,5) == 1 % unlabeled
        lambda = mu(k)* u{k}(i) * v{k}(i); 
        for j=1:(i-1)
            lambda =  lambda + K0(k)*w(k)*exp(-w(k)*(Data(i,3)-Data(j,3)))*...
            exp((-(Data(i,1)-Data(j,1))^2-(Data(i,2)-Data(j,2))^2)/(2*sigma(k)^2))/(2*pi*sigma(k)^2);          
        end
        ll_ems(k)= ll_ems(k) + log(lambda);

    else %labeled
        lambda = mu(k)* u{k}(i) * v{k}(i); 
        for j=1:(i-1)
            lambda =  lambda + K0(k)*w(k)*exp(-w(k)*(Data(i,3)-Data(j,3)))*...
            exp((-(Data(i,1)-Data(j,1))^2-(Data(i,2)-Data(j,2))^2)/(2*sigma(k)^2))/(2*pi*sigma(k)^2);          
        end
        ll_op(k)= ll_op(k) + log(lambda) ;

    end
end 

for k = 1:4
    ll_ems(k) = ll_ems(k) - K0(k)*sum(sum(new_P{k}(Data(:,5)==1,:))) - sum(diag(new_P{k}(Data(:,5)==1,:)));
    ll_op(k) = ll_op(k) - K0(k)*sum(sum(new_P{k}(Data(:,5)==2,:))) - sum(diag(new_P{k}(Data(:,5)==2,:)));
end



% unlabeled data only, no groups involved
ll2_ems = 0;
EMS = sortrows(EMS, 3);
K0 = .5; w = .1; mu = 10; sigma = .01;
[K0, w, mu, sigma, u,v] = EM(EMS, K0, w, mu, sigma, 50, 0.1, 100);
for i = 1:size(EMS,1)
     lambda = mu* u(i) * v(i); 
     for j=1:(i-1)
         lambda =  lambda + K0*w*exp(-w*(EMS(i,3)-EMS(j,3)))*...
                exp((-(EMS(i,1)-EMS(j,1))^2-(EMS(i,2)-EMS(j,2))^2)/(2*sigma^2))/(2*pi*sigma^2);          
     end 
     ll2_ems = ll2_ems + log(lambda);
end
ll2_ems = ll2_ems - K0*size(EMS,1) - mu;
% labeled data only
ll3_op = zeros(4,1);
Op = sortrows(Opioid, 3);
K0 = ones(4,1)*.5; w = ones(4,1)*.1; mu = ones(4,1)*.1; sigma = ones(4,1)*.01;
for k = 1:4
    datatmp = Op(Op(:,4)==k,:);
    [K0(k), w(k), mu(k), sigma(k), u, v] = ...
        EM(datatmp, K0(k), w(k), mu(k), sigma(k), 50, 0.1, 100);
    for i = 1:size(datatmp,1)
       lambda = mu(k)* u(i) * v(i); 
            for j=1:(i-1)
                lambda =  lambda + K0(k)*w(k)*exp(-w(k)*(datatmp(i,3)-datatmp(j,3)))*...
                exp((-(datatmp(i,1)-datatmp(j,1))^2-(datatmp(i,2)-datatmp(j,2))^2)/(2*sigma(k)^2))/(2*pi*sigma(k)^2);          
            end
        ll3_op(k) = ll3_op(k) + log(lambda);
    end
    ll3_op(k) = ll3_op(k) - K0(k)*size(datatmp,1) - mu(k);
end

ll1(rep,:) = [sum(ll_op), sum(ll_ems)];
ll2(rep,:) = sum(ll2_ems);
ll3(rep,:) = sum(ll3_op);
end

scatter(10:10:90, ll1(:,2),'o');
hold on
scatter(10:10:90, ll2,'+');
legend('ll(A) proposed model','ll(A) baseline model','Location','northwest');
ylabel('Log-likelihood of the model');
xlabel('% of size of A over all events');
hold off
set(gcf, 'PaperPosition', [0 0 3.5 3.5]); 
set(gcf, 'PaperSize', [3.5 3.5]); 
saveas(gcf, 'llA', 'pdf')

scatter(10:10:90, ll1(:,1), 'o');
hold on
scatter(10:10:90, ll3, '+');
ylabel('Log-likelihood of the model');
xlabel('% of size of A over all events');
legend('ll(B) proposed model','ll(B) baseline model');
hold off
set(gcf, 'PaperPosition', [0 0 3.5 3.5]); 
set(gcf, 'PaperSize', [3.5 3.5]); 
saveas(gcf, 'llB', 'pdf')

