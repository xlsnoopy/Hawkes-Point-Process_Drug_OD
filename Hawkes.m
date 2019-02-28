%% Data pre-processing
Opioid = readtable('/Users/xlsnoopy/Downloads/EMS_0524/indy_opioid-geo5.csv');
EMS = readtable('/Users/xlsnoopy/Downloads/EMS_0618/EMS_overdose.csv');
% Opioid = readtable('indy_opioid-geo5.csv');
% EMS = readtable('EMS_overdose.csv');

Lat_op = Opioid.Latitude;
Long_op = Opioid.Longitude;
Topic_op = Opioid.NMF4;
Time_op = datetime(Opioid.Date, 'InputFormat', 'MM/dd/yyyy');
t_op = daysact(min(Time_op), Time_op);
N_op = size(t_op, 1);
% % Week_op = zeros(N_op,1);
% % for i = 1:N_op
% %     Week_op(i) = weekday(Opioid.Date(i),'mm/dd/yyyy');
% % end
% Year_op = year(Opioid.Date);
% Month_op = month(Opioid.Date);
% Day_op = zeros(N_op,1);
Label = ones(N_op,1)*2;
Opioid = [Lat_op, Long_op, t_op, Topic_op, Label];
% Opioid = sortrows(Opioid,3);
probs_op=zeros(N_op,4);
for i=1:N_op
    probs_op(i,Opioid(i,4))=1;
    %probs_op(i,:) = 1;
end
Opioid = [Opioid, probs_op];% Day_op, Week_op, Month_op, Year_op];
%Opioid = sortrows(Opioid,3);

Lat_ems = EMS.XCOORD;
Long_ems = EMS.YCOORD;
Time_ems = datetime(EMS.DATE_TIME, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
t_ems = daysact(min(Time_op), Time_ems);
N_ems = size(t_ems,1);
% Week_ems = zeros(N_ems,1);
% for i = 1:N_ems
%     Week_ems(i) = weekday(EMS.DATE_TIME(i),'mm/dd/yyyy');
% end
% Year_ems = year(EMS.DATE_TIME);
% Month_ems = month(EMS.DATE_TIME);
% Day_ems = t_ems-floor(t_ems);
Topic_ems = zeros(N_ems,1);

p_op = zeros(4,1); % proportion of events for each topic
for i = 1:4
    p_op(i) = sum(Opioid(:,4) == i)/ N_op;
end

for i = 1:N_ems
    Topic_ems(i) = sum(rand >= cumsum(p_op))+1;
end
Label = ones(N_ems,1);

EMS = [Lat_ems, Long_ems, t_ems, Topic_ems, Label]; %Day_ems, Week_ems, Month_ems, Year_ems
ind = t_ems >= 0;
EMS = EMS(ind,:); % remove events that occur before Op
N_ems = size(EMS,1);
% N_ems = size(EMS(1:5:end,:), 1);

probs_ems=ones(N_ems,4);
for i = 1:N_ems
    probs_ems(i,:) = ones(4,1)*.25;%p_op; % p_op = [0.3670; 0.2443; 0.1461; 0.2426]
end
EMS = [EMS, probs_ems];% Day_ems(ind), Week_ems(ind), Month_ems(ind), Year_ems(ind)];


% use every 5th point of EMS
% Data = [EMS(1:20:end,:); Opioid]; N_ems = size(EMS(1:20:end,:),1);
Data = [EMS;Opioid];

N = N_ems + N_op;
% remove Opioid points that are from areas not in EMS data
r = 0.01;
j = [];
for i = 1:N
    if Data(i,5) == 2
    xv = [Data(i,1)+r, Data(i,1)+r/2, Data(i,1)-r/2, Data(i,1)-r, Data(i,1)-r/2, Data(i,1)+r/2, Data(i,1)+r];
    yv = [Data(i,2), Data(i,2)+0.866*r, Data(i,2)+0.866*r, Data(i,2), Data(i,2)-0.866*r, Data(i,2)-0.866*r, Data(i,2)];
    in = inpolygon(EMS(:,1),EMS(:,2),xv,yv);
    if sum(in) == 0
        j = [j;i];
    end
    end
end
Data(j,:) = [];
N_op = N_op - length(j);
% standardize dataset
Data(:,1) = (Data(:,1)-min(Data(:,1)))/(max(Data(:,1))-min(Data(:,1))); % Lat
Data(:,2) = (Data(:,2)-min(Data(:,2)))/(max(Data(:,2))-min(Data(:,2))); % Long
Data(:,3) = Data(:,3)-min(Data(:,3)); % Time
% Data(:,10) = Data(:,10)-min(Data(:,10)); % Time of day
% Data(:,11) = Data(:,11)-min(Data(:,11)); % Day of week
% Data(:,12) = Data(:,12)-min(Data(:,12)); % Month of year
% Data(:,13) = Data(:,13)-min(Data(:,13)); % Year 

Data = sortrows(Data,3);


% h = histogram(EMS(:,4),24);
% h.Normalization = 'probability';
% prob_h = h.Values;
% h = cumsum([0, h.Values]);
% 
% d = histogram(Data(:,5),7);
% d.Normalization = 'probability';
% prob_d = d.Values;
% d = cumsum([0, d.Values]);
% 
% m = histogram(Data(:,6), 12);
% m.Normalization = 'probability';
% prob_m = m.Values;
% m = cumsum([0, m.Values]);
% 
% y = histogram(Data(:,7), 7);
% y.Normalization = 'probability';
% prob_y = y.Values;
% y = cumsum([0, y.Values]);
%% Initialization
N = N_ems+N_op;
tic
K0 = ones(4,1)*.5; w = ones(4,1)*.1; mu1 = ones(4,1)*50; mu2 = ones(4,1)*50; sigma = ones(4,1)*.01;
[K0, w, mu1,mu2, sigma, u, v, new_P] = EMH(Data, K0, w, mu1, mu2, sigma, 40, 0.1, 100);
toc


group = zeros(N,1);
for i = 1:N
    if Data(i,5)==1
%         [~,group2(i)] = max([sum(new_P{1}(i,:)), sum(new_P{2}(i,:)), sum(new_P{3}(i,:)), sum(new_P{4}(i,:))]);
        group(i) = sum(rand>=cumsum([0, sum(new_P{1}(i,:)), sum(new_P{2}(i,:)), sum(new_P{3}(i,:)), sum(new_P{4}(i,:))]));
    else
%         group2(i) = Data(i,4);
        group(i) = Data(i,4);
    end
end


ll_op = zeros(4,1);
ll_ems = zeros(4,1);
lambda = zeros(N,4);
for i = 1:N
    k = group(i);
    if Data(i,5) == 1 % unlabeled
        lambda(i,k) = mu(k)* u{k}(i) * v{k}(i); 
        for j=1:(i-1)
            lambda(i,k) =  lambda(i,k) + K0(k)*w(k)*exp(-w(k)*(Data(i,3)-Data(j,3)))*...
            exp((-(Data(i,1)-Data(j,1))^2-(Data(i,2)-Data(j,2))^2)/(2*sigma(k)^2))/(2*pi*sigma(k)^2);          
        end
        ll_ems(k)= ll_ems(k) + log(lambda(i,k));


    else %labeled
        lambda(i,k) = mu(k)* u{k}(i) * v{k}(i); 
        for j=1:(i-1)
            lambda(i,k) =  lambda(i,k) + K0(k)*w(k)*exp(-w(k)*(Data(i,3)-Data(j,3)))*...
            exp((-(Data(i,1)-Data(j,1))^2-(Data(i,2)-Data(j,2))^2)/(2*sigma(k)^2))/(2*pi*sigma(k)^2);          
        end
        ll_op(k)= ll_op(k) + log(lambda(i,k)) ;

    end
end 
for k = 1:4
    ll_ems(k) = ll_ems(k) - K0(k)*sum(sum(new_P{k}(Data(:,5)==1,:))) - mu1(k);
    ll_op(k) = ll_op(k) - K0(k)*sum(sum(new_P{k}(Data(:,5)==2,:))) - mu2(k);
end

tic
ll2_ems = 0;
lambda_ems = zeros(N_ems,1);
EMS = sortrows(Data(Data(:,5)==1,:), 3);
K0 = .5; w = .1; mu1 = 10; mu2 = 10; sigma = .01;
[K0, w, mu1, mu2, sigma, u,v, EMS_old_P] = EM(EMS, K0, w, mu1, mu2, sigma, 50, 0.01, 5);
for i = 1:size(EMS,1)
     lambda_ems(i) = mu1* u(i) * v(i); 
      
     for j=1:(i-1)
         lambda_ems(i) =  lambda_ems(i) + K0*w*exp(-w*(EMS(i,3)-EMS(j,3)))*...
                exp((-(EMS(i,1)-EMS(j,1))^2-(EMS(i,2)-EMS(j,2))^2)/(2*sigma^2))/(2*pi*sigma^2);          
     end 
     ll2_ems = ll2_ems + log(lambda_ems(i));
end
ll2_ems = ll2_ems - K0*N_ems - mu1;
toc

% labeled data only
ll3_op = zeros(4,1);
Op = Data(Data(:,5)==2,:);
K0 = ones(4,1)*.5; w = ones(4,1)*.1; mu1 = ones(4,1)*10; mu2 = ones(4,1)*10; sigma = ones(4,1)*.01;
for k = 1:4
    datatmp = Op(Op(:,4)==k,:);
    [K0(k), w(k), mu1(k),mu2(k) sigma(k), u, v, OP_old_P{k}] = ...
        EM(datatmp, K0(k), w(k), mu1(k), mu2(k), sigma(k), 50, 0.01, 5);
    for i = 1:size(datatmp,1)
       lambda = mu2(k)* u(i) * v(i); 
            for j=1:(i-1)
                lambda =  lambda + K0(k)*w(k)*exp(-w(k)*(datatmp(i,3)-datatmp(j,3)))*...
                exp((-(datatmp(i,1)-datatmp(j,1))^2-(datatmp(i,2)-datatmp(j,2))^2)/(2*sigma(k)^2))/(2*pi*sigma(k)^2);          
            end
        ll3_op(k) = ll3_op(k) + log(lambda);
    end
    ll3_op(k) = ll3_op(k) - K0(k)*size(datatmp,1) - mu2(k);
end

ll1(rep,:) = [sum(ll_op), sum(ll_ems)];
ll2(rep,:) = sum(ll2_ems);
ll3(rep,:) = sum(ll3_op);