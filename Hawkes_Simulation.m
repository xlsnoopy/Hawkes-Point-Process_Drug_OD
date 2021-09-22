function [times x y p]=Hawkes_Simulation(mu,k0,w,T,sig, bg)
%this program simulates event times, called "times", according
%to a self-exciting point process with exponential triggering kernel and
%parameters mu (const. background rate)
%k0 (branching ratio) and w (exp parameter) on the time interval [0,T]
%it also simulates event locations, called "x" and "y",
%with a background vector (bg)


times=zeros(5000,1);
x=zeros(5000,1);
y=zeros(5000,1);
%first simulate "background" events
%this is done by picking p points where p is Poisson with parameter mu*T
%and then distributing the points uniformly in the interval [0,T]
p=pois(mu*T);
times(1:p,1)=rand(p,1)*T;
for i = 1:p
    rx = rand;
    if rx <= (bg(1)+bg(3))/sum(bg)
        x(i,1)= .5*rand;
    else
        x(i,1) = .5 + .5*rand;
    end
    ry = rand;
    if ry <= (bg(3)+bg(4))/sum(bg)
        y(i,1) = .5*rand;
    else
        y(i,1) = .5 + .5*rand;
    end

end

counts=1;
countf=p;

%Next loop through every event and simulate the "offspring"
%even the offspring events can generate their own offspring

while((countf-counts)>-1)
p0=pois(k0);  %each event generates p offspring according to a Poisson r.v. with parameter k0
for j=1:p0
    temp=times(counts)-log(rand())/w; % this generates an exponential r.v. on [t_counts,infty]
    temp2=x(counts)+sig*randn(); % inter-point distances are gaussian
    temp3=y(counts)+sig*randn();
    if(temp<T)                        % we only keep this time if it is in [t_counts,T]
        countf=countf+1;
        times(countf)=temp;
        x(countf)=temp2;
        y(countf)=temp3;
    end
end
counts=counts+1;
end
data=[times(1:countf) x(1:countf) y(1:countf)];
data=sortrows(data,1); %sort by times
times=data(:,1);
x=data(:,2);
y=data(:,3);

end


