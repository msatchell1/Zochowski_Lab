function [ AMD ] = AMD_ij( S1,S2 )
%AMD_ij return the AMD from two time series
%   Inputs are assumed to be sorted vectors indicatind the time of the
%   spikes

%calculating average minimum distance
% 1/N *sum(delta t_k i)_k=1 Ni

% Spike procedure
% --------------------------
% walk through the spikes in S1 first.  Find the first spike in S2 greater than that value.  
% Compare the first spike greater and the last spike before to the spike
% in S1 and pick the one which has the smaller temporal distance.


% find the distance from 1==>2
% -----------------------------------------
AMD1=0;
%first step is to walkthrough the spikes in S1 and find the closest spike in S2
ind = 1; %tracks the spike in S2
for ii=1:length(S1)
    t_k=S1(ii);
    t_j=S2(ind);
    while t_k>t_j %increase the indice in S2 till its the spike occurs after t_k
        if ind>=length(S2) %check that the indice doesn't exceed S2
            break % if it does exit this loop
        end
        ind=ind+1; % move to the next spike in S2
        t_j=S2(ind);
    end
    
    if ind==1 %special condition because  ind-1 doesn't exist in the vector
        delta_t=abs(t_j-t_k);%no comparison necessary if the first spike in j is greater than that spike in k
        AMD1=AMD1+delta_t;
        continue
    end
    
    % compare nearest spike on each side, pick the closer one
    if abs(t_j-t_k) < abs(t_k-S2(ind-1))
        delta_t=abs(t_j-t_k);
    else
        delta_t=abs(t_k-S2(ind-1));
    end    
    AMD1=AMD1+delta_t;
end
AMD1=AMD1/length(S1);


% find the distance from 2==>1
% ------------------------------------------
AMD2=0;
%first step is to walkthrough the spikes in S1 and find the closest spike in S2
ind = 1; %tracks the spike in S2
for ii=1:length(S2)
    t_k=S2(ii);
    t_j=S1(ind);
    while t_k>t_j %increase the indice in S2 till its the spike occurs after t_k
        if ind>=length(S1) %check that the indice doesn't exceed S2
            break % if it does exit this loop
        end
        ind=ind+1; % move to the next spike in S2
        t_j=S1(ind);
    end
    
    if ind==1 %special condition because  ind-1 doesn't exist in the vector
        delta_t=abs(t_j-t_k);%no comparison necessary if the first spike in j is greater than that spike in k
        AMD2=AMD2+delta_t;
        continue
    end
        
    if abs(t_j-t_k) < abs(t_k-S1(ind-1))
        delta_t=abs(t_j-t_k);
    else
        delta_t=abs(t_k-S1(ind-1));
    end    
    AMD2=AMD2+delta_t;
end
AMD2=AMD2/length(S2);


% return both values
% --------------------------
AMD=[AMD1, AMD2];
end

