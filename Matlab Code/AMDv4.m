function [AMD_out] = AMDv4( Data_in, fpath, fname, fNum)
%AMDV4 Complete Average Minimum Distance Analysis
%   Compares the distance of the spike trains to the expected distance.

% Input Structure: 
% Data_in: two column data file,
% first column cell ID, second column spike times.
% The 1st column is sorted by cell ID, the spikes of the corresponing IDS
% are sorted from earliest to latest.  The begining and end of file
% correspond to the t_min and t_max of the data
% f_in is the filename in, append '_AMD' to it before outputting.  Assuming
% it is a .dat or .txt input(the code subtract off last foour characters of the name
% before appending _v4_AMDv4.dat to it)

 
% Program call: AMDv4(file1,fname_example.dat)
% Output: fname_v4_AMD.dat
% this program wil typically be called by a program that automatically
% reads in data files( using the command dir)
 
% Output structure for files with N cells: 
%  row 1 - Cell IDs
%  row 2 - Spike Counts
%  row 3 - ISI mean
%  row 4 - ISI std.dev
%  row 5 - Poison mean and Std. dev
%  rows 6-(6-1+n) - Actual average minimum distances
%  row6, column1 - t_max-t_min

% Program flow:
% 1)Read in file, get cell IDs and spike counts posion values, time
% duration
% 2)Loop through spikes to find expected ISI values
% 3)Send pairs of spike times into the function AMD_ij
% 4)output file

% Expected Values ij(from spikes in train i to nearest spike in train j) 
% are independent of i.  Dependent only on the distribution j
% Poison:
% expected poison mean and width ij = (1/2)*(Total_t/N_j) - an averge 
% interval, T/N_j, should be sampled at 1/4, the poison correction accounts for larger intervals
% ISI:
% expected ISI_ij = sum((1/2)*(ISI/2)*(ISI/Total_t))
% expected ISI^2_ij = sum((1/3)*(ISI/2)^2*(ISI/total_t))
% for end intervals replace ISI/2 by ISI since there is no spike on one end


% 1)Read in file, get cell IDs, time duartion, spike counts, posion values
% -------------------------------------------------------------------------
IDs=unique(Data_in(:,1));% this sorts it, but the input should already be sorted 
N=length(IDs);
AMD_out=zeros(N+5,N);
AMD_out(1,:)=IDs;

% calculate spike counts
% -----------------------------
% Create array that indicates where the new cell start firing
tmp=diff(Data_in(:,1)); %not 0 indicates a new cell firing
%len
Spikes=find(tmp~=0);
Spikes=Spikes+1; %shift to account for differencing
Spikes=[1; Spikes; 1+size(Data_in,1)];
counts=diff(Spikes);
AMD_out(2,:)=counts;

% calculate t_max, t_min
t_max=max(Data_in(:,2));
t_min=min(Data_in(:,2));
Total_t=t_max-t_min;

% calculate poison values
Poison=(1/2)*Total_t ./(counts+1); %???? count or count+1 ????? which represent the ISIs?
AMD_out(5,:)=Poison;


% 2)Loop through spikes to find expected ISI values
% -------------------------------------------------------------------
ISI_mean=zeros(1,N);
ISI_width=zeros(1,N);
%spikes keeps track of where a new indice beings
for ii=1:N
   ISI=Data_in(Spikes(ii):(Spikes(ii+1)-1),2);
   
   %account for ISIs on the end part a)
   ISI_1=ISI(1)-t_min;
   ISI_2=t_max-ISI(end);
   ISI_mean(ii)=ISI_mean(ii)+(1/(2*Total_t))*ISI_1^2;
   ISI_mean(ii)=ISI_mean(ii)+(1/(2*Total_t))*ISI_2^2;
   ISI_width(ii)=ISI_width(ii)+(1/(3*Total_t))*ISI_1^3;
   ISI_width(ii)=ISI_width(ii)+(1/(3*Total_t))*ISI_2^3;

   % account for ISI's between the spikes
   ISI=diff(ISI);
   if( sum(tmp<0) ~=0)
      'warning: spikes are not being ordered correctly, algorithm failure' 
   end %this isn't needed, but is a sanity check that its working
   ISI_mean(ii)= ISI_mean(ii)+ sum((1/(4*Total_t))*ISI.^2);
   ISI_width(ii)= ISI_width(ii)+ sum((1/(12*Total_t))*ISI.^3);
end
ISI_width=sqrt(ISI_width-ISI_mean.^2); % std.dev= <x^2>-<x>^2
AMD_out(3,:)=ISI_mean;
AMD_out(4,:)=ISI_width;


% 3)Send pairs of spike times into the function AMD_ij
% -------------------------------------------------------------------
Real_Dist=zeros(N,N);
for ii=1:N
    S1=Data_in(Spikes(ii):(Spikes(ii+1)-1),2);
    Real_Dist(ii,ii)=NaN; % change -------------------5/19/2022-------------------

    for jj=(ii+1):N
        S2=Data_in(Spikes(jj):(Spikes(jj+1)-1),2);
        tmp=AMD_ij(S1,S2);% calculate average time difference
        
        Real_Dist(ii,jj)=tmp(1);
        Real_Dist(jj,ii)=tmp(2);
        if (AMD_out(2,ii)<3) || (AMD_out(2,jj)<3)
         %Real_Dist(ii,jj)=0; % ge rid of nurons with 2 spikes or less   
         %Real_Dist(jj,ii)=0;   
         Real_Dist(ii,jj)=NaN; % ge rid of nurons with 2 spikes or less   % change -------------------5/19/2022-------------------
         Real_Dist(jj,ii)=NaN;   % change -------------------5/19/2022-------------------
        end
    end
    
end
AMD_out(6:end,:)=Real_Dist;

%save Total_t in empty space in AMD_out
AMD_out(6,1)=Total_t;


% 4)output file
% -------------------------------------------------------------------
%str1=[fpath fname(1:end-4),'_',num2str(fNum), '_v4_AMD.dat'];
%dlmwrite(str1,AMD_out,'delimiter','\t','precision','%.5f');
end



% counts=zeros(1,N);
% for ii=1:N
%     counts(ii)=sum(Data_in(:,1)==IDs(ii);
% end

%    kk_max=counts(jj)-1;
%    for kk=1:kk_max
%        
%    end

% Poison -  Poison_mean[i] = total_t/(2*(spikes[i+1]-spikes[i]+1));//same
% as the mean for poison distribution
% ISI mean, and <x^2>
% 	  temp +=.25*ISI*ISI/total_t  ;// 1/4 ISI gives the average on the interval, ISI/Total_t is the chance of falling on that interval, total_t i
% 	  temp2 += (1/double(12))*ISI*ISI*ISI/(total_t);
%       
%       on ends
%             // add the values for the first and last ISIs, they may be zero if the spike corresponds to a min or max
%       ISI = T[spikes[ii]]-min_t;
%       if ( ISI <0){cout << "first ISI = " << ISI << ",   ";} // just checks to see if this fails somehow, should be able to remove
%       temp +=.5*ISI*ISI/total_t;//.5 on end ISI's because no spike capping one of the ends
%       temp2 += (1/double(3))*ISI*ISI*ISI/total_t;// (1/4) of the two-sided intervals
%       ISI = max_t-T[spikes[ii+1]-1];
%       if ( ISI <0) {cout << "last ISI = " << ISI << ",   ";}
%       temp +=.5*ISI*ISI/total_t  ;
%       temp2 += (1/double(3))*ISI*ISI*ISI/total_t;
%       
%       Rescaling_values[ii]=temp;
%       ISI_widths[ii]=temp2;
