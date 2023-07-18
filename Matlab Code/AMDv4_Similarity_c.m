function [ Similarity,Z_score_1,Z_score_2] = AMDv4_Similarity_c( AMD1,AMD2,r_min,r_max )
%UNTITLED Calculate the similarity between two AMD output files
%   Program Flow - 
%     1) use unique to generate the list of unique cells between the two files
%     2) construct the matries to compare by using list of cells to order the locations
%     3) Remove the diagonals, which will result in vectors.
%     4) compare as needed.

% Program Inputs  - AMD output files are (n+2)Xn sized matrices.  The first
% row is the cell ID's, the second is the number of spikes for that cell.
% The bottom nXn matrix are the pairwise values.  The diagonal is 0 except
% for the first entry which is the time of the simulation

% values computed - 
% 1a) (1-vals) sqrt(n) corrected ISI scores
% 1b) (1-vals) non corrected ISI dot product
% 1c) (1-vals) RMSD
% 1d) (1-vals) RMSD of deviation from mean(of that slice)
%

% 2a-d) Poison values


%version c, only compares incdices in a given range


% AMDv4_values to compute - 
% stabilities for ISI z-scores
% a)cosine similarity
% b) cosine similarity no sqrt(n) correction
% c) RMSD - root mean squared difference(this value should go to 0 if similar)
% d) RMSD in deviation from mean at that time

Cells1=AMD1(1,:);
Counts1 = AMD1(2,:);
ISI_means1 = AMD1(3,:);
ISI_widths1 = AMD1(4,:);
% Poison_means1=AMD1(5,:);
% Time1=AMD1(6,1);
Pairwise1=AMD1((6:end),:);
Pairwise1(1,1)=0;
% Pairwise1(1,1)=NaN; % change -------------------5/19/2022-------------------

%) Trim values to selected range
I=(Cells1>= r_min & Cells1<= r_max);% indices to keep
if (sum(I)==0) %no cells in given range
   Similarity=[nan]; %return nan
   return
end
Cells1=Cells1(I);
Counts1 =Counts1(I);
ISI_means1=ISI_means1(I);
ISI_widths1=ISI_widths1(I);
Pairwise1=Pairwise1(I,I);
n=length(Cells1);%this gives dimension of matrix) 
%create the ISI z-scores
I_vals1=(repmat(ISI_means1,n,1)-Pairwise1)./repmat(ISI_widths1,n,1);%no sqrt(n) correction
IN_vals1=I_vals1.*repmat(sqrt(Counts1'),1,n); %sqrt n correction


Cells2=AMD2(1,:);
Counts2 = AMD2(2,:);
ISI_means2 = AMD2(3,:);
ISI_widths2 = AMD2(4,:);
% Poison_means2=AMD2(5,:);
% Time2=AMD2(6,1);
Pairwise2=AMD2((6:end),:);
Pairwise2(1,1)=0;
% Pairwise2(1,1)=NaN; % change -------------------5/19/2022-------------------
%) Trim values to selected range
I=(Cells2>= r_min & Cells2<= r_max);%
if (sum(I)==0) %no cells in given range
   Similarity=[nan];
   return
end
Cells2=Cells2(I);
Counts2 =Counts2(I);
ISI_means2=ISI_means2(I);
ISI_widths2=ISI_widths2(I);
Pairwise2=Pairwise2(I,I);
n=length(Cells2);%this gives dimension of matrix) 
%create the ISI z-scores
I_vals2=(repmat(ISI_means2,n,1)-Pairwise2)./repmat(ISI_widths2,n,1);%no sqrt(n) correction
IN_vals2=I_vals2.*repmat(sqrt(Counts2'),1,n); %sqrt n correction




% make a list of the unique values in Cells
% -----------------------------------------------------
Cells = unique([Cells1, Cells2]); %returns the order locations
n=length(Cells);%this gives dimension of matrix) 


% create mapping from values in Cells1 to Cells
% -----------------------------------------------------
map1=zeros(1,length(Cells1));
for ii=1:length(Cells1)
    map1(ii)=find(Cells==Cells1(ii));
end
%Fill the vals matrix
IVals1=zeros(length(Cells));
IVals1_noN=zeros(length(Cells));
IVals1(map1,map1)=IN_vals1; %place the values from the matrix in accordingly
IVals1(1:(n+1):end)=[];%reshape by deleting diagonal elements
IVals1_noN(map1,map1)=I_vals1; %place the values from the matrix in accordingly
IVals1_noN(1:(n+1):end)=[];%reshape by deleting diagonal elements


%create mapping from values in Cells2 to Cells
% -----------------------------------------------------
map2=zeros(1,length(Cells2));
for ii=1:length(Cells2)
    map2(ii)=find(Cells==Cells2(ii));
end
%Fill the vals matrix
IVals2=zeros(length(Cells));
IVals2_noN=zeros(length(Cells));
IVals2(map2,map2)=IN_vals2; %place the values from the matrix in accordingly
IVals2(1:(n+1):end)=[];
IVals2_noN(map2,map2)=I_vals2; %place the values from the matrix in accordingly
IVals2_noN(1:(n+1):end)=[];


% ---------------------------------------------------------
%  Compute desired comparisons between the matrices.
% ---------------------------------------------------------
%) stabilities for ISI z-scores
% 1a)cosine similarity
% 1b) cosine similarity no sqrt(n) correction
% 1c) RMSD - root mean squared difference(this value should go to 0 if similar)
% 1d) RMSD of Deviation from 



Z_score_1=IN_vals1;
Z_score_2=IN_vals2;
out11 = dot(IVals1,IVals2)/sqrt(dot(IVals1,IVals1)*dot(IVals2,IVals2));
out12 = dot(IVals1_noN,IVals2_noN)/sqrt(dot(IVals1_noN,IVals1_noN)*dot(IVals2_noN,IVals2_noN));
tmp = IVals1-IVals2;
out13=sqrt(mean(tmp.^2));%sqrt of the mean of the squared differences
tmp = (repmat(mean(IVals1),size(IVals1))-IVals1)-(repmat(mean(IVals2),size(IVals2))-IVals2);
out14=sqrt(mean(tmp.^2));

% Similarity 1:4 [ ISI Z-score similarity, ISI Z-score no sqrt(n) correction, ISI Z-scores RMSD, RMSD of devaition from average at that time]
Similarity = [ out11, out12, out13, out14];
end


% Hist_vals=[reshape(I_vals1,[],1) reshape(IN_vals1,[],1) reshape(P_vals1,[],1) reshape(PN_vals1,[],1)]; %put all values of AMD1 into a matrix
% n=length(Cells1);%set n back to the size of the matrix of AMD1
% Hist_vals([1:(n+1):end],:)=[];%remove the diagonals of each matrix all in 1 step
% %current setup will miss the values from the very last file.  Could augment
% %this by passing in more variables or returning extra data.
