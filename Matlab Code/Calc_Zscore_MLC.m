function [Zscore_mats, norm_elem_dotps] = Calc_Zscore_MLC(directory, i_fig_dir)
% -------------------------------
% This function calculates the dot product element by element for Zscores
% of each phase in learning, and is adjusted from the original function
% Calc_Zscore() to take data from the multiple learning cycle (MLC)
% simulations which have no post-NREM testing phase. - Michael
%
% Inputs:
% directory - directory of data files
% i_fig_dir - directory to save individual simulation figures
%
% NOTE it is essential that the data files in directory are ordered
% correctly. The order must follow: BB1 test phase 1, BB2 test phase 1, BB1
% test phase 2, ...


dir_struct = dir(strcat(directory,'*testphase*')); % Grabs info on all files 
% in directory with 'testphase' in their names. 

%names_bb1 = {dir_struct_bb1.name}; names_bb2 = {dir_struct_bb2.name};

filenames = {dir_struct.name}; % Combines bb1 and bb2 filenames
 
num_files = length(filenames);
Zscore_mats = cell(num_files, 1); % List for holding the Zscore matrices  
data_in_cell={};
for i = 1:num_files % Loops through each file, saving data.
    file = filenames{i};
    
    newstr = erase(file,'_spikes'); % Removes _spikes from end of filename and uses this for titling and saving the figures.
    phase = strrep(newstr,'_',' '); % Replaces any remaining '_' with a space
    michaels_spikes_cell = importdata(strcat(directory,file)); % Imports file data

    m_s_pre=michaels_spikes_cell{1};

    m_s_pre=m_s_pre(2:end-1);
    startIndex = regexp(m_s_pre,'\[(.*?)\]');

    
%     it=0;
%     Data_in_2=[];
%     for p =1:length(startIndex)-1
%       
%         it=it+1;
%         my_cur_times= m_s_pre(startIndex(p)+1:startIndex(p+1)-4);%m_s_pre_cell{p}(1:end-3);
%         my_cur_times=str2num(regexprep(my_cur_times,',',''));
%         if isempty(my_cur_times)% change -------------------5/19/2022-------------------
%           my_cur_times=[NaN] ; % Cause of Zscore matrix problem when my_cur_times=[0]. 
%         end
%         first_ind_2=ones(1,length(my_cur_times))*it;
%         a=size(first_ind_2);
%         b=size(my_cur_times);
%         data_add_2=[first_ind_2',my_cur_times'];
% 
%         Data_in_2=[Data_in_2;data_add_2];
% 
%         h=it;
%            
%       
%     end
it=0;
Data_in_2=[];
for p =1:length(startIndex)
    it=it+1;
    %below is code to read it out of michaels format
    if p<length(startIndex)
    my_cur_times= m_s_pre(startIndex(p)+1:startIndex(p+1)-4);%grabs each set of times seperated by brackets in michaels file
    else
    my_cur_times= m_s_pre(startIndex(p)+1:end-1);%grabs the last bracket at the end
    end
   
    my_cur_times=str2num(regexprep(my_cur_times,',',''));
          if isempty(my_cur_times)% change -------------------5/19/2022-------------------
          my_cur_times=[NaN] ; % Cause of Zscore matrix problem when my_cur_times=[0]. 
        end
    first_ind_2=ones(1,length(my_cur_times))*it; % this puts the index in the first column as we describbed was needed
    a=size(first_ind_2);
    b=size(my_cur_times);
    data_add_michael=[first_ind_2',my_cur_times'];
 
    Data_in_2=[Data_in_2;data_add_michael]; %this is the data in the correct out put
   
    h=it;
end
data_in_cell{i}=Data_in_2;
    %% 
    AMD_OUT_pre=AMDv4(Data_in_2);% calculates amd and isi's
    AMD_OUT=AMD_OUT_pre(6:end,:);% isolates amd matrix
%     AMD_OUT(1,1)=0;
    AMD_OUT(1,1)=NaN;% change -------------------5/19/2022-------------------

    

%     A=AMD_OUT;
%     figure
%     lowestValue = min(A(A(:)>0));
%     highestValue = max(A(:));
%     imagesc(A);
%     cmap = jet(256);
%     colormap(cmap);
%     caxis(gca,[lowestValue-2/256, highestValue]);
%     % Make less than lowest value black:
%     cmap(1,:)=[0,0,0];
%     colormap(cmap)
%     colorbar
%     title(strcat("AMD ", phase))
%     xlabel("Neuron ID")
%     ylabel("Neuron ID")
% %     empty_loc=find(A==0);
%     empty_loc=find(isnan(A)); % change -------------------5/19/2022-------------------
% %     saveas(figure(4),strcat(directory,'AMD_',phase,'.png'));

    AMD_1=AMD_OUT_pre;
    AMD_2=AMD_OUT_pre;% below calculates z score i reused cod typically amd 1 and amd are different
    [stability,z1,z2] = AMDv4_Similarity_c(AMD_1, AMD_2, min([min(AMD_1(1,:)) min(AMD_2(1,:))]), max([max(AMD_1(1,:)) max(AMD_2(1,:))]));
    % % % figure(2)
    % % % imagesc(z1)
    % % % title("Zscore")
    % % % colormap(cmap)
    % % % colorbar

    

    figure(1)
    B=z1; 
    empty_loc=find(isnan(B)); % Finds location of all NaN elements. change -------------------5/19/2022-------------------

    %B(empty_loc)=0
    % lowestValue = min(B(B(:)>0))
    %lowestValue = min(B(:))
    %highestValue = max(B(:))
    highestValue = 12; % Fixed values for max and min colorbar values to make them consistent between plots.
    lowestValue = -4;
    B(empty_loc)=-252; % Changes all NaN to -252 so they can be plotted black.
    
    imagesc(B);
    cmap = jet(256);
    colormap(cmap);
    caxis(gca,[lowestValue-2/256, highestValue]);
    % Make less than lowest value black:
    cmap(1,:)=[0,0,0];
    colormap(cmap)
    colorbar

    title(strcat("Zscore ", phase))
    xlabel("Neuron ID")
    ylabel("Neuron ID")
    saveas(figure(1),strcat(i_fig_dir,'Zscore_',phase,'.png'));
    
    B(B==-252)=0;% Changes all -252 elements in matrix to zero. change -------------------5/19/2022-------------------
    Zscore_mats{i} = B; % z1 is the matrix of Zscore data. This appends this matrix to the larger array for storage

end
%%

% Below accesses all the Zscore matrices produced above and takes the dot
% product between BB1 and BB2 matrices of the same phase in order to
% quantify orthogonality.

% % This relies on the order of appending to Zscore_mats. It is the same
% % order as the strings in filenames. Note that accessing Zscore_mats using
% % {} instead of () ensures that we access the matrices themselves, not the
% % cell objects that compose Zscore_mats. Using the function dot() on the
% % matrices takes the dot product between the corresponding column vectors
% % of the matrices. Thus the output is a vector with the dotp for each
% % column. To get the total value, I sum all these elements of the vector.
% pretest_columns = dot(Zscore_mats{1} , Zscore_mats{2});
% NREMtest_columns = dot(Zscore_mats{3} , Zscore_mats{4});
% posttest_columns = dot(Zscore_mats{5} , Zscore_mats{6});
% 
% % Summing elements of the vectors:
% pretest_dotp = sum(pretest_columns);
% NREMtest_dotp = sum(NREMtest_columns);
% posttest_dotp = sum(posttest_columns);



% Maybe here is a better way to do the dotp, element by element:

elem_dotps = zeros(length(filenames)/2 , 1); % to hold dot product sums for each phase
norm_elem_dotps = zeros(length(filenames)/2 , 1); % to hold normalized dot product sums for each phase
elem_squared = zeros(length(filenames) , 1); % to hold the sum of squared elements for each BB in each phase

mat_indices = 1:2:length(filenames); % Indicies for looping matrices in pairs
for i = 1:length(mat_indices)
    mat_i = mat_indices(i); % index for matrices
    
    for BB_row = 1:length(Zscore_mats{mat_i})
        for BB_col = 1:length(Zscore_mats{mat_i})
            % Diagonal elements of the matrices are already zero, so they
            % should not contribute to the dot product. If they are ever
            % not zero, it will be important to adjust code not to sum the
            % diagonal elements.
            BB1_elem = Zscore_mats{mat_i}(BB_row,BB_col); % Element of BB1 matrix
            BB2_elem = Zscore_mats{mat_i+1}(BB_row,BB_col); % Corresponding element of BB2 matrix
                    
            elem_dotps(i) = elem_dotps(i) + BB1_elem*BB2_elem;
            elem_squared(mat_i) = elem_squared(mat_i) + BB1_elem^2 ;
            elem_squared(mat_i+1) = elem_squared(mat_i+1) + BB2_elem^2 ;
                    
                    
        end
    end
    % Normalizzes dot products by the magnitudes of the matrices involved
    % and stores in norm_elem_dotps.
    norm_elem_dotps(i) = elem_dotps(i)/(sqrt(sum(elem_squared(mat_i)))*sqrt(sum(elem_squared(mat_i+1)))) ; 
    
end


figure(2)
%bp = bar(zipped_dotps);
bp = bar(norm_elem_dotps);
xticks(1:length(filenames)/2);

% This section creates labels for the xticks
phaselabels = {};
formatSpec = "Test Phase %d %s";
for phase = 1:length(filenames)/2
    phaselabels{phase} = sprintf(formatSpec, phase);
end

xticklabels(phaselabels);

title('MLC Zscore Matrix Dot Products')
%xlabel('Phase')
ylabel('Normalized Dot Products')
%ylim([-0.02,0.15]);
%legend({'Unnormalized','Normalized'});
%legend('Location','best')
saveas(figure(2),strcat(i_fig_dir,'MLC Zscore Matrices Dotp.png'));
% saveas(figure(2),strcat(i_fig_dir,'Zscore Matrices Dotp.pdf'));
