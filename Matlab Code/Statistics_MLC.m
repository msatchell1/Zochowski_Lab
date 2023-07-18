function [] = Statistics_MLC()
% This function will generate statistical plots for multiple learning cycles (MLC).
% For now, it only does generates the plot for Zscore dot products.

file_directory = 'C:\Users\micha\OneDrive - brandeis.edu\Umich Stuff\Zochowski Lab\Results\Multiple Learning Cycles\My Runs\vary num_lrn_cycles';
file_directory = strcat(file_directory, '\');

s_fig_directory = 'C:\Users\micha\OneDrive - brandeis.edu\Umich Stuff\Zochowski Lab\Results\Multiple Learning Cycles\Figures\My Runs\vary num_lrn_cycles\statistical plots'; % Directory to store statistical figures.
s_fig_directory = strcat(s_fig_directory, '\');

i_fig_directory = 'C:\Users\micha\OneDrive - brandeis.edu\Umich Stuff\Zochowski Lab\Results\Multiple Learning Cycles\Figures\My Runs\vary num_lrn_cycles\indv sim plots'; % Directory to store individual simulation figures.
i_fig_directory = strcat(i_fig_directory, '\');

s_filenames = dir(file_directory); % Holds the names of subfolders within file_directory. These should be the folders named 'seed = X'.
i_filenames = dir(i_fig_directory); 

s_filenames = s_filenames([3:end]); % Removes the '.' and '..' folders from the structure array.
i_filenames = i_filenames([3:end]);


% for i = 1:length(s_filenames) % Shows the names of the folders in the data directory.
%    s_filenames(i).name
% end
  

s_Zscore_mats = cell(length(s_filenames),1); % For storing data from all sims.
%s_Zscore_mats = double.empty(length(s_filenames),0); % I should change
%these cell arrays over to numerical "double" arrays when I have the
%chance. This is better practice.
s_mat_diffs = cell(length(s_filenames),1);

% s_norm_dotps = zeros(length(i_filenames)/2,length(s_filenames));


for i = 1:length(s_filenames) % Loops through files, calculating Zscore matrices from each one and storing it in s_Zscore_mats.
    
    file = s_filenames(i).name;
    i_file = i_filenames(i).name;
    disp(strcat(file_directory,file,'\')) % Prints from which file the data is being loaded.
    
    [Zscore_mats, norm_elem_dotps] = Calc_Zscore_MLC(strcat(file_directory,file,'\'), strcat(i_fig_directory,i_file,'\')); % Extracts Zscore matrix for each sim, as well as dotps.
    s_Zscore_mats{i} = Zscore_mats; % Stores Zscore matrix for each sim.
    s_norm_dotps(i,:) = norm_elem_dotps;
    
%     [BB1_diffs, BB2_diffs] = Zscore_hist_MLC(Zscore_mats, strcat(i_fig_directory,i_file,'\')); % Saves change in Zscore histograms and returns data used in the plots.
%     
%     s_mat_diffs{i} =  [BB1_diffs, BB2_diffs]; % Stores Zscore matrix phase differences for each sim. 

end




% This code calculates statistics on error for the Zscore dotp plot, and
% then plots it.

s_mean_dotps = mean(s_norm_dotps, 2); % Averages each row of the matrix, giving the average values for each dotp for all sims.
deviations = s_norm_dotps - s_mean_dotps; % Deviations of each value from the means.
dev_sqrd = deviations.^2; % Square these values.
sum_dev_sqrd = sum(dev_sqrd, 2)./(length(s_filenames)-1); % Sum values from each sim and divide by number of sims-1.
std_devs = sqrt(sum_dev_sqrd); % Takes square root of each element. This gives the standard deviations.
std_error_ofmeans = std_devs/sqrt(length(s_filenames)); % The standard error of the mean.


figure(5)

%bp = bar(zipped_dotps);
bp = bar(s_mean_dotps);
xticks([1,2,3]);
xticklabels({'Pre-Learning','After NREM','After NREM + REM'});
title('s Zscore Matrix Dot Products')
%xlabel('Phase')
ylabel('Normalized Dot Products')
%ylim([-0.02,0.15]);
%legend({'Unnormalized','Normalized'});
%legend('Location','best')

hold on

x_dummy = 1:length(s_mean_dotps); % Dummy x vals for errorbar plotting.
er = errorbar(x_dummy, s_mean_dotps, [], std_error_ofmeans); % Errorbar plots a line on its own.
er.Color = ([0,0,0]);
er.LineStyle = 'none'; % Removes line from errorbar()
er.CapSize = 20;

saveas(figure(5),strcat(s_fig_directory,'s Zscore Matrices Dotp.png'));
saveas(figure(5),strcat(s_fig_directory,'s Zscore Matrices Dotp.pdf'));

hold off
