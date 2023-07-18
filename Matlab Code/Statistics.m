%function [] = Statistics()

% This function will generate statistical plots, specifically the summation
% of Zscore histograms from many different simulations. 

file_directory = 'C:\Users\satchelm\OneDrive - Umich\Desktop\Zochowski Lab\STATISTICS Data\Seed = 4\vary Idrive_LE\NREM-REM'; % Directory of data files.
file_directory = strcat(file_directory, '\');

s_fig_directory = 'C:\Users\satchelm\OneDrive - Umich\Desktop\Zochowski Lab\STATISTICS Data\Seed = 4\vary Idrive_LE\Figures\NREM-REM'; % Directory to store statistical figures.
s_fig_directory = strcat(s_fig_directory, '\');

i_fig_directory = 'C:\Users\satchelm\OneDrive - Umich\Desktop\Zochowski Lab\STATISTICS Data\Seed = 4\vary Idrive_LE\Figures\NREM-REM\Individual Sim Plots'; % Directory to store individual simulation figures.
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

s_norm_dotps = zeros(3,length(s_filenames));


for i = 1:length(s_filenames) % Loops through files, calculating Zscore matrices from each one and storing it in s_Zscore_mats.
    
    file = s_filenames(i).name;
    i_file = i_filenames(i).name;
    disp(strcat(file_directory,file,'\')) % Prints from which file the data is being loaded.
    
    [Zscore_mats, norm_elem_dotps] = Calc_Zscore(strcat(file_directory,file,'\'), strcat(i_fig_directory,i_file,'\')); % Extracts Zscore matrix for each sim, as well as dotps.
    s_Zscore_mats{i} = Zscore_mats; % Stores Zscore matrix for each sim.
    s_norm_dotps(:,i) = norm_elem_dotps;
    
    [BB1_diffs, BB2_diffs] = Zscore_hist(Zscore_mats, strcat(i_fig_directory,i_file,'\')); % Saves change in Zscore histograms and returns data used in the plots.
    
    s_mat_diffs{i} =  [BB1_diffs, BB2_diffs]; % Stores Zscore matrix phase differences for each sim. 

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





% This code deals with the histogram plots.

for i = 1:size(s_mat_diffs, 1) % Loops through length of third dimension of diffs_array, which should be equal to length(s_filenames). 
    i_diffs_cell = s_mat_diffs{i}; % Matrices for one sim.
    
    test_double_array = reshape( cat(3, i_diffs_cell{:}), [80,80,6]); % Turns cell array into double array.
    s_mat_diffs{i} = test_double_array; % Puts back in cell array.
        
end

s_double_diffs = reshape( cat(4, s_mat_diffs{:}), [80,80,6,length(s_filenames)]); % Turns whole cell array into double array.




% Plots statistical histogram data for all sims combined.

figure(6)

subplot(2,3,1);                                                                                                                              
histogram(s_double_diffs(:,:,1,:), Facecolor=[0,0,1]); % NREMtest - pretest for BB1
title('NREM BB1');
xlabel('Difference in Zscore');
% ylim([0,1000]);
% xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,2);
histogram(s_double_diffs(:,:,2,:), Facecolor=[0,0,1]);
title('REM BB1');
% ylim([0,1000]);
% xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,3);
histogram(s_double_diffs(:,:,3,:), Facecolor=[0,0,1]);
title('NREM + REM BB1');
% ylim([0,1000]);
% xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,4);
histogram(s_double_diffs(:,:,4,:), Facecolor=[0,0.5,0]);
title('NREM BB2');
% ylim([0,1000]);
% xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,5);
histogram(s_double_diffs(:,:,5,:), Facecolor=[0,0.5,0]);
title('REM BB2');
% ylim([0,1000]);
% xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,6);
histogram(s_double_diffs(:,:,6,:), Facecolor=[0,0.5,0]);
title('NREM + REM BB2');
% ylim([0,1000]);
% xlim([-5,12]);
xline(0,'--r', LineWidth=1);

sgtitle('s Change in Zscore For All LE Neurons'); % Title for whole figure

saveas(figure(6),strcat(s_fig_directory,'s Zscore Histograms.png'));
saveas(figure(6),strcat(s_fig_directory,'s Zscore Histograms.pdf'));


% Another figure, this one only using the poritons of the Zscore matrix
% that contain blue and green LE neurons.
figure(7)

subplot(2,3,1);
histogram(s_double_diffs(1:40,1:40,1,:), Facecolor=[0,0,1]); % NREMtest - pretest for BB1
title('NREM BB1');
xlabel('Difference in Zscore');
%ylim([0,1000]);
%xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,2);
histogram(s_double_diffs(1:40,1:40,2,:), Facecolor=[0,0,1]);
title('REM BB1');
%ylim([0,1000]);
%xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,3);
histogram(s_double_diffs(1:40,1:40,3,:), Facecolor=[0,0,1]);
title('NREM + REM BB1');
%ylim([0,1000]);
%xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,4);
histogram(s_double_diffs(1:40,1:40,4,:), Facecolor=[0,0.5,0]);
title('NREM BB2');
%ylim([0,1000]);
%xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,5);
histogram(s_double_diffs(1:40,1:40,5,:), Facecolor=[0,0.5,0]);
title('REM BB2');
%ylim([0,1000]);
%xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,6);
histogram(s_double_diffs(1:40,1:40,6,:), Facecolor=[0,0.5,0]);
title('NREM + REM BB2');
%ylim([0,1000]);
%xlim([-5,12]);
xline(0,'--r', LineWidth=1);

sgtitle('s Change in Zscore For Blue/Green LE Neurons');

saveas(figure(7),strcat(s_fig_directory,'s Zscore Histograms BG Only.png'));
saveas(figure(7),strcat(s_fig_directory,'s Zscore Histograms BG Only.pdf'));


