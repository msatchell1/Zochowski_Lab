function [BB1_diffs, BB2_diffs] = Zscore_hist(Zscore_mats, i_fig_dir)
% Plot the difference in Zscores between testing phases for each backbone
% as a histogram of the differences between each of the corresponding elements of the matrices.
% Zscore_mats is a list of matrices, each containing the Zscore values for
% a backbone and testing phase. The order of the list is as such: 
% [pretest_BB1, pretest_BB2, NREMtest_BB1, NREMtest_BB2, posttest_BB1,
% posttest_BB2].

BB1_mats = Zscore_mats(1:2:end); % Grabs all BB1 Zscore matrices.
BB2_mats = Zscore_mats(2:2:end); % Grabs all BB2 Zscore matrices.

BB1_diffs = {BB1_mats{2}-BB1_mats{1}, BB1_mats{3}-BB1_mats{2}, BB1_mats{3}-BB1_mats{1}}; % Calculates differences between phases and stores in a cell array.
BB2_diffs = {BB2_mats{2}-BB2_mats{1}, BB2_mats{3}-BB2_mats{2}, BB2_mats{3}-BB2_mats{1}};


figure(3)

subplot(2,3,1);
histogram(BB1_diffs{1}, Facecolor=[0,0,1]); % NREMtest - pretest for BB1
title('NREM BB1');
xlabel('Difference in Zscore');
% ylim([0,1000]);
% xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,2);
histogram(BB1_diffs{2}, Facecolor=[0,0,1]);
title('REM BB1');
% ylim([0,1000]);
% xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,3);
histogram(BB1_diffs{3}, Facecolor=[0,0,1]);
title('NREM + REM BB1');
% ylim([0,1000]);
% xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,4);
histogram(BB2_diffs{1}, Facecolor=[0,0.5,0]);
title('NREM BB2');
% ylim([0,1000]);
% xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,5);
histogram(BB2_diffs{2}, Facecolor=[0,0.5,0]);
title('REM BB2');
% ylim([0,1000]);
% xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,6);
histogram(BB2_diffs{3}, Facecolor=[0,0.5,0]);
title('NREM + REM BB2');
% ylim([0,1000]);
% xlim([-5,12]);
xline(0,'--r', LineWidth=1);

sgtitle('Change in Zscore For All LE Neurons'); % Title for whole figure

saveas(figure(3),strcat(i_fig_dir,'Zscore Histograms.png'));
% saveas(figure(3),strcat(i_fig_dir,'Zscore Histograms.pdf'));


% Another figure, this one only using the poritons of the Zscore matrix
% that contain blue and green LE neurons.
figure(4)

subplot(2,3,1);
histogram(BB1_diffs{1}(1:40,1:40), Facecolor=[0,0,1]); % NREMtest - pretest for BB1
title('NREM BB1');
xlabel('Difference in Zscore');
%ylim([0,1000]);
%xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,2);
histogram(BB1_diffs{2}(1:40,1:40), Facecolor=[0,0,1]);
title('REM BB1');
%ylim([0,1000]);
%xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,3);
histogram(BB1_diffs{3}(1:40,1:40), Facecolor=[0,0,1]);
title('NREM + REM BB1');
%ylim([0,1000]);
%xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,4);
histogram(BB2_diffs{1}(1:40,1:40), Facecolor=[0,0.5,0]);
title('NREM BB2');
%ylim([0,1000]);
%xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,5);
histogram(BB2_diffs{2}(1:40,1:40), Facecolor=[0,0.5,0]);
title('REM BB2');
%ylim([0,1000]);
%xlim([-5,12]);
xline(0,'--r', LineWidth=1);

subplot(2,3,6);
histogram(BB2_diffs{3}(1:40,1:40), Facecolor=[0,0.5,0]);
title('NREM + REM BB2');
%ylim([0,1000]);
%xlim([-5,12]);
xline(0,'--r', LineWidth=1);

sgtitle('Change in Zscore For Blue/Green LE Neurons');

saveas(figure(4),strcat(i_fig_dir,'Zscore Histograms BG Only.png'));
% saveas(figure(4),strcat(i_fig_dir,'Zscore Histograms BG Only.pdf'));
