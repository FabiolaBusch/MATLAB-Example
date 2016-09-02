%% Data_Analysis.m
% 
% This is the type of scripts I work with. In this example, I implemented an
% analysis to identify an analyze specific epoch transitions. 
% 
% Some function and variable names have been modified to the original
% script. The additionally used functions for data processing and mean cal-
% culation were made by colleagues.
%
% % This script analyses Data from a Stimulus wich shows two different epochs.
% %
% % The data which is examined (called 'rats' or 'ratio') is a flourescent 
% % calcium Signal, measured by a 2-photon microscope.
% %

% %
% % Fabiola 8/2016


clear all;
close all;
clc;

% Change this variable to your pData path
pdatapath='C:\path\to\processed\data';
addpath(pdatapath);
cd(pdatapath);

% For the code
addpath('C:\path\to\code');

% fetch Information from Masterfoldersummary
database_select_samples;

%% Processing Data 


% find matching pData files
matching_data = find(used_stimtype .* interesting_genotype .*~i_moving);

% create structures by neuron
data_struct = create_neuron_structure_all(matching_data);

% load all the data!
loaded_data = load_neuron_data10Hz(data_struct,pdatapath);

% This function groups the different interesting epoch combinations in a
% Matrix, stored in combination_storage
grouped_data = aggregate_Data(loaded_data);

%% Printing ALL Traces

THRESHOLD = 1;
fprintf('THRESHOLD for Signal sum: %d \n',THRESHOLD );

% take a mean of all assigned values, ignore 0
crossing_rats = grouped_data.rats;
crossing_rats(crossing_rats==0) = NaN;

% Mean over all ROI's (Region Of Interest -> Neurons)
mean_val = squeeze(mean(crossing_rats,1,'omitnan'));

plots = size(grouped_data.rats,2);

% Modiefied from original
titles = {'1';'2';'3';'4';'5';'6''7';'8';'9''10';'11';'12''13';'14';'15''16';'17';'18'};


% Watch over 20 seconds period
cur_t = linspace(0,20,length(mean_val(1,:)))';

for ii = 1:plots
    
    if abs(sum(mean_val(ii,:))) > THRESHOLD
        fprintf('Signal sum %d: %d \n',ii,abs(sum(mean_val(:,ii))));
        figure(ii);
        hold on;
        
        epochRats = squeeze(crossing_rats(:,ii,:));
        [x,m,e] = mean_cat_full(epochRats,1,grouped_data.flyID); %m=mean; e=std
        h1 = plot_err_patch_v2(cur_t,m,e,[0 0 1],[0.8 0.8 1]);
        
        xlabel('time (s)');
        ylabel('dF/F - Calcium Signal');
        ylim([-0.4 1.4]);
        niceaxes;
        title(titles(ii));
        
        line([0 20],[0 0],'color',[0 0 0]);
        line([0 10 10 20],[0.7 0.7 0.85 0.85],'color',[0 0 0]);
  
        hold off;
    else
        fprintf('Case %d did not cross THRESHOLD (%s).\n',ii,titles{ii});
    end
end
%% Evaluating the correltions between same contrast differences
close all;

num_bars = size(mean_val,1);

maximum = zeros(1,num_bars);
minimum = zeros(1,num_bars);
max_data = zeros(1,num_bars);
min_data = zeros(1,num_bars);
err_max = zeros(1,num_bars);
err_min = zeros(1,num_bars);

% Search max-Signal per ROI and take the mean of all max"s per ROI signals
for kk=1:num_bars
    % find maximum in each trace
    max_data = max(crossing_rats(:,kk,:));
    min_data = min(crossing_rats(:,kk,:));
    
    % take the mean of all max's
    maximum(kk) = mean(max_data,'omitnan');
    minimum(kk) = mean(min_data,'omitnan');
    
    % Standard error
    err_max(kk) = std(max_data,'omitnan')/sqrt(length(max_data));
    err_min(kk) = std(min_data,'omitnan')/sqrt(length(min_data));
end

abs_max = [maximum, minimum];
abs_max = abs(abs_max);

% When does the full change to stim_one?
stim_one = [1 1 1 0 1 1 1 0 0 1 1 0 0 0 1 0 0 0];
stim_two = ~ stim_one;

% Plot ALL Maxima/ Minima-------------------------------------------------
bar_all = figure;
grid on;

% Plot the maxima only for steps_two
y1 = maximum(stim_two);
x1 = find(stim_two);
e1 = err_max(stim_two);

hold on;
bar(x1,y1,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.0);
errorbar(x1,y1,e1,'rx');

% Plot the minima only for stim_one
y2 = abs(minimum(find(stim_one)));
x2 = find(stim_one);
e2 = err_min(find(stim_one));

bar(x2,y2,'FaceColor',[.5 .5 0],'EdgeColor',[.9 .9 0],'LineWidth',1.0);
legend({'Stim two - Maxima','Error','Stim One - |Minima|'},'Location','northwest');
errorbar(x2,y2,e2,'rx');

title('All Stim - Signal strength');

% Add titles
ax = gca;
ax.XTick = 1:18;
ax.XTickLabels = {'1','2','3','4','5','6''7','8','9''10','11','12','13','14','15','16','17','18'};
ax.XTickLabelRotation = 45;

hold off;
%% Plot by category - 1 Step ---------------------------------------------
two_steps = [0 1 0 0 0 1 0 1 0 0 1 0 1 0 0 0 1 0];
three_steps = [0 0 1 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0];
one_step = ~ two_steps & ~ three_steps;

bar_1 = figure;
grid on;
y1 = maximum(stim_two & one_step);
x1 = find(stim_two & one_step);
e1 = err_max(stim_two & one_step);

hold on;
w = 0.5;
bar(x1,y1,w,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.0);
errorbar(x1,y1,e1,'rx');

y2 = abs(minimum(find(stim_one & one_step)));
x2 = find(stim_one & one_step);
e2 = err_min(stim_one & one_step);

bar(x2,y2,w,'FaceColor',[.5 .5 0],'EdgeColor',[.9 .9 0],'LineWidth',1.0);
errorbar(x2,y2,e2,'rx');

legend({'Stim two - Maxima','Error','Stim one - |Minima|'},'Location','northwest');
title('One Step - Signal strength');

ax = gca;
ax.XTick = sort([x1, x2]);
ax.XTickLabels = {'1','3','5','8','14','18'};
ax.XTickLabelRotation = 45;

hold off;

%% 2 Steps --------------------------------------------------------------
bar_2 = figure;
grid on;
y1 = maximum(stim_two & two_steps);
x1 = find(stim_two & two_steps);
e1 = err_max(stim_two &two_steps);

hold on;
w = 0.5;
bar(x1,y1,w,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.0);
errorbar(x1,y1,e1,'rx');

y2 = abs(minimum(find(stim_one & two_steps)));
x2 = find(stim_one & two_steps);
e2 = err_min(stim_one & two_steps);

bar(x2,y2,w,'FaceColor',[.5 .5 0],'EdgeColor',[.9 .9 0],'LineWidth',1.0);
errorbar(x2,y2,e2,'rx');
legend({'Stim two - Maxima','Error','Stim one - |Minima|'},'Location','northwest');
title('Two Steps - Signal strength');

ax = gca;
ax.XTick = sort([x1, x2]);
ax.XTickLabels = {'2','6','8','10','11','13','15'};
ax.XTickLabelRotation = 45;

hold off;
%% 3 Steps --------------------------------------------------------------

bar_3 = figure;
grid on;
y1 = maximum(stim_two & three_steps);
x1 = find(stim_two & three_steps);
e1 = err_max(stim_two & three_steps);

hold on;
w = 0.5;
bar(x1,y1,w,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.0);
errorbar(x1,y1,e1,'rx');

y2 = abs(minimum(find(stim_one & three_steps)));
x2 = find(stim_one & three_steps);
e2 = err_min(stim_one & three_steps);

bar(x2,y2,w,'FaceColor',[.5 .5 0],'EdgeColor',[.9 .9 0],'LineWidth',1.0);
errorbar(x2,y2,e2,'rx');

legend({'Stim two - Maxima','Error','Stim one - |Minima|'},'Location','northwest');
title('Three Steps - Signal strength');

ax = gca;
ax.XTick = sort([x1, x2]);
ax.XTickLabels = {'4','6','12','17'};
ax.XTickLabelRotation = 45;

hold off;
%% Just if needed: Printing - Sorted by FlyID (there will be lots of plots!)
close all;

num_flies = length(unique(grouped_data.flyID));
fly_ID = unique(grouped_data.flyID);

% Do the same like above, just now for each fly
for ii=1:num_flies
    this_fly = fly_ID(ii);
    
    temp = zeros(size(grouped_data.rats));
    counter = 1;
    
    % Select data recorded from this fly
    for jj=1:length(grouped_data.rats)
        
        if(grouped_data.flyIDpN(jj)==this_fly)
            temp(:,:,counter) = crossing_rats(:,:,jj);
            counter = counter +1;
        end
    end
    
    % Mean
    temp = temp(:,:,1:counter);
    temp = mean(temp,3,'omitnan');
    
    % Prepare plots
    plots = length(temp(1,:));
    fly_string = strcat('Fly #',num2str(this_fly),' ');
    
    % Plot for all interesting epoch combinations
    for hh = 1:plots
   
        if abs(sum(mean_val(:,ii))) > THRESHOLD
            fprintf('Signal sum %d: %d \n',ii,abs(sum(mean_val(:,ii))));
            figure(ii*plots + hh);
          %  hold on;
            
            epochRats = squeeze(cur_mat(:,hh,:));
            [x,m,e] = mean_cat_full(epochRats,2,grouped_data.flyIDpN); %m=mean; e=std
            h1 = plot_err_patch_v2(cur_t,m,e,[0 0 1],[0.8 0.8 1]);
            
            xlabel('time (s)');
            ylabel('dF/F - Calcium Signal');
            ylim([-0.4 1]);
            niceaxes;
            title([fly_string,titles(hh)]);
            
            line([0 20],[0 0],'color',[0 0 0]);
            flash = line([10 10],[-0.4 1],'color',[0.5 0.5 1], 'LineWidth',4,'LineStyle',':');
            legend(flash,'Stimulus flash','location','northeast');
          %  hold off;
        else
            fprintf('Case %d did not cross THRESHOLD (%s).\n',ii,titles{ii});
        end
        
    end

end

