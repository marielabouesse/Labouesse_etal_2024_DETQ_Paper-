% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT
% Name: OptoDoseResponse_PooledOverall_2023
% Cohort: PFC dLight1.3b with opto 

% May 2023 - Marie Labouesse, marie.labouesse@gmail.com

% 1- POOLED GRAPHS AND MEASUREMENTS
% loads FP data previously analyzed (with: OptoDoseResponse_DataExtraction then OptoDoseResponse_PooledMatrices) in matlab space "PooledAllMice.mat" from multiple animals (select folder containing multiple animals) 
% you need to input what trials you want to work on
% optionally: calls the script OptoDoseResponse_PooledGraphs_2023 and makes pooled graphs (trace over time, average curves, heatmaps)
% optionally: calls the script OptoDoseResponse_PooledQuantification_2023 calculates AUC, Peak Maxima, DecayHalfTime

% edit relevant parameters in %% SETUP PARAMETERS

% functions needed: OptoDoseResponse_PooledGraphs_2023, OptoDoseResponse_PooledQuantification_2023

% inputs needed: if you want to run graphs, quantifications or both, what trials to work on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATIONS
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% Load the data
% Setup the path
PATH2DATA = uigetdir('select folder'); %path 2 folder above folder of different groups
PATH2SAVEPOOL = PATH2DATA;
mkdir([PATH2SAVEPOOL,'\pooled figures\']);
load([PATH2DATA,'\PooledAllMice.mat']);

% SESSIONS ...............(only 1 session)
s=1;

% PARAMETERS 
show_plot = 1; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
reanalysis = 1; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed.
overwrite = 0; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not          
Collapse_sessions = 0; % If 0, analysis only of individual sessions. If 1, analysis of the group data (pooled)   
done = 0;

% NUM TRIALS
v=1;d=1;k=1;nummice=1;
for pow = 1:length(Opto_Powers) % if multiple powers tested ;not the case in most projects
    numtrials.(Opto_Powers{pow}) = size(PooledINDIV.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}),1);
end

% dFF TO USE FOR QUANTIFICATION
datatype_2use_for_graphs = 2; % 2 if you want to use dFFwithin within 'datatypes'. Check this is the case for each project.

% pooledtype
pooledtype={'raw','baselinecorr'}
animals = AnimalIDs;


%% Determine which trials to analyze and which drug (note: all trials; all drugs area available in the data we just loaded)
% Determine which trials to analyze and which drug (note: all trials; all drugs area available in the data we just loaded)
% N.B. TrialDay = {'VEH','DETQ'}; 
% info about what will be cropped
trialcrop_at_start_VEH = 1; % if we are chopping at the start of VEH: default is start at 1
trialcrop_at_start_DETQ = 0; % if we are chopping at the start of DETQ: default is start at 1
trialcrop_at_end_VEH = 0; % if we are chopping at the end of VEH: default is end at 6
trialcrop_at_end_DETQ = 1; % if we are chopping at the end of DETQ: default is end at 18
% info about trials to crop
TrialDay_2analyze = {'VEH','DETQ'}
if length(TrialDay_2analyze) == 1
    Trials_2analyze.DETQ = [1:14]
elseif length(TrialDay_2analyze) == 2
    Trials_2analyze.VEH = [2:6]
    Trials_2analyze.DETQ = [1:14]    
end
% update time_optotrials_aligned_simple
if length(TrialDay_2analyze) == 1
    time_optotrials_aligned_simple_cropped.(TrialDay_2analyze{v}) = time_optotrials_aligned_simple.(TrialDay_2analyze{v})(Trials_2analyze.DETQ,:)
elseif length(TrialDay_2analyze) == 2
    for v=1:length(TrialDay_2analyze)
        time_optotrials_aligned_simple_cropped.(TrialDay_2analyze{v}) = time_optotrials_aligned_simple.(TrialDay_2analyze{v})(Trials_2analyze.(TrialDay_2analyze{v}),:)
    end
end

% Possible other variables missing
for v=1:length(TrialDay_2analyze)
    mice_list_virus.(TrialDay_2analyze{v}) = mice_list_virus.(TrialDay{v});
end

%% New path for specific trials
if length(TrialDay_2analyze) == 1
    PATH2SAVEPOOL_SELECT = [PATH2DATA,'\SELECTD',' D',num2str(Trials_2analyze.DETQ(1)),'to',num2str(Trials_2analyze.DETQ(end))];    
else
    PATH2SAVEPOOL_SELECT = [PATH2DATA,'\SELECTD', ' V',num2str(Trials_2analyze.VEH(1)),'to',num2str(Trials_2analyze.VEH(end)),...
    ' D',num2str(Trials_2analyze.DETQ(1)),'to',num2str(Trials_2analyze.DETQ(end))];
end
mkdir([PATH2SAVEPOOL_SELECT,'\pooled figures\']);


%% Pool graphs
% WHAT GRAPHS I WANT: 1 if YES, 0 IF NO.
todo_streamsaligned = 0;
todo_averagecurves_allgroups = 0;
todo_indivcurves_onemouse = 0;
todo_heatmap_alltrials_allmice = 0;
todo_heatmap_alltrials_allmice_VEHnormalized = 0;
todo_heatmap_alltrials_onemouse = 0;
todo_heatmap_avetrials_allmice = 1;
todo_vertical_avetrials_allmice = 0;
todo_vertical_alltrials_1mouse = 0;
todo_overlay_alltrials_tempwindow_allmice = 0;
runasindependentscript = 0; % if this variable is defined, will allow to skip the 'setup parameters' section in the OptoDoseResponse_PooledGraphs_2023 script
show_plot = 1;
save_plot = 1;
done = 0;
reanalysis = 1;
overwrite = 1;

% RUN SCRIPT
OptoDoseResponse_PooledGraphs_2023
% Variables i would need if were to make a function instead of a script
% (todo_streamsaligned,todo_averagecurves_allgroups,todo_indivcurves_onemouse,todo_heatmap_alltrials_allmice,...
%      todo_heatmap_alltrials_onemouse,todo_heatmap_avetrials_allmice,todo_vertical_avetrials_allmice,todo_vertical_alltrials_1mouse,Trials_2analyze,TrialDay_2analyze,...
%      mice_list_virus,dFF_names,Opto_Powers,datatype,pooledtype,animals,SessionIDs,PooledINDIV,PooledAVE,PooledINDIV_VEHnormalized,PooledAVE_VEHnormalized,PooledSEM,PooledSEM_VEHnormalized,...
%        TrialDay,AnimalIDs,dt_ds,fullstreams_aligned,show_plot,save_plot,done,reanalysis,overwrite,PATH2SAVEPOOL_SELECT,sampling_rate_ds,datatype_2use_for_graphs,t_trials,TRANGE,runasindependentscript);


%% Run quantifications
% WINDOW DURING WHICH TO CALCULATE AUC AND OTHER VARIABLES: 
windowcalc_duration = 60; %60 seconds
prestim_epoch_start = -60;
prestim_epoch_stop = -1; 
prestim_duration = -(prestim_epoch_start - prestim_epoch_stop); % 9 seconds, ie -10 to -1 sec
data4measurements = 0; % 0 if regular, 1 if normalized;
runasindependentscript = 1;
alreadycroppeddata = 1; % 1 if already ran the cropping of the PooledINDIV etc. based on deciding with trials to analyze using OptoDoseResponse_PooledGraphs_2023; 0 otherwise
offdecay_todo = 0; % 0 if no, 1 if on individual trials, 2 if on the average

% RUN SCRIPT
OptoDoseResponse_PooledQuantification_2023
% FYI results you get are for example: AbsMaxima; AbsMinima; DegChange; DecayHalfTime; Baseline
% for instance: absmax = Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMaxima


%% Extract peak maxima for baseline_corr, dFFwithin 
p=2; d=1; pow=1; k=2; s=1;
color2plot = {'m','r','b','k','c'};

% select the absolute maxima
for v=1:length(TrialDay_2analyze)
    AbsMaxima_all.(TrialDay_2analyze{v}) = ones(length(AnimalIDs),size(Trials_2analyze.(TrialDay_2analyze{v}),2))*nan;
end
for v=1:length(TrialDay_2analyze)
    for nummice=1:length(AnimalIDs)    % all mice
        AbsMaxima_all.(TrialDay_2analyze{v})(nummice,:) = Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMaxima';
    end
end

clear time_trials time_trials_text
% keeping track of time
for v=1:length(TrialDay_2analyze)
    time_trials.(TrialDay_2analyze{v}) = (Trials_2analyze.(TrialDay_2analyze{v})-1)*5;
    time_trials_text.(TrialDay_2analyze{v}) = {};
    for i=1:length(time_trials.(TrialDay_2analyze{v}))
        time_trials_text.(TrialDay_2analyze{v}){i} = num2str(time_trials.(TrialDay_2analyze{v})(i))
    end
end

% plot
figure; clf; hold on
for v=1:length(TrialDay_2analyze)
    tmp_avg_pool = nanmean(AbsMaxima_all.(TrialDay_2analyze{v}),1); 
    tmp_error_pool = nanstd(AbsMaxima_all.(TrialDay_2analyze{v}),1,1)./sqrt(sum(~isnan(AbsMaxima_all.(TrialDay_2analyze{v})(:,1))));
    error_area(time_trials.(TrialDay_2analyze{v}),tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector

    % Formatting
    xlim([1 time_trials.DETQ(end)]);
    ax = gca;
    ax.FontSize = 20; 
    xlabel('Trials (min)','FontSize', 20)
    ylabel('Absolute Maxima','FontSize', 20)
    box off
    ylimits = get(gca,'YLim');
    xlimits = get(gca,'XLim');
    set(gca,'TickLength',[0 0])
end                           
% Legend
L=legend(TrialDay_2analyze{1},' ',TrialDay_2analyze{2},' '); L.Location = 'Best';                
sgtitle([pooledtype{p},' ',datatype{k}],'FontSize',20); % for groups of subplots
%saveplot or not
saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\absmaxima',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.tif']);
saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\absmaxima',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.fig']);
saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\absmaxima',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.pdf']);                  
                            
                            
                            
                            
                            
                            