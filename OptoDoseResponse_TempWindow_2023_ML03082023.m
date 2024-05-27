% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT
% Name: OptoDoseResponse_TempWindow_2023
% Cohort: PFC dLight1.3b with opto 

% May 2023 - Marie Labouesse, marie.labouesse@gmail.com

% 1- CALCULATE TEMP WINDOW
% loads FP data previously analyzed (with: OptoDoseResponse_DataExtraction then OptoDoseResponse_PooledMatrices) in matlab space "PooledAllMice.mat" from multiple animals (select folder containing multiple animals) 
% calls the script OptoDoseResponse_PooledQuantification_2023 calculates AUC, Peak Maxima, DecayHalfTime on whole session (you need to input what trials you want to work on)
% plots Peak Max for all opto-stims over time, finds the maxima, and estimates timepoints when the peak is within 15% of the maxima (or there is an option to use zscore);
% allowing to calculate the temporal window for 'stable' imaging for all
% temporal window number is produced --> calculates it for all + generates data to be plotted in GraphPad

% 2- POOLED GRAPHS AND MEASUREMENTS ON TEMP WINDOW TRIALS
% calls the script OptoDoseResponse_PooledGraphs_2023 and makes pooled graphs (trace over time, average curves, heatmaps) within the temporal window (you need to input what trials you want to work on)
% calls the script OptoDoseResponse_PooledQuantification_2023 calculates AUC, Peak Maxima, DecayHalfTime within the temporal window (you need to input what trials you want to work on)

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

% SESSIONS (only 1 session)
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
datatype_2use_for_graphs = 2; % 2 if you want to use dFFwithin within 'datatypes'

% pooledtype
pooledtype={'raw','baselinecorr'}
animals = AnimalIDs;


%% Determine which trials to analyze and which drug (note: all trials; all drugs area available in the data we just loaded)
% Determine which trials to analyze and which drug (note: all trials; all drugs area available in the data we just loaded)
% note TrialDay = {'VEH','DETQ'}; 
% info about what will be cropped
trialcrop_at_start_VEH = 0; % if we are chopping at the start of VEH: default is start at 1
trialcrop_at_start_DETQ = 0; % if we are chopping at the start of DETQ: default is start at 1
trialcrop_at_end_VEH = 0; % if we are chopping at the end of VEH: default is end at 6
trialcrop_at_end_DETQ = 0; % if we are chopping at the end of DETQ: default is end at 18
% info about trials to crop
TrialDay_2analyze = {'VEH','DETQ'}
if length(TrialDay_2analyze) == 1
    Trials_2analyze.DETQ = [1:18] % here easier to just not edit. so that later the indexes are on the full data and the indexes don't get mixed up
elseif length(TrialDay_2analyze) == 2
    Trials_2analyze.VEH = [1:6] % here easier to just not edit. so that later the indexes are on the full data
    Trials_2analyze.DETQ = [1:18]    % here easier to just not edit. so that later the indexes are on the full data
end
% update time_optotrials_aligned_simple
if length(TrialDay_2analyze) == 1
    time_optotrials_aligned_simple_cropped.(TrialDay_2analyze{v}) = time_optotrials_aligned_simple.(TrialDay_2analyze{v})(Trials_2analyze.DETQ,:)
elseif length(TrialDay_2analyze) == 2
    for v=1:length(TrialDay_2analyze)
        time_optotrials_aligned_simple_cropped.(TrialDay_2analyze{v}) = time_optotrials_aligned_simple.(TrialDay_2analyze{v})(Trials_2analyze.(TrialDay_2analyze{v}),:)
    end
end

% Possible variables missing
for v=1:length(TrialDay_2analyze)
    mice_list_virus.(TrialDay_2analyze{v}) = mice_list_virus.(TrialDay{v});
end

%% New path for specific trials
if length(TrialDay_2analyze) == 1
    PATH2SAVEPOOL_SELECT = [PATH2DATA,'\TEMPWINDOW',' D',num2str(Trials_2analyze.DETQ(1)),'to',num2str(Trials_2analyze.DETQ(end))];    
else
    PATH2SAVEPOOL_SELECT = [PATH2DATA,'\TEMPWINDOW', ' V',num2str(Trials_2analyze.VEH(1)),'to',num2str(Trials_2analyze.VEH(end)),...
    ' D',num2str(Trials_2analyze.DETQ(1)),'to',num2str(Trials_2analyze.DETQ(end))];
end
mkdir([PATH2SAVEPOOL_SELECT,'\pooled figures\']);



%% Calculate temporal window

%% Run quantifications
% WINDOW DURING WHICH TO CALCULATE AUC AND OTHER VARIABLES: 
windowcalc_duration = 10; %10 seconds
prestim_epoch_start = -60;
prestim_epoch_stop = -1; 
prestim_duration = -(prestim_epoch_start - prestim_epoch_stop); % 9 seconds, ie -10 to -1 sec
data4measurements = 1; % 0 if regular, 1 if VEH normalized;
runasindependentscript = 1; % if this variable is defined, will allow to skip the 'setup parameters' section in the OptoDoseResponse_PooledGraphs_2023 script
alreadycroppeddata = 0; % 1 if already ran the cropping of the PooledINDIV etc. based on deciding with trials to analyze using OptoDoseResponse_PooledGraphs_2023; 0 otherwise
offdecay_todo = 0; % 0 if no, 1 if on individual trials, 2 if on the average

% RUN SCRIPT
OptoDoseResponse_PooledQuantification_2023


%% Extract peak maxima over time
p=2; pow=1; s=1; d=1; k=2;  % datatype {'dFF'}    {'dFFwithin'}    {'ZScoredFF_within'}    {'ZScoredFF'}
color2plot = {'k','r','g','b','m','y','c','#7E2F8E'};

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

% plot average
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
    ylim([0 250])
    set(gca,'TickLength',[0 0])
end                           
% Legend
L=legend(TrialDay_2analyze{1},' ',TrialDay_2analyze{2},' '); L.Location = 'Best';                
sgtitle([pooledtype{p},' ',datatype{k}],'FontSize',20); % for groups of subplots
%saveplot or not
saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\absmaxima',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.tif']);
saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\absmaxima',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.fig']);
saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\absmaxima',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.pdf']);                  

figure; clf; hold on
color2plot = {'k','r','g','b','m','y','c','#7E2F8E'};
for v=1:length(TrialDay_2analyze)
    % plot individual trials
    for oi = 1:size(AbsMaxima_all.(TrialDay_2analyze{v}),1)
        plot(time_trials.(TrialDay_2analyze{v}),AbsMaxima_all.(TrialDay_2analyze{v})(oi,:),'Color',color2plot{oi},'LineWidth',v); hold on;
    end
    
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
    ylim([0 350])
end                           
% Legend
% L=legend(TrialDay_2analyze{1},' ',TrialDay_2analyze{2},' '); L.Location = 'Best';                
sgtitle([pooledtype{p},' ',datatype{k}],'FontSize',20); % for groups of subplots
%saveplot or not
saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\absmaxima indiv',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.tif']);
saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\absmaxima indiv',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.fig']);
saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\absmaxima indiv',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.pdf']);                  

%% Determine the maxima and the trials within temp windows, then the temp window of stability
% Determine the maximal peak maxima, average and average value at the plateau (5-20min). This is to calculate the plateau.
timepoint5min = find(time_optotrials_aligned_simple_cropped.DETQ(:,1) == 5);
timepoint10min = find(time_optotrials_aligned_simple_cropped.DETQ(:,1) == 10);
timepoint15min = find(time_optotrials_aligned_simple_cropped.DETQ(:,1) == 15);
timepoint20min = find(time_optotrials_aligned_simple_cropped.DETQ(:,1) == 20);
for v=1:length(TrialDay_2analyze)
    for nummice=1:length(AnimalIDs)    % all mice
        MAX_AbsMaxima_all.(TrialDay_2analyze{v})(nummice) = max(AbsMaxima_all.(TrialDay_2analyze{v})(nummice,:),[],2);
        AVE_AbsMaxima_all.(TrialDay_2analyze{v})(nummice) = mean(AbsMaxima_all.(TrialDay_2analyze{v})(nummice,:));
        EMPIRIMAX_AbsMaxima_all.(TrialDay_2analyze{v})(nummice) = mean([AbsMaxima_all.(TrialDay_2analyze{v})(nummice,timepoint5min),AbsMaxima_all.(TrialDay_2analyze{v})(nummice,timepoint10min),...
                                                                AbsMaxima_all.(TrialDay_2analyze{v})(nummice,timepoint15min),AbsMaxima_all.(TrialDay_2analyze{v})(nummice,timepoint20min)]);
        Zscore.(TrialDay_2analyze{v})(nummice,:) = zscore(AbsMaxima_all.(TrialDay_2analyze{v})(nummice,:)); % in case you want to use the zscore to estimate the stable temp window
    end
end

% Determine the trials within 15% of the 10 min timepoint (aternatively, calculate the zscore for the data, then take the top 33% of the data
for v=1:length(TrialDay_2analyze)
    WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v}) = ones(length(AnimalIDs),size(AbsMaxima_all.(TrialDay_2analyze{v})(nummice,:),2))*nan;
end
for v=1:length(TrialDay_2analyze)
    for nummice=1:length(AnimalIDs)    % all mice
        for ju = 1:length(AbsMaxima_all.(TrialDay_2analyze{v})(nummice,:))
            if AbsMaxima_all.(TrialDay_2analyze{v})(nummice,ju) > 0.85 * EMPIRIMAX_AbsMaxima_all.(TrialDay_2analyze{v})(nummice) && AbsMaxima_all.(TrialDay_2analyze{v})(nummice,ju) < 1.15 * EMPIRIMAX_AbsMaxima_all.(TrialDay_2analyze{v})(nummice)% version with %
%             if Zscore.(TrialDay_2analyze{v})(nummice,ju) > 0.44 % version with zscore. >0.44 gts the top 33%
                WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,ju) = 1;
            else
                WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,ju) = 0;
            end
        end
    end
end
WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})

% Find indexes (written a bit redundant below, could simplify the code)
for v=2; %1:length(TrialDay_2analyze)
    WINSTABL_IX.(TrialDay_2analyze{v}) = ones(length(AnimalIDs),2)*nan;
end
for v=2; %1:length(TrialDay_2analyze)
    for nummice=1:length(AnimalIDs)    % all mice
        % first ix, its going to be either 5 or 10 min (not before 5 min: the kinetics are too fast, it would not be stable enough or reliable)
        if WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,2) == 1
            WINSTABL_IX.(TrialDay_2analyze{v})(nummice,1) = 2; % 5 min timepoint
        else
            WINSTABL_IX.(TrialDay_2analyze{v})(nummice,1) = 3; % 10 min timepoint
        end
        % second ix, its going to be anyting as of 10 min        
        if  WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,4) == 0
            WINSTABL_IX.(TrialDay_2analyze{v})(nummice,2) = 3; % 10 min timepoint
        elseif  WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,5) == 0
            WINSTABL_IX.(TrialDay_2analyze{v})(nummice,2) = 4; % 15 min timepoint
        elseif  WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,6) == 0
            WINSTABL_IX.(TrialDay_2analyze{v})(nummice,2) = 5; % 20 min timepoint
        elseif  WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,7) == 0
            WINSTABL_IX.(TrialDay_2analyze{v})(nummice,2) = 6; % 25 min timepoint
        elseif  WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,8) == 0
            WINSTABL_IX.(TrialDay_2analyze{v})(nummice,2) = 7; % 30 min timepoint            
        elseif  WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,9) == 0
            WINSTABL_IX.(TrialDay_2analyze{v})(nummice,2) = 8; % 35 min timepoint            
        elseif  WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,10) == 0
            WINSTABL_IX.(TrialDay_2analyze{v})(nummice,2) = 9; % 40 min timepoint
        elseif  WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,11) == 0
            WINSTABL_IX.(TrialDay_2analyze{v})(nummice,2) = 10; % 45 min timepoint
        elseif find(WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,:))
            WINSTABL_IX.(TrialDay_2analyze{v})(nummice,2) = 15; % beyond 10 up until last timepoint            
        end
    end
end
All_Indiv_TempWindowStab = WINSTABL_IX.(TrialDay_2analyze{v})

% Average
Average_start_TempWindowStab = mean(All_Indiv_TempWindowStab(:,1))
Average_end_TempWindowStab = mean(All_Indiv_TempWindowStab(:,2))

%% Determine the trials where performance is substantially superior, at least 150% of VEH; was not implemented for now
run_window_supperf = 0;
if run_window_supperf == 1;
% Determine the trials where DETQ > ave VEH by at least 50%
WINSUPPERF_AbsMaxima_all.DETQ = ones(length(AnimalIDs),size(AbsMaxima_all.(TrialDay_2analyze{v})(nummice,:),2))*nan;
for nummice=1:length(AnimalIDs)    % all mice
    for ju = 1:length(AbsMaxima_all.(TrialDay_2analyze{v})(nummice,:))
        if AbsMaxima_all.DETQ(nummice,ju) > 1.5 * AVE_AbsMaxima_all.VEH(nummice)
            WINSUPPERF_AbsMaxima_all.DETQ(nummice,ju) = 1;
        else
            WINSUPPERF_AbsMaxima_all.DETQ(nummice,ju) = 0;
        end
    end
end
WINSUPPERF_AbsMaxima_all.DETQ

% Find indexes
WINSUPPERF_IX.DETQ = ones(length(AnimalIDs),2)*nan;
for nummice=1:length(AnimalIDs)    % all mice
    % first ix, its going to be either 5 or 10 min
    if WINSUPPERF_AbsMaxima_all.DETQ(nummice,1) == 1
        WINSUPPERF_IX.DETQ(nummice,1) = 2; % 5 min timepoint
    else
        WINSUPPERF_IX.DETQ(nummice,1) = 3; % 10 min timepoint
    end
    % second ix, its going to be anyting as of 10 min        
    if  WINSTABL_AbsMaxima_all.(TrialDay_2analyze{v})(nummice,4) == 0
        WINSUPPERF_IX.DETQ(nummice,2) = 3; % 10 min timepoint
    elseif  WINSTABL_AbsMaxima_all.DETQ(nummice,5) == 0
        WINSUPPERF_IX.DETQ(nummice,2) = 4; % 15 min timepoint
    elseif  WINSTABL_AbsMaxima_all.DETQ(nummice,6) == 0
        WINSUPPERF_IX.DETQ(nummice,2) = 5; % 20 min timepoint
    elseif  WINSTABL_AbsMaxima_all.DETQ(nummice,7) == 0
        WINSUPPERF_IX.DETQ(nummice,2) = 6; % 25 min timepoint
    elseif  WINSTABL_AbsMaxima_all.DETQ(nummice,8) == 0
        WINSUPPERF_IX.DETQ(nummice,2) = 7; % 30 min timepoint            
    elseif  WINSTABL_AbsMaxima_all.DETQ(nummice,9) == 0
        WINSUPPERF_IX.DETQ(nummice,2) = 8; % 35 min timepoint            
    elseif  WINSTABL_AbsMaxima_all.DETQ(nummice,10) == 0
        WINSUPPERF_IX.DETQ(nummice,2) = 9; % 40 min timepoint
    elseif  WINSTABL_AbsMaxima_all.DETQ(nummice,11) == 0
        WINSUPPERF_IX.DETQ(nummice,2) = 10; % 45 min timepoint
    elseif  WINSTABL_AbsMaxima_all.DETQ(nummice,12) == 0
        WINSUPPERF_IX.DETQ(nummice,2) = 11; % 50 min timepoint
    elseif  WINSTABL_AbsMaxima_all.DETQ(nummice,13) == 0
        WINSUPPERF_IX.DETQ(nummice,2) = 12; % 55 min timepoint
    elseif  WINSTABL_AbsMaxima_all.DETQ(nummice,14) == 0
        WINSUPPERF_IX.DETQ(nummice,2) = 13; % 60 min timepoint
    elseif  WINSTABL_AbsMaxima_all.DETQ(nummice,15) == 0
        WINSUPPERF_IX.DETQ(nummice,2) = 14; % 65 min timepoint
    elseif find(WINSTABL_AbsMaxima_all.DETQ(nummice,:))
        WINSUPPERF_IX.DETQ(nummice,2) = 15; % last timepoint            
    end
end
WINSUPPERF_IX.DETQ

end

%% Polyfit; was not implemented for now
runpolyfit = 0
if runpolyfit == 1
    % polyfit of the absolute maxima only as of trial 10 min
    timepoint10min = find(time_optotrials_aligned_simple_cropped.DETQ(:,1) == 10);
    for v=2; % 1:length(TrialDay_2analyze)
        PolyfitAbsMaxima_all.(TrialDay_2analyze{v}) = ones(length(AnimalIDs),size(Trials_2analyze.(TrialDay_2analyze{v}),2)-1)*nan;
    end
    for v=2; %1:length(TrialDay_2analyze)
        figure; 
        for nummice=1:length(AnimalIDs)    % all mice       
            % Step 1: Polyfit
            x_poly = time_trials.(TrialDay_2analyze{v})(2:end);
            to_fit = AbsMaxima_all.(TrialDay_2analyze{v})(nummice,2:end);
            p_fit = polyfit(x_poly,to_fit,1); %fitting control Dv2 against signal Dv1 using least squares fit of degree 1
            f_fit = polyval(p_fit,x_poly);
            subplot(1,length(AnimalIDs),nummice); plot(x_poly,to_fit); hold on; plot(x_poly,f_fit); 
            a_fitted = to_fit(2)+(to_fit-f_fit)./f_fit;
            subplot(1,length(AnimalIDs),nummice); plot(x_poly,a_fitted); hold on; xlabel('Time (min'); ylabel(['Mouse ',num2str(nummice),' (',AnimalIDs{nummice},')']); ylim([0 15])
            PolyfitAbsMaxima_all.(TrialDay_2analyze{v})(nummice,:) = a_fitted;
        end
        L=legend('AbsMaxima','polyfit','fitted');
    end
    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\absmaxima_polyfit',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.tif']);
    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\absmaxima_polyfit',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.fig']);
    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\absmaxima_polyfit',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.pdf']);                  

end


%% %% %% %% WINDOW OF STABILITY

%% New path for specific trials
if length(TrialDay_2analyze) == 1
    PATH2SAVEPOOL_SELECT = [PATH2DATA,'\TEMPWINDOWSTABILITY',' D',num2str(Trials_2analyze.DETQ(1)),'to',num2str(Trials_2analyze.DETQ(end))];    
else
    PATH2SAVEPOOL_SELECT = [PATH2DATA,'\TEMPWINDOWSTABILITY', ' V',num2str(Trials_2analyze.VEH(1)),'to',num2str(Trials_2analyze.VEH(end)),...
    ' D',num2str(Trials_2analyze.DETQ(1)),'to',num2str(Trials_2analyze.DETQ(end))];
end
mkdir([PATH2SAVEPOOL_SELECT,'\pooled figures\']);


%% Window of stability
for v=2; %:length(TrialDay_2analyze)
    % Create a new variable with the temporal window of stability for each animal where DETQ responses are within 15% of the max
    tempwindow_stability_ix = ones(length(mice_list_virus.(TrialDay_2analyze{v})),2)*nan;
    for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
        tempwindow_stability_ix(nummice,1) = WINSTABL_IX.(TrialDay_2analyze{v})(nummice,1); % replace 
        tempwindow_stability_ix(nummice,2) = WINSTABL_IX.(TrialDay_2analyze{v})(nummice,2); % replace
    end

    % Plot the temporal window
        color2plot = {'k','r','g','b','m','y','c','#7E2F8E'};        
        figure; clf; hold on
        for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
            plot([tempwindow_stability_ix(nummice,1) tempwindow_stability_ix(nummice,2)],[nummice nummice],'Color',color2plot{nummice}); %t_trials is the time vector
        end
        % Formatting
        ylim([0 length(mice_list_virus.(TrialDay_2analyze{v}))+1]);
        xlim([1 16]);
        ax = gca;
        ax.FontSize = 20; 
        ylabel('Mice','FontSize', 20)
        xlabel('Temp Window (trials)','FontSize', 20)
        box off
        ylimits = get(gca,'YLim');
        xlimits = get(gca,'XLim');
        xticks([2 4 6 8 10 12 14]);
        yticks([1 2 3 4 5 6 7]);        
        set(gca,'TickLength',[0 0])
        sgtitle([TrialDay_2analyze{v},' TempWindow Stability'],'FontSize',20); % for groups of subplots
        %saveplot or not
        saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\tempwindowstability',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.tif']);
        saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\tempwindowstability',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.fig']);
        saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\tempwindowstability',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.pdf']);                  
end

% Get Average temp window % will be for DETQ. Uses closest integer
tempwindow_stability_ave_ix1 = round(nanmean(tempwindow_stability_ix(:,1),1))
tempwindow_stability_ave_ix2 = round(nanmean(tempwindow_stability_ix(:,2),1))
tempwindow_stability_ave_allmice = [tempwindow_stability_ave_ix1:tempwindow_stability_ave_ix2]


%% Determine which trials to now get graphs and quantif from; here specifically on the temp window trials
trialcrop_at_start_VEH = 1; % if we are chopping at the start of VEH: default is start at 1
trialcrop_at_start_DETQ = 1; % if we are chopping at the start of DETQ: default is start at 1
trialcrop_at_end_VEH = 0; % if we are chopping at the end of VEH: default is end at 6
trialcrop_at_end_DETQ = 1; % if we are chopping at the end of DETQ: default is end at 18
% info about trials to crop
TrialDay_2analyze = {'VEH','DETQ'}
if length(TrialDay_2analyze) == 1
    Trials_2analyze.DETQ = tempwindow_stability_ave_allmice; % tempwindow_stability_ave_allmice or tempwindow_supperf_ave_allmice
elseif length(TrialDay_2analyze) == 2
    Trials_2analyze.VEH = [2:6]
    Trials_2analyze.DETQ = tempwindow_stability_ave_allmice; % tempwindow_stability_ave_allmice or tempwindow_supperf_ave_allmice
end
% update time_optotrials_aligned_simple
if length(TrialDay_2analyze) == 1
    time_optotrials_aligned_simple_cropped.(TrialDay_2analyze{v}) = time_optotrials_aligned_simple.(TrialDay_2analyze{v})(Trials_2analyze.DETQ,:)
elseif length(TrialDay_2analyze) == 2
    for v=1:length(TrialDay_2analyze)
        time_optotrials_aligned_simple_cropped.(TrialDay_2analyze{v}) = time_optotrials_aligned_simple.(TrialDay_2analyze{v})(Trials_2analyze.(TrialDay_2analyze{v}),:)
    end
end


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
runasindependentscript = 1; % if this variable is defined, will allow to skip the 'setup parameters' section in the OptoDoseResponse_PooledGraphs_2023 script
show_plot = 1;
save_plot = 1;
done = 0;
reanalysis = 1;
overwrite = 1;

% RUN SCRIPT
OptoDoseResponse_PooledGraphs_2023


%% Run quantifications: needs to be run after OptoDoseResponse_PooledGraphs where cropping is happening
% WINDOW DURING WHICH TO CALCULATE AUC AND OTHER VARIABLES: 
windowcalc_duration = 60; %60 seconds
prestim_epoch_start = -60; %-60 seconds
prestim_epoch_stop = -1; %-1 seconds
prestim_duration = -(prestim_epoch_start - prestim_epoch_stop); % 9 seconds, ie -10 to -1 sec
data4measurements = 1; % 0 if regular, 1 if normalized;
runasindependentscript = 1;
alreadycroppeddata = 1; % 1 if already ran the cropping of the PooledINDIV etc. based on deciding with trials to analyze using OptoDoseResponse_PooledGraphs_2023; 0 otherwise
offdecay_todo = 0; % 0 if no, 1 if on individual trials, 2 if on the average

% RUN SCRIPT
OptoDoseResponse_PooledQuantification_2023


%%
if run_window_supperf == 1;

%% %% %% %% WINDOW OF SUPERIOR PERFORMANCE

%% New path for specific trials
if length(TrialDay_2analyze) == 1
    PATH2SAVEPOOL_SELECT = [PATH2DATA,'\TEMPWINDOWSUPPERF',' D',num2str(Trials_2analyze.DETQ(1)),'to',num2str(Trials_2analyze.DETQ(end))];    
else
    PATH2SAVEPOOL_SELECT = [PATH2DATA,'\TEMPWINDOWSUPPERF', ' V',num2str(Trials_2analyze.VEH(1)),'to',num2str(Trials_2analyze.VEH(end)),...
    ' D',num2str(Trials_2analyze.DETQ(1)),'to',num2str(Trials_2analyze.DETQ(end))];
end
mkdir([PATH2SAVEPOOL_SELECT,'\pooled TEMPWINDOWSUPPERF\']);


%% Window of superior performance
for v=1:length(TrialDay_2analyze)
    % Create a new variable with the temporal window of superior performance for each animal where DETQ > VEH by at least 50% 
    tempwindow_supperf_ix = ones(length(mice_list_virus.(TrialDay_2analyze{v})),2)*nan;
    for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
        tempwindow_supperf_ix(nummice,1) = WINSUPPERF_IX.DETQ(nummice,1); % replace 
        tempwindow_supperf_ix(nummice,2) = WINSUPPERF_IX.DETQ(nummice,2); % replace
    end

    % Plot the temporal window
    color2plot = {'m','r','b','k','c','y','g'};
    figure; clf; hold on
    for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
        plot([tempwindow_stability_ix(nummice,1) tempwindow_stability_ix(nummice,2)],[nummice nummice],'Color',color2plot{nummice}); %t_trials is the time vector
    end
    % Formatting
    ylim([0 length(mice_list_virus.(TrialDay_2analyze{v}))+1]);
    xlim([1 16]);
    ax = gca;
    ax.FontSize = 20; 
    ylabel('Mice','FontSize', 20)
    xlabel('Temp Window (trials)','FontSize', 20)
    box off
    ylimits = get(gca,'YLim');
    xlimits = get(gca,'XLim');
    xticks([2 4 6 8 10 12 14]);
    yticks([1 2 3 4 5 6 7]);        
    set(gca,'TickLength',[0 0])
    sgtitle([TrialDay_2analyze{v},' TempWindow Stability'],'FontSize',20); % for groups of subplots
    %saveplot or not
    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\tempwindowstability',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.tif']);
    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\tempwindowstability',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.fig']);
    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\tempwindowstability',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.pdf']);                  
end

% Get Average temp window sup perf
tempwindow_supperf_ave_ix1 = round(nanmean(tempwindow_supperf_ix(:,1),1))
tempwindow_supperf_ave_ix2 = round(nanmean(tempwindow_supperf_ix(:,2),1))
tempwindow_supperf_ave_allmice = [tempwindow_supperf_ave_ix1:tempwindow_supperf_ave_ix2]


%% Determine which trials to now get graphs and quantif from; here specifically on the temp window trials
trialcrop_at_start_VEH = 0; % if we are chopping at the start of VEH: default is start at 1
trialcrop_at_start_DETQ = 1; % if we are chopping at the start of DETQ: default is start at 1
trialcrop_at_end_VEH = 0; % if we are chopping at the end of VEH: default is end at 6
trialcrop_at_end_DETQ = 1; % if we are chopping at the end of DETQ: default is end at 18
% info about trials to crop
TrialDay_2analyze = {'VEH','DETQ'}
if length(TrialDay_2analyze) == 1
    Trials_2analyze.DETQ = tempwindow_supperf_ave_allmice; % tempwindow_stability_ave_allmice or tempwindow_supperf_ave_allmice
elseif length(TrialDay_2analyze) == 2
    Trials_2analyze.VEH = [2:6]
    Trials_2analyze.DETQ = tempwindow_supperf_ave_allmice; % tempwindow_stability_ave_allmice or tempwindow_supperf_ave_allmice
end
% update time_optotrials_aligned_simple
if length(TrialDay_2analyze) == 1
    time_optotrials_aligned_simple_cropped.(TrialDay_2analyze{v}) = time_optotrials_aligned_simple.(TrialDay_2analyze{v})(Trials_2analyze.DETQ,:)
elseif length(TrialDay_2analyze) == 2
    for v=1:length(TrialDay_2analyze)
        time_optotrials_aligned_simple_cropped.(TrialDay_2analyze{v}) = time_optotrials_aligned_simple.(TrialDay_2analyze{v})(Trials_2analyze.(TrialDay_2analyze{v}),:)
    end
end


%% Pool graphs
% WHAT GRAPHS I WANT: 1 if YES, 0 IF NO.
todo_streamsaligned = 1;
todo_averagecurves_allgroups = 1;
todo_indivcurves_onemouse = 1;
todo_heatmap_alltrials_allmice = 1;
todo_heatmap_alltrials_onemouse = 1;
todo_heatmap_avetrials_allmice = 1;
todo_vertical_avetrials_allmice = 1;
todo_vertical_alltrials_1mouse = 1;
todo_overlay_alltrials_tempwindow_allmice = 1;
runasindependentscript = 1;
show_plot = 1;
save_plot = 1;
done = 0;
reanalysis = 1;
overwrite = 1;

% RUN SCRIPT
OptoDoseResponse_PooledGraphs_2023


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

end
