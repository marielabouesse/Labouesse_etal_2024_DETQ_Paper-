% OptoDoseResponse_PooledGraphs_2023

% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT
% Name: OptoDoseResponse_PooledGraphs_2023
% Cohort: PFC dLight1.3b with opto 

% May 2023 - Marie Labouesse, marie.labouesse@gmail.com

% 1- POOLED GRAPHS
% loads FP data previously analyzed (with: OptoDoseResponse_DataExtraction then OptoDoseResponse_PooledMatrices) in matlab space "PooledAllMice.mat" from multiple animals (select folder containing multiple animals) 
% or if used as an independent script the data was loaded previously within OptoDoseResponse_PooledOverall_2023 or OptoDoseResponse_TempWindow_2023
% generates average graphs curves 
% generates heatmaps (1 average trace/mouse or all trials/all mice in one big heatmap - separated by trial type)
% generates graphs with trials over time 

% edit relevant parameters in %% SETUP PARAMETERS or within the code running above this one

% functions needed: error_area_onlyrectangle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('runasindependentscript') % then we don't need the below as the other scripts already initialized everything we need
    %% INITIALIZATIONS  --> run this if you want to use this code independently.
    close all; clear variables; clc;
    set(0,'defaultfigurecolor',[1 1 1])

    % SETUP PARAMETERS
    PATH2DATA = uigetdir('select folder'); %select the overall folder, 
    PATH2SAVEPOOL = PATH2DATA;
    mkdir([PATH2SAVEPOOL,'\pooled figures\']);

    % load data
    load([PATH2DATA,'\PooledAllMice.mat']);

    % dFF TO USE FOR QUANTIFICATION
    datatype_2use_for_graphs = 2; % 2 if you want to use dFFwithin within 'datatypes' as the main datatype to plot; 1 for dFF, 4 for Zscore (3 for Zscore dFF within which is not in use)

    % PARAMETERS 
    show_plot = 1; % If 0, plots are not displayed
    save_plot = 1; % If 0, plots are not saved
    reanalysis = 1; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed. 
    overwrite = 0; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not          
    Collapse_sessions = 0; % If 0, analysis only of individual sessions. If 1, analysis of the group data (pooled)   
    done = 0;

    % pooledtype
    pooledtype={'raw','baselinecorr'}
    animals = AnimalIDs;

    % WHAT TRIALS TO ANALYZE FOR THE AVERAGE CURVE OR THE HEATMAP CURVE; OR WHAT TO SHOW FOR THE OTHER GRAPHS WITH INDIVIDUAL DATA
    % N.B. TrialDay = {'VEH','DETQ'}; 
    TrialDay_2analyze = {'VEH','DETQ'}
    if length(TrialDay_2analyze) == 1
        Trials_2analyze.DETQ = [2:7]
    elseif length(TrialDay_2analyze) == 2
        Trials_2analyze.VEH = [2:4]
        Trials_2analyze.DETQ = [2:7]    
    end
    % update time_optotrials_aligned_simple
    if length(TrialDay_2analyze) == 1
        time_optotrials_aligned_simple_cropped.(TrialDay_2analyze{v}) = time_optotrials_aligned_simple.(TrialDay_2analyze{v})(Trials_2analyze.DETQ,:)
    elseif length(TrialDay_2analyze) == 2
        for v=1:length(TrialDay_2analyze)
            time_optotrials_aligned_simple_cropped.(TrialDay_2analyze{v}) = time_optotrials_aligned_simple.(TrialDay_2analyze{v})(Trials_2analyze.(TrialDay_2analyze{v}),:)
        end
    end

    % New path for specific trials
    if length(TrialDay_2analyze) == 1
        PATH2SAVEPOOL_SELECT = [PATH2DATA,'\SELECTD',' D',num2str(Trials_2analyze.DETQ(1)),'to',num2str(Trials_2analyze.DETQ(end))];    
    else
        PATH2SAVEPOOL_SELECT = [PATH2DATA,'\SELECTD', ' V',num2str(Trials_2analyze.VEH(1)),'to',num2str(Trials_2analyze.VEH(end)),...
        ' D',num2str(Trials_2analyze.DETQ(1)),'to',num2str(Trials_2analyze.DETQ(end))];
    end
    mkdir([PATH2SAVEPOOL_SELECT,'\pooled figures\']);

    % Possible variables missing
    for v=1:length(TrialDay_2analyze)
        mice_list_virus.(TrialDay_2analyze{v}) = mice_list_virus.(TrialDay{v});
    end
    
    % which graphs to I want to plot
    todo_streamsaligned = 1;
    todo_averagecurves_allgroups = 1;
    todo_indivcurves_onemouse = 1;
    todo_heatmap_alltrials_allmice = 1;
    todo_heatmap_alltrials_allmice_VEHnormalized = 0;
    todo_heatmap_alltrials_onemouse = 1;
    todo_heatmap_avetrials_allmice = 1;
    todo_vertical_avetrials_allmice = 1;
    todo_vertical_alltrials_1mouse = 1;
    todo_overlay_alltrials_tempwindow_allmice = 1;

end

% ABOVE THIS LINE NOT IN USE WHEN USED AS A INDEPENDENT SCRIPT CALLED OUT FROM OTHER SCRIPTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SESSIONS - only 1 session
s=1;

%% Info about what will be cropped; we only want 5 trials for VEH and 14 for DETQ (or less if we are focusing on defining the temporal window) 
if length(TrialDay_2analyze) == 1
    if Trials_2analyze.DETQ(1) == 1; trialcrop_at_start_DETQ = 0; else trialcrop_at_start_DETQ = 1; end % if we are chopping at the start of DETQ: default is start at 1
    if Trials_2analyze.DETQ(end) == 18; trialcrop_at_end_DETQ = 0; else trialcrop_at_end_DETQ = 1; end % if we are chopping at the end of DETQ: default is end at 18
else
    if Trials_2analyze.VEH(1) == 1; trialcrop_at_start_VEH = 0; else trialcrop_at_start_VEH = 1; end % if we are chopping at the start of VEH: default is start at 1
    if Trials_2analyze.VEH(end) == 6; trialcrop_at_end_VEH = 0; else trialcrop_at_end_VEH = 1; end % if we are chopping at the end of VEH: default is end at 6

    if Trials_2analyze.DETQ(1) == 1; trialcrop_at_start_DETQ = 0; else trialcrop_at_start_DETQ = 1; end % if we are chopping at the start of DETQ: default is start at 1
    if Trials_2analyze.DETQ(end) == 18; trialcrop_at_end_DETQ = 0; else trialcrop_at_end_DETQ = 1; end % if we are chopping at the end of DETQ: default is end at 18
end

%% CROP THE DATA TRIALS NOT TO ANALYZE IN THIS RUN
if trialcrop_at_start_VEH == 1 || trialcrop_at_start_DETQ == 1 || trialcrop_at_end_VEH == 1 || trialcrop_at_end_DETQ == 1
    % crop PooledINDIV
    for v=1:length(TrialDay_2analyze)
        trials_2_keep = [Trials_2analyze.(TrialDay_2analyze{v})];
        for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
            for d=1:length(dFF_names)
                for pow = 1:length(Opto_Powers)
                    for k = 1:size(datatype,2)
                        for p = 1:length(pooledtype)
                             % PooledINDIV
                             PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s}) = ...
                                 PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s})(trials_2_keep,:);
                             % PooledINDIV_VEHnormalized
                             PooledINDIV_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s}) = ...
                                 PooledINDIV_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s})(trials_2_keep,:);
                        end
                    end
                end
            end
        end
    end

    % Calculate averages again
    clear PooledAVE PooledSEM PooledAVE_VEHnormalized PooledSEM_VEHnormalized
    for v=1:length(TrialDay)
        for p=1:length(pooledtype)
            for d=1:length(dFF_names)
                for pow = 1:length(Opto_Powers)
                    for k = 1:size(datatype,2)
                        for nummice=1:length(mice_list_virus.(TrialDay{v}))
                            for s = 1:length(SessionIDs)
                                data2average = PooledINDIV.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s});
                                PooledAVE.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = ...
                                    nanmean(data2average,1);
                                PooledSEM.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = ...
                                    nanstd(data2average,1,1)./sqrt(size(data2average,1));
                                clear data2average
                                data2average = PooledINDIV_VEHnormalized.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s});
                                PooledAVE_VEHnormalized.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = ...
                                    nanmean(data2average,1);
                                PooledSEM_VEHnormalized.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = ...
                                    nanstd(data2average,1,1)./sqrt(size(data2average,1));                            
                           end
                        end
                    end
                end
            end
        end
    end
end


%% GRAPHS
%% Streams aligned on a graph; ALL DRUGS TOGETHER ON SAME GRAPH, all animals
% todo_streamsaligned = 1
if todo_streamsaligned == 1
    % Chop off the last 30 to 60 sec (potential artefacts)
    chop_off_ix.VEH = round(60./dt_ds); % 60 seconds
    chop_off_ix.DETQ = chop_off_ix.VEH;
    
    % Also chop trials at the end that are not selected for this analysis (see above which trials we need)
    if trialcrop_at_end_VEH == 1
        numtrials_to_chop_at_end_VEH = 6 - Trials_2analyze.VEH(end);
        chop_off_ix.VEH = chop_off_ix.VEH + (numtrials_to_chop_at_end_VEH * 5 * round(60./dt_ds)); % numtrials_to_chop_at_end * 5 min * round(60./dt_ds) to be in seconds and then convert to indexes
    end
    if trialcrop_at_end_DETQ == 1
        numtrials_to_chop_at_end_DETQ = 18 - Trials_2analyze.DETQ(end);        
        chop_off_ix.DETQ = chop_off_ix.DETQ + (numtrials_to_chop_at_end_DETQ * 5 * round(60./dt_ds)); % numtrials_to_chop_at_end * 5 min * round(60./dt_ds) to be in seconds and then convert to indexes
    end
    
    % Also chop trials at the beginning that are not selected for this analysis
    if trialcrop_at_start_VEH == 1
        numtrials_to_chop_at_start_VEH = Trials_2analyze.VEH(1) -1;
        chop_start_ix.VEH = (numtrials_to_chop_at_start_VEH * 5 * round(60./dt_ds)); % numtrials_to_chop_at_end * 5 min * round(60./dt_ds) to be in seconds and then convert to indexes
    end
%     if trialcrop_at_start_DETQ == 1
        numtrials_to_chop_at_start_DETQ = Trials_2analyze.DETQ(1) -1;        
        chop_start_ix.DETQ = (numtrials_to_chop_at_start_DETQ * 5 * round(60./dt_ds)); % numtrials_to_chop_at_end * 5 min * round(60./dt_ds) to be in seconds and then convert to indexes
%     end
    
    % create empty pooled array
    nummice = 1; k=2; d=1; pow = 1; s=1; 
    for v=1:length(TrialDay_2analyze)    
        for d=1:length(dFF_names)
            for pow = 1:length(Opto_Powers)
                for k = 1:size(datatype,2)
                    length_fullstreams_aligned.(TrialDay_2analyze{v}) = ...
                    size(find(~isnan(fullstreams_aligned.(TrialDay_2analyze{v}).(AnimalIDs{nummice}).(SessionIDs{s}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}))),2)       

                    fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}) = ...
                        ones(length(mice_list_virus.(TrialDay_2analyze{v})),length_fullstreams_aligned.(TrialDay_2analyze{v})-chop_off_ix.(TrialDay_2analyze{v})+1)*nan;
                    length_updated.(TrialDay_2analyze{v}) = length(fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}));
                end
            end
        end

        %STORE STREAMS ALIGNED TO OPTO ONSET #1 with the average of each mouse in each row (eg 3 mice --> 3 rows)
        for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
            for d=1:length(dFF_names)
                for pow = 1:length(Opto_Powers)
                    for k = 1:size(datatype,2)
                        fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(nummice,1:length_updated.(TrialDay_2analyze{v})) = ...
                            fullstreams_aligned.(TrialDay_2analyze{v}).(AnimalIDs{nummice}).(SessionIDs{s}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(1:length_updated.(TrialDay_2analyze{v}));
                        
                            
                    end
                end
            end
        end
    end

    % Make new time vectors and also delete the NaNs in the VEH
    for v=1:length(TrialDay_2analyze)
        for d=1:length(dFF_names)
            for pow = 1:length(Opto_Powers)
                for k = 1:size(datatype,2)
                     length_data_aligned{v} = size(find(~isnan(fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(1,:))),2);
                     time_vect_aligned_special.(TrialDay_2analyze{v}) = 0:dt_ds*1:dt_ds*(length_data_aligned{v}-1);
                     fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(:,length_data_aligned{v}+1:end) = [] ;                   
                end
            end
        end
    end
    
    % Now chop off at the beginning any trials not selected for this analysis
    if trialcrop_at_start_VEH == 1 || trialcrop_at_start_DETQ == 1
        % create empty pooled array
        for v=1:length(TrialDay_2analyze)    
            for d=1:length(dFF_names)
                for pow = 1:length(Opto_Powers)
                    for k = 1:size(datatype,2)
                        length_fullstreams_aligned_pooled.(TrialDay_2analyze{v}) = ...
                            size(fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),2)           
                        fullstreams_aligned_pooled_v2.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}) = ...
                            ones(length(mice_list_virus.(TrialDay_2analyze{v})),length_fullstreams_aligned_pooled.(TrialDay_2analyze{v})-chop_start_ix.(TrialDay_2analyze{v})+1)*nan;
                        length_updated_v2.(TrialDay_2analyze{v}) = length(fullstreams_aligned_pooled_v2.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}))
                    end
                end
            end
        end
        % fill
        for v=1:length(TrialDay_2analyze)            
            for d=1:length(dFF_names)
                for pow = 1:length(Opto_Powers)
                    for k = 1:size(datatype,2)
                        data2incl = fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})...
                            (:,chop_start_ix.(TrialDay_2analyze{v})+1:end);
                        fullstreams_aligned_pooled_v2.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}) = data2incl;                   
                    end
                end
            end
        end
        fullstreams_aligned_pooled = fullstreams_aligned_pooled_v2;
        
        % adjust time vector
        for v=1:length(TrialDay_2analyze)
            length_data_aligned{v} = size(find(~isnan(fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(1,:))),2);
            time_vect_aligned_special.(TrialDay_2analyze{v}) = 0:dt_ds*1:dt_ds*(length_data_aligned{v}-1);
        end         
    end
    
    % Graphs
    % Create plot with nanmean and nanstd, Save graph 
    color2plot = {'m','r','b','k','c','y','g'};

    %% Plot 
    for k = 1; %:size(datatype,2) % dFF only k=1
        for d = 1; %:length(dFF_names)
            for p = 2; %1:length(pooledtype) % baselinecorr only p=2            
                % figure 1: overlap
                figure; clf; hold on
                if show_plot == 0
                    set(gcf,'visible','off')
                end

                for pow = 1:length(Opto_Powers)
                    for v=1:length(TrialDay_2analyze)
                        % average plot
                        tmp_avg_pool = nanmean(fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),1); 
                        tmp_error_pool = nanstd(fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),1,1)./...
                        sqrt(sum(~isnan(fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(:,1))));
                        error_area(time_vect_aligned_special.(TrialDay_2analyze{v}),tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector
% %                                                
                        % all individual trials; in case want to see individual trials
%                         for jkl=1:size(fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),1);
%                             data2plotnow = fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(jkl,:);
%                             plot(time_vect_aligned_special.(TrialDay_2analyze{v}),data2plotnow,color2plot{jkl}); hold on;
%                         end
                    end
                end

                % Formatting
                ylim([-5 15]);                    
                ax = gca;
                ax.FontSize = 20; 
                xlabel('Time (s)','FontSize', 20)
                ylabel('dFF','FontSize', 20)
                box off
                set(gca,'TickLength',[0 0])
                if length(TrialDay_2analyze)==1
                    L=legend(TrialDay_2analyze{1},' '); L.Location = 'Best';                
                elseif length(TrialDay_2analyze)== 2
                    L=legend(TrialDay_2analyze{1},' ',TrialDay_2analyze{2},' '); L.Location = 'Best';
                end
                %saveplot or not
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\streams_aligned_overlap',pooledtype{p},' ',datatype{k},'.tif']);
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\streams_aligned_overlap',pooledtype{p},' ',datatype{k},'.fig']);
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\streams_aligned_overlap',pooledtype{p},' ',datatype{k},'.pdf']);
                end          
            end
        end
    end

    %% Downsample data to transfer the graph into GraphPad 
    for v=1:length(TrialDay_2analyze)
        for d=1:length(dFF_names)
            for pow = 1:length(Opto_Powers)
                for k = 1; %:size(datatype,2) % here only dFF cos its longitudinal so can't do dFF within
                    for nummice=1:size(fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),1)
                        data2plot=fullstreams_aligned_pooled.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(nummice,:);
                        
                         
                        ttt = timetable(data2plot','SampleRate',sampling_rate_ds);
                        tttAvg = retime(ttt,'regular','mean','SampleRate',1); %1fps
    %                     figure; plot(tttAvg.Time,tttAvg.Var1)
                        fullstreams_aligned_pooled_ds.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(nummice,:) = tttAvg.Var1;
                    end
    %                 % option 1
    %                 time_vect_aligned_ds.(Opto_Powers{pow}) = [seconds(tttAvg.Time)] + time_vect_aligned.(Opto_Powers{pow})(1); %also add the first timepoint from the previous timevector allowing to shift in time
    %                 time_vect_aligned_ds_minutes.(Opto_Powers{pow}) = time_vect_aligned_ds.(Opto_Powers{pow})/60;
                    % option 2
                    time_vect_aligned_ds.(TrialDay_2analyze{v}) = [seconds(tttAvg.Time)] + time_vect_aligned_special.(TrialDay_2analyze{v})(1); %also add the first timepoint from the previous timevector allowing to shift in time                
                    time_vect_aligned_ds_minutes.(TrialDay_2analyze{v}) = time_vect_aligned_ds.(TrialDay_2analyze{v})/60;
                end
            end
        end
    end

    % Create plot with nanmean and nanstd, Save graph 
    color2plot = {'b','g','r','k','c','b','k','c'};

    for k = 1; %:size(datatype,2) % dFF only k=2
        for d=1:length(dFF_names)
            for p = 2; %1:length(pooledtype) % baselinecorr only p=2            
                % figure 1: overlap
                figure; clf; hold on
                if show_plot == 0
                    set(gcf,'visible','off')
                end

                for pow = 1:length(Opto_Powers)
                    for v=1:length(TrialDay_2analyze)
                        tmp_avg_pool = nanmean(fullstreams_aligned_pooled_ds.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),1); 
                        tmp_error_pool = nanstd(fullstreams_aligned_pooled_ds.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),1,1)./...
                        sqrt(sum(~isnan(fullstreams_aligned_pooled_ds.(TrialDay_2analyze{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(:,1))));
                        error_area(time_vect_aligned_ds.(TrialDay_2analyze{v}),tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector
                    end
                end                       

                % Formatting
                ylim([-5 15]);                    
                ax = gca;
                ax.FontSize = 20; 
                xlabel('Time (s)','FontSize', 20)
                ylabel('dFF','FontSize', 20)
                box off
                set(gca,'TickLength',[0 0])          
                if length(TrialDay_2analyze)==1
                    L=legend(TrialDay_2analyze{1},' '); L.Location = 'Best';                
                elseif length(TrialDay_2analyze)== 2
                    L=legend(TrialDay_2analyze{1},' ',TrialDay_2analyze{2},' '); L.Location = 'Best';
                end
                title('downsampled version for GraphPad','FontSize',14);

                %saveplot or not
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\streams_aligned_overlap_AVE1sec',pooledtype{p},' ',datatype{k},'.tif']);
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\streams_aligned_overlap_AVE1sec',pooledtype{p},' ',datatype{k},'.fig']);
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\streams_aligned_overlap_AVE1sec',pooledtype{p},' ',datatype{k},'.pdf']);
                end          
            end
        end
    end
end


%% Individual trials- average curves: ALL GROUPS TOGETHER ON SAME GRAPH
% todo_averagecurves_allgroups=1
if todo_averagecurves_allgroups == 1
    % Create plot with nanmean and nanstd, Save graph
    color2plot = {'m','r','b','k','c'};
    TrialType_Number = length(Opto_Powers);

    for k = datatype_2use_for_graphs; %:size(datatype,2) 
        for d=1:length(dFF_names)
            for p = 2; %1:length(pooledtype) % baselinecorr only p=2
                figure; clf; hold on
                if show_plot == 0
                    set(gcf,'visible','off')
                end

                for v=1:length(TrialDay_2analyze)
                    for pow = 1:length(Opto_Powers)
                        tmp_avg_pool = nanmean(PooledAVE.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1); 
                        tmp_error_pool = nanstd(PooledAVE.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1,1)./...
                        sqrt(sum(~isnan(PooledAVE.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(:,1))));
                        error_area(t_trials,tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector
                        hold on
                    end
                end

                % Formatting
                xlim([TRANGE(1) TRANGE(2)]);
                ax = gca;
                ax.FontSize = 20; 
                xlabel('Time (s)','FontSize', 20)
                ylabel('dFF','FontSize', 20)
                box off
                ylimits = get(gca,'YLim');
                xlimits = get(gca,'XLim');
                set(gca,'TickLength',[0 0])
                error_area_onlyrectangle([0 2],ylimits,[10*ylimits(2)-10*ylimits(1) 10*ylimits(2)-10*ylimits(1)],[0.0500, 0.0250, 0.0980],0.1); %(X,Y,barlength,color,alpha for shading of SEM). Barlength should be Ylimit_max - Ylimitmin
                plot([0 0],ylimits,'k','LineWidth',0.5,'LineStyle','-.')
                yline(0,'-.k');
                ylim([-5 10])
                % Legend
                for v=1:length(TrialDay_2analyze)
                    text4legend = TrialDay_2analyze{v};
                    text([xlimits(1)+1.5 xlimits(1)+1.5],[(ylimits(2)-0.1)-v*2 (ylimits(2)-0.1)-v*2],text4legend,'Color',color2plot{v},'FontSize',14)
                end

                %saveplot or not
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\AVE_alldays_alltrials_plots',pooledtype{p},' ',datatype{k},'.tif']);
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\AVE_alldays_alltrials_plots',pooledtype{p},' ',datatype{k},'.fig']);
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\AVE_alldays_alltrials_plots',pooledtype{p},' ',datatype{k},'.pdf']);

                end
            end
        end
    end
end


%% Individual trials- average curves: ALL GROUPS TOGETHER ON SAME GRAPH - PooledINDIV_VEHnormalized
% todo_averagecurves_allgroups=1
if todo_averagecurves_allgroups == 1
    % Create plot with nanmean and nanstd, Save graph
    color2plot = {'m','r','b','k','c'};
    TrialType_Number = length(Opto_Powers);

    for k = datatype_2use_for_graphs; %:size(datatype,2) 
        for d=1:length(dFF_names)
            for p = 2; %1:length(pooledtype) % baselinecorr only p=2
                figure; clf; hold on
                if show_plot == 0
                    set(gcf,'visible','off')
                end

                for v=1:length(TrialDay_2analyze)                
                    for pow = 1:length(Opto_Powers)
                        tmp_avg_pool = nanmean(PooledAVE_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1); 
                        tmp_error_pool = nanstd(PooledAVE_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1,1)./...
                        sqrt(sum(~isnan(PooledAVE_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(:,1))));
                        error_area(t_trials,tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector
                        hold on
                    end
                end

                % Formatting
                xlim([TRANGE(1) TRANGE(2)]);
                ax = gca;
                ax.FontSize = 20; 
                xlabel('Time (s)','FontSize', 20)
                ylabel('Normalized dFF responses','FontSize', 20)
                box off
                ylimits = get(gca,'YLim');
                xlimits = get(gca,'XLim');
                set(gca,'TickLength',[0 0])
                error_area_onlyrectangle([0 2],ylimits,[10*ylimits(2)-10*ylimits(1) 10*ylimits(2)-10*ylimits(1)],[0.0500, 0.0250, 0.0980],0.1); %(X,Y,barlength,color,alpha for shading of SEM). Barlength should be Ylimit_max - Ylimitmin
                plot([0 0],ylimits,'k','LineWidth',0.5,'LineStyle','-.')
                yline(0,'-.k');
                ylim([-20 200])
                % Legend
                for v=1:length(TrialDay_2analyze)
                    text4legend = TrialDay_2analyze{v};
                    text([xlimits(1)+1.5 xlimits(1)+1.5],[(ylimits(2)-0.1)-v*10 (ylimits(2)-0.1)-v*10],text4legend,'Color',color2plot{v},'FontSize',14)
                end

                %saveplot or not
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\AVE_alldays_alltrials_plots VEHnormalized',pooledtype{p},' ',datatype{k},'.tif']);
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\AVE_alldays_alltrials_plots VEHnormalized',pooledtype{p},' ',datatype{k},'.fig']);
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\AVE_alldays_alltrials_plots VEHnormalized',pooledtype{p},' ',datatype{k},'.pdf']);

                end
            end
        end
    end
end

%% %% Individual trials- curves: ALL GROUPS TOGETHER ON SAME GRAPH - ONLY MOUSE OF CHOICE
% todo_indivcurves_onemouse = 1
if todo_indivcurves_onemouse == 1
    %Edit nummice to mouse of choice below, eg nummiceX = 1 for 527_2 and trial 2 select if you're plotting individual trials
    nummiceX = 4;
    trial2select = 5;
    graphtype_indiv = 2; %1 if average+SEM for mouse, 2 if only 1 trial
    
    % Create plot with nanmean and nanstd, Save graph
    color2plot = {'m','r','b','k','c'};
    TrialType_Number = length(Opto_Powers);

    for k = datatype_2use_for_graphs; %:size(datatype,2) 
        for d=1:length(dFF_names)
            for p = 2; %1:length(pooledtype) % baselinecorr only p=2
                figure; clf; hold on
                if show_plot == 0
                    set(gcf,'visible','off')
                end
                
                % if you want the mouse's average + SEM
                if graphtype_indiv == 1
                    for v=1:length(TrialDay_2analyze)                 
                        for pow = 1:length(Opto_Powers)
                            tmp_avg_pool = PooledAVE.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummiceX,:); 
                            tmp_error_pool = PooledSEM.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummiceX,:);
                            error_area(t_trials,tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector
                        end
                    end
                elseif graphtype_indiv == 2
                % if you want the individual trial from a mouse
                    for v=1:length(TrialDay_2analyze)
                        for pow = 1:length(Opto_Powers)
                            data4plot = PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummiceX}).(SessionIDs{s})(trial2select,:); 
                            plot(t_trials,data4plot,color2plot{v}); %t_trials is the time vector
                        end
                    end
                end
                % Formatting
                xlim([TRANGE(1) TRANGE(2)]);
                ylim([-5 5]);
                ax = gca;
                ax.FontSize = 20; 
                xlabel('Time (s)','FontSize', 20)
                ylabel('dFF','FontSize', 20)
                title([AnimalIDs{nummiceX}],'FontSize', 12);
                box off
                ylimits = get(gca,'YLim');
                xlimits = get(gca,'XLim');
                set(gca,'TickLength',[0 0])
                %plot the opto area
                error_area_onlyrectangle([0 2],ylimits,[ylimits(2)-ylimits(1) ylimits(2)-ylimits(1)],[0.0500, 0.0250, 0.0980],0.1); %(X,Y,barlength,color,alpha for shading of SEM). Barlength should be Ylimit_max - Ylimitmin
                plot([0 0],ylimits,'k','LineWidth',0.5,'LineStyle','-.')
                yline(0,'-.k');
                % Legend, self-made - with text only 
                for v=1:length(TrialDay_2analyze)
                    text4legend = TrialDay_2analyze{v};
                    text([xlimits(1)+1.5 xlimits(1)+1.5],[(ylimits(2)-0.1)-v*0.8 (ylimits(2)-0.1)-v*0.8],text4legend,'Color',color2plot{v},'FontSize',14)
                end
                %saveplot or not
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\indivcurves_onemouse',AnimalIDs{nummiceX},' ',pooledtype{p},' ',datatype{k},'.tif']);
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\indivcurves_onemouse',AnimalIDs{nummiceX},' ',pooledtype{p},' ',datatype{k},'.fig']);
                    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\indivcurves_onemouse',AnimalIDs{nummiceX},' ',pooledtype{p},' ',datatype{k},'.pdf']);
                end
            end
        end
    end
end

%% Individual trials- heatmaps; one line per trial; all mice ; subplots VEH vs DETQ
% todo_heatmap_alltrials_allmice = 1;
if todo_heatmap_alltrials_allmice == 1
    figure;
    if show_plot == 0
        set(gcf,'visible','off')
    end
    hold on
    orderedmap = 0; % 0 to plot in normal order, 1 if you want to plot in the order of deg stim by sorting e.g. on the minimal or max value within defined period
    for v=1:length(TrialDay_2analyze)      
        for p = 2; %1:length(pooledtype)
            for d=1:length(dFF_names)
                for k =  datatype_2use_for_graphs; %1:size(datatype,2)
                    subplot(length(TrialDay_2analyze),1,v);
                    for pow = 1 %:length(Opto_Powers)
                        alltrials_array = [];
                        for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
                            animals = PooledAnimalID_only.(TrialDay_2analyze{v});
                            for s = 1:length(SessionIDs)
                                indivdata = PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s});
                                numrows = size(indivdata,1);
                                if isempty(alltrials_array)                               
                                    alltrials_array(1:numrows,:) = indivdata;
                                else
                                   numcurrentrows = size(alltrials_array,1);
                                   alltrials_array(numcurrentrows+1:numcurrentrows+numrows,:) = indivdata;
                                end
                            end
                        end
                        if orderedmap == 0
                            numrows = size(alltrials_array,1);
                            order = [1:numrows]';
                        % if you want each animal to have its trials sorted by order in terms of deg increase
                        elseif orderedmap == 1
                            [~,order] = sort(abs(max(alltrials_array(:,t_trials > 0 & t_trials < stim_duration),[],2))); %if A is a matrix, then min(A,[],2) is a column vector containing the minimum value of each row
                            %B = sort(A) sorts the elements of A in ascending order + [B,I] = sort(___) also returns a collection of index vectors
                        % if you want to put animals in order in order of how good they are (but their trials stay in order)
                        elseif orderedmap == 2
                           A= max(alltrials_array(:,t_trials > 0 & t_trials < stim_duration+3),[],2); 
                           numrows = size(alltrials_array,1);
                           numtrials = size(indivdata,1);
                           for i=1:length(mice_list_virus.(TrialDay_2analyze{v}))
                               if i==1
                                   mouse_score(i) = mean(A(1:numtrials))
                               else 
                                   mouse_score(i) = mean(A((numtrials*(i-1))+1:(numtrials*(i-1))+numtrials))
                               end
                           end
                           [~,order_overall] = sort(mouse_score);
                           order = ones(numrows,1)*nan;
                           for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
                               if nummice==1
                                    order(1:numtrials) = ((order_overall(nummice)-1)*numtrials)+1:((order_overall(nummice)-1)*numtrials)+numtrials;
                               else
                                    order((numtrials*(nummice-1))+1:(numtrials*(nummice-1))+numtrials) = ((order_overall(nummice)-1)*numtrials)+1:((order_overall(nummice)-1)*numtrials)+numtrials;
                               end
                           end
                        end
                        %plot
                        imagesc(t_trials,1,alltrials_array(order,:))
                        Color_scale = [-4 10];
                        colormap(jet); % colormap(flipud(autumn)) colormap(parula)
    %                     colorbar('Ticks',[-5,-2,1,4,7],'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'})
                        colorbar('Ticks',[0,3,6])
                        c = colorbar;
                        c.Label.String = '\DeltaF/F0';
                        if ~isempty(Color_scale)
                            c_limits = Color_scale;
                            caxis(c_limits)
                        end
                        xlim([-5 60])
                        numtrials_final.(TrialDay_2analyze{v}) = Trials_2analyze.(TrialDay_2analyze{v})(end) - Trials_2analyze.(TrialDay_2analyze{v})(1) +1;
                        % option 1
                        if length(TrialDay_2analyze) == 2 
                            if v==1
                            set(gca,'YTick',numtrials_final.(TrialDay_2analyze{v}):numtrials_final.(TrialDay_2analyze{v}):size(alltrials_array,1));
                            yticklabels({'1','2','3','4','5'})
                            ylabel('Mouse')
                        else
                            set(gca,'YTick',numtrials_final.(TrialDay_2analyze{v}):numtrials_final.(TrialDay_2analyze{v}):size(alltrials_array,1));
                            yticklabels({'1','2','3','4','5'})
                            ylabel('Mouse')                        
                        end
                        % ylines
                        if v==1
                           for i=1:length(mice_list_virus.(TrialDay_2analyze{v}))-1
                                yline(i*numtrials_final.(TrialDay_2analyze{v}),':k','Linewidth',1);
                           end
                        else
                            for i=1:length(mice_list_virus.(TrialDay_2analyze{v}))-1
                                yline(i*numtrials_final.(TrialDay_2analyze{v}),':k','Linewidth',1);
                            end
                        end
                        else % only DETQ
                            set(gca,'YTick',numtrials_final.(TrialDay_2analyze{v}):numtrials_final.(TrialDay_2analyze{v}):size(alltrials_array,1));
                            yticklabels({'1','2','3','4','5'})
                            ylabel('Mouse')
                            for i=1:length(mice_list_virus.(TrialDay_2analyze{v}))-1
                                yline(i*numtrials_final.(TrialDay_2analyze{v}),':k','Linewidth',1);
                            end
                        end
                        ax = gca;
                        ax.FontSize = 14; 
                        box off
                        ylim([1 size(alltrials_array,1)])
                        ylimits = get(gca,'YLim');
                        xlimits = get(gca,'XLim');
                        set(gca,'TickLength',[0 0])
                        xline(0,':k','Linewidth',2);
                        xline(stim_duration,':k','Linewidth',2);
                        %title
                        text4legend = TrialDay_2analyze{v};
                        t = title(text4legend);
                        t.FontSize = 14;
                        box off
                        set(findall(gcf,'-property','FontSize'),'FontSize',12)

                    end
                    if save_plot == 1
                        saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Heatmaps of dFF all mice all trials',pooledtype{p},' ',datatype{k},'.tif']);
                        saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Heatmaps of dFF all mice all trials',pooledtype{p},' ',datatype{k},'.fig']);
                        saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Heatmaps of dFF all mice all trials',pooledtype{p},' ',datatype{k},'.pdf']);
                    end
                end
            end
        end
    end
    xlabel('Time (s)','FontSize', 12)
end                
set(findall(gcf,'-property','FontSize'),'FontSize',12)
   

%% Individual trials- heatmaps; one line per trial; all mice ; subplots VEH vs DETQ - VEH normalized
% todo_heatmap_alltrials_allmice_VEHnormalized = 1;
if todo_heatmap_alltrials_allmice_VEHnormalized == 1
    f= figure; 
%     f.Position=[520 155 560 570]; % Use this for the bigger heatmaps with many trials
    if show_plot == 0
        set(gcf,'visible','off')
    end
    hold on
    orderedmap = 3; % 0 to plot in normal order, 1 if you want to plot in the order of deg stim by sorting e.g. on the minimal or max value on defined period, 
                        % 2 ordered but keep trials in order, 3 if manual set (See below)
    for v=1:length(TrialDay_2analyze) % say v=2 if only want DETQ in the heatmap
        for p = 2; %1:length(pooledtype)
            for d=1:length(dFF_names)
                for k =  datatype_2use_for_graphs; %1:size(datatype,2)
                    subplot(length(TrialDay_2analyze),1,v);
                    for pow = 1 %:length(Opto_Powers)
                        alltrials_array = [];
                        for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
                            animals = PooledAnimalID_only.(TrialDay_2analyze{v});
                            for s = 1:length(SessionIDs)
                                indivdata = PooledINDIV_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s});
                                numrows = size(indivdata,1);
                                if isempty(alltrials_array)                               
                                    alltrials_array(1:numrows,:) = indivdata;
                                else
                                   numcurrentrows = size(alltrials_array,1);
                                   alltrials_array(numcurrentrows+1:numcurrentrows+numrows,:) = indivdata;
                                end
                            end
                        end
                        if orderedmap == 0
                            numrows = size(alltrials_array,1);
                            order = [1:numrows]';
                        % if you want each animal to have its trials sorted by order in terms of deg increase
                        elseif orderedmap == 1
                            [~,order] = sort(abs(max(alltrials_array(:,t_trials > 0 & t_trials < stim_duration + 5),[],2))); %if A is a matrix, then min(A,[],2) is a column vector containing the minimum value of each row
                            %B = sort(A) sorts the elements of A in ascending order + [B,I] = sort(___) also returns a collection of index vectors
                        % if you want to put animals in order in order of how good they are (but their trials stay in order)
                        elseif orderedmap == 2
                           A= max(alltrials_array(:,t_trials > 0 & t_trials < stim_duration + 1),[],2); 
                           numrows = size(alltrials_array,1);
                           numtrials = size(indivdata,1);
                           for i=1:length(mice_list_virus.(TrialDay_2analyze{v}))
                               if i==1
                                   mouse_score(i) = mean(A(1:numtrials))
                               else 
                                   mouse_score(i) = mean(A((numtrials*(i-1))+1:(numtrials*(i-1))+numtrials))
                               end
                           end
                           [~,order_overall] = sort(mouse_score);
                           order = ones(numrows,1)*nan;
                           for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
                               if nummice==1
                                    order(1:numtrials) = ((order_overall(nummice)-1)*numtrials)+1:((order_overall(nummice)-1)*numtrials)+numtrials;
                               else
                                    order((numtrials*(nummice-1))+1:(numtrials*(nummice-1))+numtrials) = ((order_overall(nummice)-1)*numtrials)+1:((order_overall(nummice)-1)*numtrials)+numtrials;
                               end
                           end
                        elseif orderedmap == 3
                           numrows = size(alltrials_array,1);
                           numtrials = size(indivdata,1);
                           order_overall = [3 4 1 2 5]; 
                           order = ones(numrows,1)*nan;
                           for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
                               if nummice==1
                                    order(1:numtrials) = ((order_overall(nummice)-1)*numtrials)+1:((order_overall(nummice)-1)*numtrials)+numtrials;
                               else
                                    order((numtrials*(nummice-1))+1:(numtrials*(nummice-1))+numtrials) = ((order_overall(nummice)-1)*numtrials)+1:((order_overall(nummice)-1)*numtrials)+numtrials;
                               end
                           end                            
                        end
                        
                        %plot
                        imagesc(t_trials,1,alltrials_array(order,:))
                        Color_scale = [0 150];
                        xlim([-30 60])
                        colormap(parula)
                        colorbar('Ticks',[0,3,6])
                        c = colorbar;
                        c.Label.String = 'Normalized \DeltaF/F0 (%)';
                        if ~isempty(Color_scale)
                            c_limits = Color_scale;
                            caxis(c_limits)
                        end
                        numtrials_final.(TrialDay_2analyze{v}) = Trials_2analyze.(TrialDay_2analyze{v})(end) - Trials_2analyze.(TrialDay_2analyze{v})(1) +1;
                        % option 1
                        if length(TrialDay_2analyze) == 2 
                            if v==1
                            set(gca,'YTick',numtrials_final.(TrialDay_2analyze{v}):numtrials_final.(TrialDay_2analyze{v}):size(alltrials_array,1));
                            yticklabels({'1','2','3','4','5'})
                            ylabel('Mouse')
                        else
                            set(gca,'YTick',numtrials_final.(TrialDay_2analyze{v}):numtrials_final.(TrialDay_2analyze{v}):size(alltrials_array,1));
                            yticklabels({'1','2','3','4','5'})
                            ylabel('Mouse')                        
                        end
                        % ylines
                        if v==1
                           for i=1:length(mice_list_virus.(TrialDay_2analyze{v}))-1
                                yline(i*numtrials_final.(TrialDay_2analyze{v}),':k','Linewidth',1);
                           end
                        else
                            for i=1:length(mice_list_virus.(TrialDay_2analyze{v}))-1
                                yline(i*numtrials_final.(TrialDay_2analyze{v}),':k','Linewidth',1);
                            end
                        end
                        else % only DETQ
                            set(gca,'YTick',numtrials_final.(TrialDay_2analyze{v}):numtrials_final.(TrialDay_2analyze{v}):size(alltrials_array,1));
                            yticklabels({'1','2','3','4','5'})
                            ylabel('Mouse')
                            for i=1:length(mice_list_virus.(TrialDay_2analyze{v}))-1
                                yline(i*numtrials_final.(TrialDay_2analyze{v}),':k','Linewidth',1);
                            end
                        end
                        ax = gca;
                        ax.FontSize = 14; 
                        box off
                        ylim([1 size(alltrials_array,1)])
                        ylimits = get(gca,'YLim');
                        xlimits = get(gca,'XLim');
                        set(gca,'TickLength',[0 0])
                        xline(0,':k','Linewidth',2);
                        xline(stim_duration,':k','Linewidth',0.5);
                        % keeping track of time
                        time_trials = (Trials_2analyze.(TrialDay_2analyze{v})-1)*5;
                        time_trials_text = {};
                        for i=1:length(time_trials)
                            time_trials_text{i} = num2str(time_trials(i));
                        end
                        %title
                        %option 1
                        if v== 1; text4legend = ['Veh',' (Trials ',time_trials_text{1},'-',time_trials_text{end},' min post-drug)']; 
                        elseif v==2; text4legend = ['DETQ',' (Trials ',time_trials_text{1},'-',time_trials_text{end},' min post-drug)']; end;
                        % option 2 in case want to plot sthg different (for final fig for paper)
%                         if v== 1; text4legend = ['Veh']; elseif v==2; text4legend = ['Individual responses, DETQ']; end;
                        t = title(text4legend);
                        t.FontSize = 12;
                        t.FontWeight = 'normal';
                        ax.FontSize = 12; 
                        c.FontSize = 12;
                        box off
                    end
                end
            end
        end
    end
    xlabel('Time (s)','FontSize', 12)
end                   
if save_plot == 1
    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Heatmaps of dFF all mice all trials VEH normalized',pooledtype{p},' ',datatype{k},'.tif']);
    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Heatmaps of dFF all mice all trials VEH normalized',pooledtype{p},' ',datatype{k},'.fig']);
    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Heatmaps of dFF all mice all trials VEH normalized',pooledtype{p},' ',datatype{k},'.pdf']);
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% Average trials- heatmaps; one line per trial; average of all mice ; subplots VEH vs DETQ- VEH normalized
% todo_heatmap_avetrials_allmice = 1;
if todo_heatmap_avetrials_allmice == 1
    f= figure; 
    f.Position=[520 155 560 570]; % Use this for the bigger heatmaps with many trials
    if show_plot == 0
        set(gcf,'visible','off')
    end
    hold on
    orderedmap = 0; % 0 to plot in normal order, 1 if you want to plot in the order of deg stim by sorting e.g. on the minimal or max value in determined time period
    for v=1:length(TrialDay_2analyze) % say v=2 if only want DETQ in the heatmap
        for p = 2; %1:length(pooledtype)
            for d=1; %:length(dFF_names)
                for k = datatype_2use_for_graphs; %1:size(datatype,2)
                    if length(TrialDay_2analyze)>1
                        subplot(length(TrialDay_2analyze),1,v);
                    end
                    for pow = 1 %:length(Opto_Powers)
                         alltrials_array = [];
                        for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
                            animals = PooledAnimalID_only.(TrialDay_2analyze{v});
                            for s = 1:length(SessionIDs)
                                VEHnormalized = 1;
                                if VEHnormalized == 1;
                                    indivdata = PooledINDIV_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s});                                
                                else    
                                    indivdata = PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s});
                                end
                                numrows = size(indivdata,1);
                                if isempty(alltrials_array)                               
                                    alltrials_array(1:numrows,:) = indivdata;
                                else
                                   numcurrentrows = size(alltrials_array,1);
                                   alltrials_array(numcurrentrows+1:numcurrentrows+numrows,:) = indivdata;
                                end
                            end
                        end
                        % get the mean
                        meantrials_array = ones(size(indivdata,1),size(indivdata,2))*nan;
                        number_animals = length(mice_list_virus.(TrialDay_2analyze{v}));
                        for ki = 1:size(indivdata,1)
                            meantrials_array(ki,:) = nanmean(alltrials_array(ki:numrows:size(alltrials_array,1),:),1);
                        end
                        % if you want mice sorted by order in terms of deg inhibition
                        if orderedmap == 1
                            [~,order] = sort(abs(min(meantrials_array(:,t_trials > 2 & t_trials < stim_duration+3),[],2))); %if A is a matrix, then min(A,[],2) is a column vector containing the minimum value of each row
                            %B = sort(A) sorts the elements of A in ascending order + [B,I] = sort(___) also returns a collection of index vectors
                        else
                            order = [1:numrows]';
                        end
                        imagesc(t_trials,1,meantrials_array(order,:))
                        %plot
                        if VEHnormalized == 1;
                            Color_scale = [0 150];
                        else    
                            Color_scale = [0 5];
                        end
                        colormap(parula)
                        colorbar('Ticks',[0,3,6])
                        c = colorbar;
                        if VEHnormalized == 1;
                            c.Label.String = 'Normalized \DeltaF/F0 (%)';
                        else    
                            c.Label.String = '\DeltaF/F0';
                        end
                        if ~isempty(Color_scale)
                            c_limits = Color_scale;
                            caxis(c_limits)
                        end
                        xlim([-30 60])
                        % keeping track of time
                        time_trials = (Trials_2analyze.(TrialDay_2analyze{v})-1)*5;
                        time_trials_text = {};
                        for i=1:length(time_trials)
                            time_trials_text{i} = num2str(time_trials(i));
                        end
                        % option 1
                        if v==1
                            set(gca,'YTick',1:1:size(alltrials_array,1));                            
                            yticklabels(time_trials_text)
                            ylabel('min post-drug')
                        else
                            set(gca,'YTick',2:2:size(alltrials_array,1)); % can edit this to 1:1 and in the line below also; if want to plot all timepoints
                            yticklabels(time_trials_text(2:2:end))
                            ylabel('min post-drug')                        
                        end
                        box off
                        ylim([1 size(alltrials_array,1)]);
                        ax = gca;
                        box off
                        ylim([0.5 size(meantrials_array,1)+0.5])
                        ylimits = get(gca,'YLim');
                        xlimits = get(gca,'XLim');
                        set(gca,'TickLength',[0 0])
                        xline(0,':k','Linewidth',2);
                        xline(stim_duration,':k','Linewidth',0.5);
                        %title
                        %option 1
                        if v== 1; text4legend = ['Veh']; elseif v==2; text4legend = ['DETQ']; end;
                        % option 2 in case want to plot sthg different (for final fig for paper)
%                         if v== 1; text4legend = ['Veh']; elseif v==2; text4legend = ['Average response, DETQ']; end;
                        t = title(text4legend);
                        t.FontSize = 12;
                        t.FontWeight = 'normal';
                        ax.FontSize = 12; 
                        c.FontSize = 12;
                    end
                end
            end
        end
    end
    xlabel('Time (s)','FontSize', 12)
end
if save_plot == 1
    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Heatmaps of dFF all mice AVE trials',pooledtype{p},' ',datatype{k},'.tif']);
    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Heatmaps of dFF all mice AVE trials',pooledtype{p},' ',datatype{k},'.fig']);
    saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Heatmaps of dFF all mice AVE trials',pooledtype{p},' ',datatype{k},'.pdf']);
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)


%% Individual trials- heatmaps; one line per trial; ONLY MOUSE OF CHOICE ; VEH VS DETQ
nummiceX = 1;
% todo_heatmap_alltrials_onemouse = 1;
if todo_heatmap_alltrials_onemouse == 1
    figure;
    if show_plot == 0
        set(gcf,'visible','off')
    end
    hold on
    orderedmap = 2; % 0 to plot in normal order, 1 if you want to plot in the order of deg stim by sorting e.g. on the min or max value
    for v=1:length(TrialDay_2analyze)   
        for p = 2; %1:length(pooledtype)
            for d=1:length(dFF_names)
                for k = datatype_2use_for_graphs; %1:size(datatype,2)
                    if length(TrialDay_2analyze)>1
                        subplot(length(TrialDay_2analyze),1,v);
                    end
                    for pow = 1 %:length(Opto_Powers)
                        alltrials_array = [];
                        for nummice=nummiceX
                            animals = PooledAnimalID_only.(TrialDay_2analyze{v});
                            for s = 1:length(SessionIDs)
                                indivdata = PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s});
                                numrows = size(indivdata,1);
                                if isempty(alltrials_array)                               
                                    alltrials_array(1:numrows,:) = indivdata;
                                else
                                   numcurrentrows = size(alltrials_array,1);
                                   alltrials_array(numcurrentrows+1:numcurrentrows+numrows,:) = indivdata;
                                end
                            end
                        end
                        %plot
                        imagesc(t_trials,1,alltrials_array)
                        Color_scale = [-4 10];
                        colormap(flipud(autumn))
                        colorbar('Ticks',[0,3,6])
                        c = colorbar;
                        c.Label.String = '\DeltaF/F0';
                        if ~isempty(Color_scale)
                            c_limits = Color_scale;
                            caxis(c_limits)
                        end
                        xlim([-5 60])
                        % option 1
                        numtrials_final.(TrialDay_2analyze{v}) = Trials_2analyze.(TrialDay_2analyze{v})(end) - Trials_2analyze.(TrialDay_2analyze{v})(1) +1;
                        if v==1
                            set(gca,'YTick',numtrials_final.(TrialDay_2analyze{v}):numtrials_final.(TrialDay_2analyze{v}):size(alltrials_array,1));
                            ylabel('Trials')
                        else
                            set(gca,'YTick',numtrials_final.(TrialDay_2analyze{v}):numtrials_final.(TrialDay_2analyze{v}):size(alltrials_array,1));
                            ylabel('Trials')                        
                        end
                        % ylines
                        if v==1
                           for i=1:length(mice_list_virus.(TrialDay_2analyze{v}))-1
%                                 yline(i*6,':k','Linewidth',1);
                           end
                        else
                            for i=1:length(mice_list_virus.(TrialDay_2analyze{v}))-1
%                                 yline(i*18,':k','Linewidth',1);
                            end
                        end
                        ylim([1 size(alltrials_array,1)])
                        ax = gca;
                        ax.FontSize = 14; 
                        box off
                        ylimits = get(gca,'YLim');
                        xlimits = get(gca,'XLim');
                        set(gca,'TickLength',[0 0])
                        xline(0,':k','Linewidth',2);
                        xline(stim_duration,':k','Linewidth',2);
                        %title
                        text4legend = TrialDay_2analyze{v};
                        t = title(text4legend);
                        t.FontSize = 14;
                        box off
                        ylimits = get(gca,'YLim');
                    end
                    if save_plot == 1
                        saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Heatmaps dFF 1 mouse all trials',AnimalIDs{nummiceX},' ',pooledtype{p},' ',datatype{k},'.tif']);
                        saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Heatmaps dFF 1 mouse all trials',AnimalIDs{nummiceX},' ',pooledtype{p},' ',datatype{k},'.pdf']);
                        saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Heatmaps dFF 1 mouse all trials',AnimalIDs{nummiceX},' ',pooledtype{p},' ',datatype{k},'.fig']);
                    end
                end
            end
        end
    end
    xlabel('Time (s)','FontSize', 14)
end 

%% Plot the results per trial types separately over time average - VERTICAL
% todo_vertical_avetrials_allmice = 1;
if todo_vertical_avetrials_allmice == 1
    startgraph = -10;
    dottedlines = [1:1:1];
    scalingfactor = -85;
   
    for v=1:length(TrialDay_2analyze)
        for p = 2; %1:length(pooledtype)
            for d=1; %:length(dFF_names)
                for k = datatype_2use_for_graphs; %1:size(datatype,2)
%                     subplot(length(TrialDay),1,v);
                    for pow = 1 %:length(Opto_Powers)
                        alltrials_array = [];
                        for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
                            animals = PooledAnimalID_only.(TrialDay_2analyze{v});
                            for s = 1:length(SessionIDs)
                                VEHnormalized = 1;
                                if VEHnormalized == 1;
                                    indivdata = PooledINDIV_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s});                                
                                else 
                                    indivdata = PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s});
                                end
                                numrows = size(indivdata,1);
                                if isempty(alltrials_array)                               
                                    alltrials_array(1:numrows,:) = indivdata;
                                else
                                   numcurrentrows = size(alltrials_array,1);
                                   alltrials_array(numcurrentrows+1:numcurrentrows+numrows,:) = indivdata;
                                end
                            end
                        end
                        % get the mean
                        meantrials_array = ones(size(indivdata,1),size(indivdata,2))*nan;
                        number_animals = length(mice_list_virus.(TrialDay_2analyze{v}));
                        for ki = 1:size(indivdata,1)
                            meantrials_array(ki,:) = nanmean(alltrials_array(ki:number_animals:end,:),1);
                            semtrials_array(ki,:) = nanstd(alltrials_array(ki:number_animals:end,:),1)./sqrt(number_animals);
                        end
                        
                        % plot
                        data2plot = meantrials_array;
                        figure; set(gcf, 'Position',  [100+(400*v), 70, 400, 45*size(data2plot,1)])
                        if show_plot == 0
                            set(gcf,'visible','off')
                        end
                        % keeping track of time
                        time_trials = (Trials_2analyze.(TrialDay_2analyze{v})-1)*5;
                        time_trials_text = {};
                        for i=1:length(time_trials)
                            time_trials_text{i} = num2str(time_trials(i));
                        end
                        for q=1:size(data2plot,1)
                            plot(t_trials,data2plot(q,:)+q*scalingfactor,'LineWidth',2);   
                            yline(q*scalingfactor,'-.k');
                            text(startgraph-4.5,q*(scalingfactor),time_trials_text{q},'Color','black','FontSize',10);
                            text(startgraph-3,q*(scalingfactor),'min','Color','black','FontSize',10);
                            hold on
                        end
                        xline(0,'-k');
                        xline(stim_duration,'-m');
                        xlim([startgraph 20])
                        xlabel('Time (s)')
                        set(gca,'YTick',[]);
                        set(gca,'Yticklabel',[]);
                        set(gca,'YLabel',[]);
                        ax=gca; ax.YAxis.Color = 'w';
                        set(gca,'Box','off');
                        yline(0,'-.k');
                        ax = gca;
                        ax.FontSize = 14; 
                        for u=1:length(dottedlines)-1
                            yline(size(data2plot,1)*scalingfactor+(scalingfactor*dottedlines(u)),'-.k');
                        end
                        sgtitle([pooledtype{p},' ',datatype{k},' ',dFF_names{d},' ',TrialDay_2analyze{v}]); % for groups of subplots
                        if save_plot == 1
                            saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Vertical dFF over time all mice AVE trials',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.tif']);
                            saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Vertical dFF over time all mice AVE trials',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.fig']);
                            saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Vertical dFF over time all mice AVE trials',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.pdf']);
                        end
                    end
                end
            end
        end
    end
end
                
%% Plot the results per trial types over time 1 mouse- only some trials not all (specify below), overlay on same graph
% todo_vertical_alltrials_1mouse = 1;
if todo_vertical_alltrials_1mouse == 1
    startgraph = -10;
    dottedlines = [1:1:1];
    scalingfactor = 100;
    nummiceX = 1;
    for v=1:length(TrialDay_2analyze)  
        for p = 2; %1:length(pooledtype)
            for d=1; %:length(dFF_names)
                for k = datatype_2use_for_graphs; %1:size(datatype,2)
%                     subplot(length(TrialDay),1,v);
                    for pow = 1 %:length(Opto_Powers)
                        alltrials_array = [];
                        for nummice=nummiceX;
                            animals = PooledAnimalID_only.(TrialDay_2analyze{v});
                            for s = 1:length(SessionIDs)
                                indivdata = PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s});
                                numrows = size(indivdata,1);
                                if isempty(alltrials_array)                               
                                    alltrials_array(1:numrows,:) = indivdata;
                                else
                                   numcurrentrows = size(alltrials_array,1);
                                   alltrials_array(numcurrentrows+1:numcurrentrows+numrows,:) = indivdata;
                                end
                            end
                        end
                        % get the mean
                        meantrials_array = ones(size(indivdata,1),size(indivdata,2))*nan;
                        number_animals = length(mice_list_virus.(TrialDay_2analyze{v}));
                        for ki = 1:size(indivdata,1)
                            meantrials_array(ki,:) = nanmean(alltrials_array(ki:number_animals:end,:),1);
                        end
                        
                        % plot
                        data2plot = meantrials_array;
                        figure; set(gcf, 'Position',  [100+(400*v), 70, 400, 45*size(data2plot,1)])
                        if show_plot == 0
                            set(gcf,'visible','off')
                        end
                         % keeping track of time
                        time_trials = (Trials_2analyze.(TrialDay_2analyze{v})-1)*5;
                        time_trials_text = {};
                        for i=1:length(time_trials)
                            time_trials_text{i} = num2str(time_trials(i));
                        end
                        for q=1:size(data2plot,1)
                            plot(t_trials,data2plot(q,:)+q*scalingfactor,'LineWidth',2);   
                            yline(q*scalingfactor,'-.k');
%                             time2plot = time_optotrials_aligned_simple.(TrialDay_2analyze{v})(q,1);
                            text(startgraph-4.5,q*(scalingfactor),time_trials_text{q},'Color','black','FontSize',8);
                            text(startgraph-3,q*(scalingfactor),'min','Color','black','FontSize',8);
                            hold on
                        end
                        xline(0,'-k');
                        xline(stim_duration,'-m');
                        xlim([startgraph 20])
                        xlabel('Time (s)')
                        ylabel(datatype{k});
                        set(gca,'YTick',[]);
                        set(gca,'Yticklabel',[]);
                        set(gca,'YLabel',[]);
                        ax=gca; ax.YAxis.Color = 'w';
                        set(gca,'Box','off');
                        yline(0,'-.k');
                        ax = gca;
                        ax.FontSize = 14; 
                        for u=1:length(dottedlines)-1
                            yline(size(data2plot,1)*scalingfactor+(scalingfactor*dottedlines(u)),'-.k');
                        end
                        sgtitle([pooledtype{p},' ',datatype{k},' ',dFF_names{d},' ',animals{nummiceX},TrialDay_2analyze{v},' ']); % for groups of subplots
                        if save_plot == 1
                            saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Vertical dFF over time 1mouse',TrialDay_2analyze{v},' ',animals{nummiceX},pooledtype{p},' ',datatype{k},'.tif']);
                            saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Vertical dFF over time 1mouse',TrialDay_2analyze{v},' ',animals{nummiceX},pooledtype{p},' ',datatype{k},'.fig']);
                            saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Vertical dFF over time 1mouse',TrialDay_2analyze{v},' ',animals{nummiceX},pooledtype{p},' ',datatype{k},'.pdf']);
                        end
                    end
                end
            end
        end
    end
end


%% Plot the results per trial types separately over time - OVERLAY. This is VEH normalized

% todo_overlay_alltrials_tempwindow_allmice = 1;
    color2plot = {'b','g','r','k','c','m','y','b','g','r','k','c','m','y','b','g','r','k','c','m','y'};
if todo_overlay_alltrials_tempwindow_allmice == 1
    % empty
    for v=1:length(TrialDay_2analyze)   
        for p = 2; %1:length(pooledtype)
            for d=1; %:length(dFF_names)
                for k = datatype_2use_for_graphs; %1:size(datatype,2)
                    for pow = 1 %:length(Opto_Powers)
                        PooledAVETrials.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(size(Trials_2analyze.(TrialDay_2analyze{v}),2),length(t_trials))*nan;
                        PooledSEMTrials.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(size(Trials_2analyze.(TrialDay_2analyze{v}),2),length(t_trials))*nan;
                    end
                end
            end
        end
    end
    
    % fill                   
    for v=1:length(TrialDay_2analyze)   
        for p = 2; %1:length(pooledtype)
            for d=1; %:length(dFF_names)
                for k = datatype_2use_for_graphs; %1:size(datatype,2)
                    for pow = 1 %:length(Opto_Powers)
                        for kh = 1:length(Trials_2analyze.(TrialDay_2analyze{v}))
                            VEHnormalized = 1;
                            if VEHnormalized == 1                                
                                for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
                                    temp_datastorage(nummice,:) = PooledINDIV_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s})(kh,:);
                                end
                            else
                                for nummice=1:length(mice_list_virus.(TrialDay_2analyze{v}))
                                    temp_datastorage(nummice,:) = PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s})(kh,:);
                                end                                
                            end
                            PooledAVETrials.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(kh,:) = ...
                                nanmean(temp_datastorage,1);
                            PooledSEMTrials.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(kh,:) = ...
                                nanstd(temp_datastorage,1,1)./sqrt(sum(~isnan(temp_datastorage(:,1))))
                            clear temp_datastorage
                        end
                    end
                end
            end
        end
    end 
        
    % figure
    figure; 
    for v=1:length(TrialDay_2analyze)
        subplot(1,2,v)
        for p = 2; %1:length(pooledtype)
            for d=1; %:length(dFF_names)
                for k = datatype_2use_for_graphs; %1:size(datatype,2)
                    for pow = 1 %:length(Opto_Powers)                        
                        % only plots
%                         for kh = 1:length(Trials_2analyze.(TrialDay_2analyze{v}))
%                             A_plot = PooledAVETrials.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(kh,:); 
%                             plot(t_trials,A_plot,'Color',color2plot{kh}); hold on; %t_trials is the time vector                           
%                         end
                        % AVE and SEM
                        for kh = 1:length(Trials_2analyze.(TrialDay_2analyze{v}))                        
                            A_plot = PooledAVETrials.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(kh,:);          
                            SEM_plot = PooledSEMTrials.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(kh,:);
                            error_area(t_trials,A_plot,SEM_plot,color2plot{kh},0.25); %t_trials is the time vector
                        end                                                
                        xlim([-10 20]);
                        if VEHnormalized == 1                                
                            ylim([-50 250]);
                        else
                            ylim([-5 10]);
                        end                        
                        xline(0,'-k');
                        xline(stim_duration,':k');
                        xlabel('Time (s)')
                        ylabel('Normalized \DeltaF/F0')
                        ax=gca; ax.YAxis.Color = 'k';
                        set(gca,'Box','off');
                        yline(0,'-.k');
                        text4legend = {};
                        for er=1:size(time_optotrials_aligned_simple_cropped.(TrialDay_2analyze{v}),1)
                            text4legend{er} = [num2str(time_optotrials_aligned_simple_cropped.(TrialDay_2analyze{v})(er,1)),' min'];
                        end
                        if length(text4legend) == 1
                            L=legend(text4legend{1},text4legend{1}); L.Location = 'Best';                
                        elseif length(text4legend) == 2
                            L=legend(text4legend{1},text4legend{1},text4legend{2},text4legend{2}); L.Location = 'Best';                
                        elseif length(text4legend) == 3
                            L=legend(text4legend{1},text4legend{1},text4legend{2},text4legend{2},text4legend{3},text4legend{3}); L.Location = 'Best';                
                        elseif length(text4legend) == 4
                            L=legend(text4legend{1},text4legend{1},text4legend{2},text4legend{2},text4legend{3},text4legend{3},text4legend{4},text4legend{4}); L.Location = 'Best';                
                        elseif length(text4legend) == 5
                            L=legend(text4legend{1},text4legend{1},text4legend{2},text4legend{2},text4legend{3},text4legend{3},text4legend{4},text4legend{4},text4legend{5},text4legend{5}); L.Location = 'Best';                
                        elseif length(text4legend) == 6
                            L=legend(text4legend{1},text4legend{1},text4legend{2},text4legend{2},text4legend{3},text4legend{3},text4legend{4},text4legend{4},text4legend{5},text4legend{5},...
                                text4legend{6},text4legend{6}); L.Location = 'Best';                
                        elseif length(text4legend) == 7
                            L=legend(text4legend{1},text4legend{1},text4legend{2},text4legend{2},text4legend{3},text4legend{3},text4legend{4},text4legend{4},text4legend{5},text4legend{5},...
                                text4legend{6},text4legend{6},text4legend{7},text4legend{7}); L.Location = 'Best';                
                        elseif length(text4legend) == 8
                            L=legend(text4legend{1},text4legend{1},text4legend{2},text4legend{2},text4legend{3},text4legend{3},text4legend{4},text4legend{4},text4legend{5},text4legend{5},...
                                text4legend{6},text4legend{6},text4legend{7},text4legend{7},text4legend{8},text4legend{8}); L.Location = 'Best';                
                        end                            
                        ax = gca;
                        ax.FontSize = 14;
                        if VEHnormalized == 1                                
                            sgtitle([pooledtype{p},' ',datatype{k},' ',dFF_names{d},' ',TrialDay_2analyze{v},' VEH normalized']); % for groups of subplots
                        else
                            sgtitle([pooledtype{p},' ',datatype{k},' ',dFF_names{d},' ',TrialDay_2analyze{v},' not normalized']); % for groups of subplots                            
                        end
                        if save_plot == 1
                            saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Overlay dFF over time all mice',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.tif']);
                            saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Overlay dFF over time all mice',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.fig']);
                            saveas(gcf,[PATH2SAVEPOOL_SELECT,'\pooled figures\Overlay dFF over time all mice',TrialDay_2analyze{v},' ',pooledtype{p},' ',datatype{k},'.pdf']);
                        end
                    end
                end
            end
        end
    end
end


