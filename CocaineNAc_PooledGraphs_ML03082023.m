% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT
% Name: CocaineNAc_PooledGraphs
% Cohort: NAc AllodLightctrl with opto and drugs (DETQ or cocaine)

% Feb 2023 - Marie Labouesse, marie.labouesse@gmail.com

% 1- POOLED DATA
% loads FP data previously analyzed (with: CocaineNAc_DataExtraction) in matlab space "IndividualData.mat" from multiple animals (select folder containing multiple animals) 
% puts individual trials into one big matrix (all mice together): PooledINDIV
% calculates the mean of individual trial trajectories aligned to specific event for each mouse: PooledAVE
% saves a workspace that can then be used for quantifications (CocaineNAc_PooledQuantification)

% 2- POOLED GRAPHS
% generates average graphs curves (this used to be in the main script)
% generates heatmaps (1 average trace/mouse or all trials/all mice in one big heatmap - separated by power)

% edit relevant parameters in %% SETUP PARAMETERS

% functions needed: N/A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       
%% INITIALIZATIONS
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% PARAMETERS 
show_plot = 1; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
reanalysis = 0; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed. This is based on the existence of the matlab space "IndividualData.mat"
overwrite = 0; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not          
pooledtype = {'raw','baselinecorr'};

%% Define the path where the data is
PATH2DATA_0 = uigetdir('select folder'); %select the overall folder, where the VEH and DETQ
TrialDay = {'VEH','DETQ'}; 
PATH2SAVEPOOL = PATH2DATA_0;
mkdir([PATH2SAVEPOOL,'\pooled figures\']);
% Opto_Powers = {'SAL','COC'}; %  % imported via matlab .mat

for v=1:length(TrialDay)
    PATH2DATA.(TrialDay{v}) = [PATH2DATA_0,'\',TrialDay{v}];
    PATH2SAVEFOLDER.(TrialDay{v}) = [PATH2DATA_0,'\',TrialDay{v}];
    mice_list_virus.(TrialDay{v}) = dir(PATH2DATA.(TrialDay{v})); %all things in this folder
    mkdir([PATH2SAVEFOLDER.(TrialDay{v}),'\pooled figures\']);
    mkdir([PATH2SAVEFOLDER.(TrialDay{v}),'\pooled data\']);
end

                

%% IDENTIFY MICE TO ANALYZE 
for v=1:length(TrialDay)
    for o = length(mice_list_virus.(TrialDay{v})):-1:1
        if mice_list_virus.(TrialDay{v})(o).isdir == 0  %remove non-folders
            mice_list_virus.(TrialDay{v})(o) = [];
        else
            if  strcmp(mice_list_virus.(TrialDay{v})(o).name,'data') == 1 || strcmp(mice_list_virus.(TrialDay{v})(o).name,'figures') == 1 ...   %remove folders with data or figures or ., ..
                || contains(mice_list_virus.(TrialDay{v})(o).name,'data') || contains(mice_list_virus.(TrialDay{v})(o).name,'figures') ...
                || contains(mice_list_virus.(TrialDay{v})(o).name,'results')  || contains(mice_list_virus.(TrialDay{v})(o).name,'other') ......
                || strcmp(mice_list_virus.(TrialDay{v})(o).name,'.') == 1 || strcmp(mice_list_virus.(TrialDay{v})(o).name,'..') == 1
                mice_list_virus.(TrialDay{v})(o) = [];
            end
        end
    end
    Nmice_virus{v} = length(mice_list_virus.(TrialDay{v}));
end
              

%% Import individual workspaces and load the data into a pooled array
for v=1:length(TrialDay)
    %AnimalIDs
    for nummice=1:length(mice_list_virus.(TrialDay{v}))
        AnimalIDs(nummice) = {['ID',mice_list_virus.(TrialDay{v})(nummice).name(end-5:end)]};
    end   
    
    for nummice=1:length(mice_list_virus.(TrialDay{v}))
        % Define the path to the mouse and find the folders to analyze: 
        path2mouse = [PATH2DATA.(TrialDay{v}),'\',mice_list_virus.(TrialDay{v})(nummice).name,'\']; 
        % if there are several sessions inside the mouse folder
        sessions = dir(path2mouse);
        %remove non relevant folders
        for o = length(sessions):-1:1
            if sessions(o).isdir == 0  %remove non-folders
                sessions(o) = [];
            elseif strcmp(sessions(o).name,'data') == 1 || strcmp(sessions(o).name,'figures') == 1 || strcmp(sessions(o).name,'results') == 1 ...
                 || strcmp(sessions(o).name,'.') == 1 || strcmp(sessions(o).name,'..') == 1 || contains(sessions(o).name,'.csv') == 1
                 sessions(o) = [];
            end
        end
        if isempty(sessions)
            sessions(1).name = []; %to have an existent session
        end
        
        %SessionIDs
        for s=1:length(sessions)
            SessionIDs(s) = {['SessionNum',num2str(s)]};
            SessionNames(s) = {sessions(s).name};   
        end       
       
        % Create a folder to save the data for this mouse and define the path
        if exist([PATH2SAVEFOLDER.(TrialDay{v}),'\',mice_list_virus.(TrialDay{v})(nummice).name,'\'],'dir') == 0
            mkdir([PATH2SAVEFOLDER.(TrialDay{v}),'\',mice_list_virus.(TrialDay{v})(nummice).name,'\']);
        end
        path2save_mouse = [PATH2SAVEFOLDER.(TrialDay{v}),'\',mice_list_virus.(TrialDay{v})(nummice).name,'\'];
        
        
        %% Loop for all the sessions for the mouse: Sessions are experimental days you want to average together (replicates of each other) 
        for s = 1:length(sessions)
            % Define the path to the session and create folder to save if needed:
            PATH2SESSION = [path2mouse,sessions(s).name];
            if exist([path2save_mouse,sessions(s).name],'dir') == 0
                mkdir([path2save_mouse,sessions(s).name])
            end
            if length(sessions) == 1
                PATH2SAVE = path2save_mouse;
            else
                PATH2SAVE = [path2save_mouse,sessions(s).name,'\'];
            end
            if exist([PATH2SAVE,'figures'],'dir') == 0
                mkdir([PATH2SAVE,'figures'])
            end
            if exist([PATH2SAVE,'results'],'dir') == 0
                mkdir([PATH2SAVE,'results'])
            end
        
            % Check if results are already saved for this session 
            done = exist([PATH2SAVE,'PooledAllMice.mat'],'file'); % ........
            if done == 0 || reanalysis == 1 %if not done or if wish to reanalyze
        
                % load the mouse/session matlab space
                load([PATH2SESSION,'\IndividualData.mat']);
 
                % Initialization of pooled structure with individual trials   %.................... ONLY ONE SESSION PER MOUSE

                if nummice == 1 && s == 1
                    for d=1:length(dFF_names)
                        for pow = 1:length(Opto_Powers)
                            for k = 1:size(datatype,2)
                                % ANIMAL ID .................... add session ID later if needed
                                PooledAnimalID.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = {};
                                % INDIVIDUAL TRIALS OF EACH MICE FOR ALL MICE
                                PooledINDIV.(TrialDay{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = ones(size(IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})))*nan;  
                                PooledINDIV.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan;  
                                % AVERAGE TRIAL OF EACH MICE (ALL SESSIONS) FOR ALL MICE, EACH MOUSE ON ONE LINE
                                PooledAVE.(TrialDay{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan;  
                                PooledAVE.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan; 
                                % SEM
                                PooledSEM.(TrialDay{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan;  
                                PooledSEM.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan; 
                                % STREAMS                
                                fullstreams.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = {};
                                % STREAMS ALIGNED                
                                fullstreams_aligned.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = {};
                            end
                        end
                    end
                end
                
                % Add data one animal at a time
                % ANIMAL ID and SESSION ID
                PooledAnimalID.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = SessionNames(s);
                PooledAnimalID_only.(TrialDay{v}) = AnimalIDs;
                % STORE INDIVIDUAL TRIALS
                for d=1:length(dFF_names)
                    for pow = 1:length(Opto_Powers)
                        for k = 1:size(datatype,2)
                            PooledINDIV.(TrialDay{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}); 
                            PooledINDIV.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}); 
                        end
                    end
                end                            
                % STORE AVERAGE
                for d=1:length(dFF_names)
                    for pow = 1:length(Opto_Powers)
                        for k = 1:size(datatype,2)
                            if length(sessions) == 1
                                PooledAVE.(TrialDay{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = nanmean(IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1); 
                                PooledAVE.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = nanmean(IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1); 
                            
                                PooledSEM.(TrialDay{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = nanstd(IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1,1)./...
                                    sqrt(size(IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1));
                                
                                PooledSEM.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = nanstd(IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1,1)./...
                                    sqrt(size(IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1));
                                                       
                            elseif length(sessions) > 1
                                % NA
                            end
                        end
                    end
                end
                %STORE STREAMS
                fullstreams.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = streams;
                           
                %STORE STREAMS ALIGNED TO OPTO ONSET #1
                fullstreams_aligned.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = streams_aligned;
 
            end
        end
    end
end

%% Save pooled data for all viruses together
% STORE INDIVIDUAL and AVERAGE TRIALS
save([PATH2SAVEPOOL,'\PooledAllMice.mat'],'PooledINDIV','PooledAVE','PooledSEM','AnimalIDs','mice_list_virus','TrialDay','PooledAnimalID','PooledAnimalID_only','dFF_names','datatype','Opto_Powers',...
    't_trials','stim_duration','TRANGE','BASELINE_WIN','Events2Align2','dt_ds','sampling_rate_ds','length_data','fullstreams','fullstreams_aligned','time_vect','pooledtype','drug_name');



%% Streams aligned on a graph; ALL DRUGS TOGETHER ON SAME GRAPH, all animals
% CHOP OFF LAST 60-180 SECONDS
chop_off_ix = round(60./dt_ds); % 60 seconds

% create empty pooled array
length_fullstreams_aligned_pooled = length(fullstreams_aligned.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}));
for v=1:length(TrialDay)
    for d=1:length(dFF_names)
        for pow = 1:length(Opto_Powers)
            for k = 1:size(datatype,2)
                fullstreams_aligned_pooled.(TrialDay{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}) = ones(nummice,length_fullstreams_aligned_pooled-chop_off_ix+1)*nan;
                length_updated = length(fullstreams_aligned_pooled.(TrialDay{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}));
            end
        end
    end
end

%STORE STREAMS ALIGNED TO OPTO ONSET #1 with the average of each mouse in each row (eg 3 mice --> 3 rows)
for nummice=1:length(mice_list_virus.(TrialDay{v}))
    for v=1:length(TrialDay)
        for d=1:length(dFF_names)
            for pow = 1:length(Opto_Powers)
                for k = 1:size(datatype,2)
                    fullstreams_aligned_pooled.(TrialDay{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(nummice,:) = ...
                        fullstreams_aligned.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(1:length_updated);
                    
                end
            end
        end
    end
end

for pow = 1:length(Opto_Powers)
    time_vect_aligned.(Opto_Powers{pow}) = time_vect_aligned.(Opto_Powers{pow})(1:length_fullstreams_aligned_pooled-chop_off_ix+1);
end
% Create plot with nanmean and nanstd, Save graph --> NOTE TO GET ONE DRUG AFTER THE OTHER DO THIS IN GRAPHAD
color2plot = {'b','g','r','k','c','b','k','c'};

for k = 1; %:size(datatype,2) % dFF only k=1
    for d=1:length(dFF_names)
        for p = 2; %1:length(pooledtype) % baselinecorr only p=2            
            % figure 1: overlap
            figure; clf; hold on
            if show_plot == 0
                set(gcf,'visible','off')
            end
            for v=1:length(TrialDay)
                subplot(length(TrialDay),1,v)
                for pow = 1:length(Opto_Powers)
                    tmp_avg_pool = nanmean(fullstreams_aligned_pooled.(TrialDay{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),1); 
                    tmp_error_pool = nanstd(fullstreams_aligned_pooled.(TrialDay{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),1,1)./...
                    sqrt(sum(~isnan(fullstreams_aligned_pooled.(TrialDay{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(:,1))));
                    error_area(time_vect_aligned.(Opto_Powers{pow}),tmp_avg_pool,tmp_error_pool,color2plot{pow},0.25); %t_trials is the time vector                  
                end
            

            % Formatting
            ylim([-10 50]); %see parameters setting at start                     
            ax = gca;
            ax.FontSize = 20; 
            xlabel('Time (s)','FontSize', 20)
            ylabel('dFF','FontSize', 20)
            box off
            set(gca,'TickLength',[0 0])
            title(TrialDay{v});
%             L=legend(drug_name{1},' ',drug_name{2},' '); L.Location = 'Best';
            end
            
            %saveplot or not
            if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\streams_aligned_overlap',pooledtype{p},' ',datatype{k},'.tif']);
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\streams_aligned_overlap',pooledtype{p},' ',datatype{k},'.fig']);
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\streams_aligned_overlap',pooledtype{p},' ',datatype{k},'.pdf']);

            end          
        end
    end
end


%% Downsample for GraphPad plotting
for v=1:length(TrialDay)
    for d=1:length(dFF_names)
        for pow = 1:length(Opto_Powers)
            for k = 1; %:size(datatype,2) % here only dFF
                for nummice=1:size(fullstreams_aligned_pooled.(TrialDay{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),1)
                    data2plot=fullstreams_aligned_pooled.(TrialDay{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(nummice,:);
                    ttt = timetable(data2plot','SampleRate',sampling_rate_ds);
                    tttAvg = retime(ttt,'regular','mean','SampleRate',1); %1fps
%                     figure; plot(tttAvg.Time,tttAvg.Var1)
                    fullstreams_aligned_pooled_ds.(TrialDay{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(nummice,:) = tttAvg.Var1;
                end
                time_vect_aligned_ds.(Opto_Powers{pow}) = [seconds(tttAvg.Time)] + time_vect_aligned.(Opto_Powers{pow})(1); %also add the first timepoint from the previous timevector allowing to shift in time
                time_vect_aligned_ds_minutes.(Opto_Powers{pow}) = time_vect_aligned_ds.(Opto_Powers{pow})/60;
                if contains(PATH2SAVE,'SCH')
                    time_vect_aligned_ds_minutes.(Opto_Powers{pow}) = time_vect_aligned_ds_minutes.(Opto_Powers{pow}) - 6; % skip first 6 minutes so that starts at 0
                elseif contains(PATH2SAVE,'HdLight-Ctrl')
                    time_vect_aligned_ds_minutes.(Opto_Powers{pow}) = time_vect_aligned_ds_minutes.(Opto_Powers{pow}) - 7; % skip first 6 minutes so that starts at 0            
                end
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
            for v=1:length(TrialDay)
                subplot(length(TrialDay),1,v)
                for pow = 1:length(Opto_Powers)
                    tmp_avg_pool = nanmean(fullstreams_aligned_pooled_ds.(TrialDay{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),1); 
                    tmp_error_pool = nanstd(fullstreams_aligned_pooled_ds.(TrialDay{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),1,1)./...
                    sqrt(sum(~isnan(fullstreams_aligned_pooled_ds.(TrialDay{v}).(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(:,1))));
                    error_area(time_vect_aligned_ds.(Opto_Powers{pow}),tmp_avg_pool,tmp_error_pool,color2plot{pow},0.25); %t_trials is the time vector
                end

            % Formatting
            ylim([-10 50]); %see parameters setting at start                     
            ax = gca;
            ax.FontSize = 20; 
            xlabel('Time (s)','FontSize', 20)
            ylabel('dFF','FontSize', 20)
            box off
            set(gca,'TickLength',[0 0])          
%             L=legend(drug_name{1},' ',drug_name{2},' '); L.Location = 'Best';
            title(TrialDay{v});

            end

            %saveplot or not
            if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\streams_aligned_overlap_AVE1sec',pooledtype{p},' ',datatype{k},'.tif']);
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\streams_aligned_overlap_AVE10sec',pooledtype{p},' ',datatype{k},'.fig']);
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\streams_aligned_overlap_AVE1sec',pooledtype{p},' ',datatype{k},'.pdf']);

            end          
        end
    end
end



%% %% Individual trials- curves: ALL GROUPS TOGETHER ON SAME GRAPH, all animals
% Create plot with nanmean and nanstd, Save graph
color2plot = {'m','r','m','k','c','b','k','c'};
TrialType_Number = length(Opto_Powers);


% for k = 1:size(Body_parts,2)
for k = 1:size(datatype,2) %Zscore oly, k=1
    for d=1:length(dFF_names)
        for p = 2; %1:length(pooledtype) % baselinecorr only p=2
            figure; clf; hold on
            if show_plot == 0
                set(gcf,'visible','off')
            end
            
            for pow = 1:length(Opto_Powers)
                for v=1:length(TrialDay)
                    tmp_avg_pool = nanmean(PooledAVE.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1); 
                    tmp_error_pool = nanstd(PooledAVE.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1,1)./...
                    sqrt(sum(~isnan(PooledAVE.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(:,1))));
                    error_area(t_trials,tmp_avg_pool,tmp_error_pool,color2plot{v*(pow+1)},0.25); %t_trials is the time vector
                    max_avg(pow) = max(tmp_avg_pool);
                end
            end
            % plot properties
            PeakAve = num2str(mean(max_avg));
            PeakMax = num2str(max(max_avg));
%             sgtitle([(pooledtype{p}),' ',datatype{k},', AvePeak: ',PeakAve,', MaxPeak: ',PeakMax],'Interpreter','none')

            % Formatting
            xlim([TRANGE(1) TRANGE(2)]);
%             xlim([-30 60]);

            ylim([-5 50]); %see parameters setting at start                     
            ax = gca;
            ax.FontSize = 20; 
            xlabel('Time (s)','FontSize', 20)
            ylabel('dFF','FontSize', 20)
%             tickarray = [0:10:size(alltrials_array,1)];
%             tickarray(1) = 1;
%             set(gca,'YTick',tickarray);
%             title('   Shock','FontSize', 12);
            box off
            ylimits = get(gca,'YLim');
            xlimits = get(gca,'XLim');
            set(gca,'TickLength',[0 0])
            %             xticks([-30 0 30 60 90]);
            %plot the shock area
            error_area_onlyrectangle([0 2],ylimits,[ylimits(2)-ylimits(1) ylimits(2)-ylimits(1)],[0.0500, 0.0250, 0.0980],0.1); %(X,Y,barlength,color,alpha for shading of SEM). Barlength should be Ylimit_max - Ylimitmin
            plot([0 0],ylimits,'k','LineWidth',0.5,'LineStyle','-.')
            yline(0,'-.k');

            % Legend, self-made - with lines (OPTION 1)
%             for v=1:length(TrialDay)
%                 text4legend = TrialDay{v};
%                 text4legend(end-1:end)=[];
%                 plot([xlimits(1)+1.5 xlimits(1)+3],[(ylimits(2)-0.1)-v*0.8 (ylimits(2)-0.1)-v*0.8],'color',color2plot{v},'LineWidth',2,'LineStyle','-')
%                 text([xlimits(1)+3.5 xlimits(1)+3.5],[(ylimits(2)-0.1)-v*0.8 (ylimits(2)-0.1)-v*0.8],text4legend,'Color','black','FontSize',14)
%             end
            
            % Legend, self-made - with text only (OPTION 2)
            for pow = 1:length(Opto_Powers)
                for v=1:length(TrialDay)
                    text4legend = [TrialDay{v},' ',Opto_Powers{pow}];
    %                 text4legend(1:4)=[];
                    text4legend = [text4legend]
    %                 plot([xlimits(1)+1.5 xlimits(1)+3],[(ylimits(2)-0.1)-v*0.8 (ylimits(2)-0.1)-v*0.8],'color',color2plot{v},'LineWidth',2,'LineStyle','-')
                    text([xlimits(1)+30*pow xlimits(1)+30*pow],[(ylimits(2)-0.1)-v*4 (ylimits(2)-0.1)-v*4],text4legend,'Color',color2plot{v*(pow+1)},'FontSize',14)
                end
            end


            %saveplot or not
            if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\AVE_alldays_alltrials_plots',pooledtype{p},' ',datatype{k},'.tif']);
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\AVE_alldays_alltrials_plots',pooledtype{p},' ',datatype{k},'.fig']);
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\AVE_alldays_alltrials_plots',pooledtype{p},' ',datatype{k},'.pdf']);

            end
        end
    end
end
% end


%% %% Individual trials- curves: ALL GROUPS TOGETHER ON SAME GRAPH - ONLY MOUSE OF CHOICE (representative trial for each drug)
%Edit nummice to mouse of choice below
nummiceX = 1; %3
s=1; % session number
TrialNumber2Show=2; % trial number or representative figure

% Create plot with nanmean and nanstd, Save graph
color2plot = {'m','r','m','k','c','b','k','c'};
TrialType_Number = length(Opto_Powers);

% for k = 1:size(Body_parts,2)
for k = 2; %:size(datatype,2) %
    for d=1:length(dFF_names)
        for p = 2; %1:length(pooledtype) % baselinecorr only p=2
            figure; clf; hold on
            if show_plot == 0
                set(gcf,'visible','off')
            end
            for pow = 1:length(Opto_Powers)
                for v=1:length(TrialDay)
                    plot(t_trials,PooledINDIV.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummiceX}).(SessionIDs{s})(TrialNumber2Show,:),'Color',color2plot{v*(pow+1)}); hold on               
                end
            end
            
            % Formatting
            xlim([TRANGE(1) TRANGE(2)]);
%             ylim(limits2plot.(pooledtype{p}).(dFF_names{d}).(datatype{k})); %see parameters setting at start               
            ylim([-5 15])
            ax = gca;
            ax.FontSize = 20; 
            xlabel('Time (s)','FontSize', 20)
            ylabel('dFF','FontSize', 20)
%             tickarray = [0:10:size(alltrials_array,1)];
%             tickarray(1) = 1;
%             set(gca,'YTick',tickarray);
%             title('   Shock','FontSize', 12);
            box off
            ylimits = get(gca,'YLim');
            xlimits = get(gca,'XLim');
            set(gca,'TickLength',[0 0])
            %             xticks([-30 0 30 60 90]);
            %plot the shock area
            error_area_onlyrectangle([0 2],ylimits,[ylimits(2)-ylimits(1) ylimits(2)-ylimits(1)],[0.0500, 0.0250, 0.0980],0.1); %(X,Y,barlength,color,alpha for shading of SEM). Barlength should be Ylimit_max - Ylimitmin
            plot([0 0],ylimits,'k','LineWidth',0.5,'LineStyle','-.')
            yline(0,'-.k');

            % Legend, self-made - with lines (OPTION 1)
%             for v=1:length(TrialDay)
%                 text4legend = TrialDay{v};
%                 text4legend(end-1:end)=[];
%                 plot([xlimits(1)+1.5 xlimits(1)+3],[(ylimits(2)-0.1)-v*0.8 (ylimits(2)-0.1)-v*0.8],'color',color2plot{v},'LineWidth',2,'LineStyle','-')
%                 text([xlimits(1)+3.5 xlimits(1)+3.5],[(ylimits(2)-0.1)-v*0.8 (ylimits(2)-0.1)-v*0.8],text4legend,'Color','black','FontSize',14)
%             end
            
            % Legend, self-made - with text only (OPTION 2)
            for pow = 1:length(Opto_Powers)

            for v=1:length(TrialDay)
                    text4legend = [TrialDay{v},' ',Opto_Powers{pow}];
%                 text4legend(1:4)=[];
                text4legend = [text4legend]
%                 plot([xlimits(1)+1.5 xlimits(1)+3],[(ylimits(2)-0.1)-v*0.8 (ylimits(2)-0.1)-v*0.8],'color',color2plot{v},'LineWidth',2,'LineStyle','-')
                    text([xlimits(1)+30*pow xlimits(1)+30*pow],[(ylimits(2)-0.1)-v*1.4 (ylimits(2)-0.1)-v*1.4],text4legend,'Color',color2plot{v*(pow+1)},'FontSize',14)
            end
            end


            %saveplot or not
            if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\AVE_alldays_alltrials_plots_specificmouse',pooledtype{p},' ',datatype{k},'.tif']);
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\AVE_alldays_alltrials_plots_specificmouse',pooledtype{p},' ',datatype{k},'.fig']);
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\AVE_alldays_alltrials_plots_specificmouse',pooledtype{p},' ',datatype{k},'.pdf']);

            end
        end
    end
end
% end


%% Create a giant heatmaps with the average for each mouse - all mice included
orderedmap = 0; % 0 to plot in normal order, 1 if you want to plot in the order of deg inhibition by sorting on the minimal value between 3 and 10 sec (for 10 sec stim...)
% OVERALL HEATMAP
for v=1:length(TrialDay)
    for p = 2; %1:length(pooledtype)
        for d=1:length(dFF_names)
            for k = 2; %1:size(datatype,2)
                if k == 1
                    Color_scale = [-1 8];
                else
                    Color_scale = [-2 4];
                end
                figure;
                if show_plot == 0
                    set(gcf,'visible','off')
                end
                for pow = 1:length(Opto_Powers)
                    subplot(length(Opto_Powers),1,pow)
                    % mean of individual trials for each mouse- into a heatmap that includes all mice
                    mean_array = PooledAVE.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k});
%                     subplot(length(IdChannel)+1,1,o)
                    % if you want mice sorted by order in terms of deg inhibition
                    if orderedmap == 1
                        [~,order] = sort(abs(min(mean_array(:,t_trials > 2 & t_trials < stim_duration),[],2))); %if A is a matrix, then min(A,[],2) is a column vector containing the minimum value of each row
                        %B = sort(A) sorts the elements of A in ascending order + [B,I] = sort(___) also returns a collection of index vectors
                    else
                        numrows = size(mean_array,1);
                        order = [1:numrows]';
                    end
%                     imagesc(t_trials,1,mean_array(order,:))
                    imagesc(t_trials,1,mean_array)
                    colormap('parula')
                    c = colorbar;
                    c.Label.String = 'dFF';
                    if ~isempty(Color_scale)
                        c_limits = Color_scale;
                        caxis(c_limits)
                    end
                    hold on
                    xlim([t_trials(1) t_trials(end)])
                    xlabel('Time (s)')
                    ylabel('Animal')
                    set(gca,'YTick',1:size(mean_array,1));
                    title([Opto_Powers{pow}]);
                    box off
                    ylimits = get(gca,'YLim');
                    plot([0 0],ylimits,'k','LineWidth',2)
                    plot([stim_duration stim_duration],ylimits,'k','LineWidth',2,'LineStyle',':')

                    %     xline(0,'-k');
%                     plot(round(t_Dipper(1,2)-t_Dipper(1,1)),ylimits,'m','LineWidth',2)
                end
                sgtitle([TrialDay{v},', ',datatype{k},', ',pooledtype{p},', ',dFF_names{d}]);
                if save_plot == 1
                    saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\Heatmaps of dFF signals',pooledtype{p},' ',datatype{k},'.tif']);
                    saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\Heatmaps of dFF signals',pooledtype{p},' ',datatype{k},'.fig']);
                end
            end
        end
    end
end






