% Marie_FP_PooledData
% Marie Labouesse, marie.labouesse@gmail.com - Dec 2023
% Cohort: NAc dLight1.3b with opto-inhibition (eNPHR)

% 1- POOLED DATA
% load matlab spaces generated in: Marie_FP_IndivData_extraction
% puts individual trials into one big matrix (all mice together): PooledINDIV
% calculates the mean of individual trial trajectories aligned to specific event for each mouse: PooledAVE
% saves a workspace that can then be used for quantifications

% 2- POOLED GRAPHS
% generates average graphs curves (this used to be in the main script)
% generates heatmaps (1 average trace/mouse or all trials/all mice in one big heatmap - separated by power)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%% INITIALIZATIONS
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% PARAMETERS 
show_plot = 1; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
reanalysis = 0; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed. This is based on the existence of the matlab space "IndividualData.mat"
overwrite = 0; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not          
Collapse_sessions = 0; % If 0, analysis only of individual sessions. If 1, analysis of the group data (pooled)          

pooledtype = {'raw','baselinecorr'};

%% Define the path where the data is
PATH2DATA_0 = uigetdir('select folder'); %select the overall folder
virus = {'eNPHR3'}; 


for v=1:length(virus)
    PATH2DATA.(virus{v}) = [PATH2DATA_0,'\',virus{v}];
    PATH2SAVEFOLDER.(virus{v}) = [PATH2DATA_0,'\',virus{v}];
    mice_list_virus.(virus{v}) = dir(PATH2DATA.(virus{v})); %all things in this folder
    mkdir([PATH2SAVEFOLDER.(virus{v}),'\pooled figures\']);
    mkdir([PATH2SAVEFOLDER.(virus{v}),'\pooled data\']);
end

                

%% IDENTIFY MICE TO ANALYZE 
for v=1:length(virus)
    for o = length(mice_list_virus.(virus{v})):-1:1
        if mice_list_virus.(virus{v})(o).isdir == 0  %remove non-folders
            mice_list_virus.(virus{v})(o) = [];
        else
            if  strcmp(mice_list_virus.(virus{v})(o).name,'data') == 1 || strcmp(mice_list_virus.(virus{v})(o).name,'figures') == 1 ...   %remove folders with data or figures or ., ..
                || contains(mice_list_virus.(virus{v})(o).name,'data') || contains(mice_list_virus.(virus{v})(o).name,'figures') ...
                || contains(mice_list_virus.(virus{v})(o).name,'results') || contains(mice_list_virus.(virus{v})(o).name,'turns') || contains(mice_list_virus.(virus{v})(o).name,'other') ...
                || strcmp(mice_list_virus.(virus{v})(o).name,'.') == 1 || strcmp(mice_list_virus.(virus{v})(o).name,'..') == 1
                mice_list_virus.(virus{v})(o) = [];
            end
        end
    end
    Nmice_virus{v} = length(mice_list_virus.(virus{v}));

end

              

%% Import individual workspaces and load the data into a pooled array
for v=1:length(virus)
    %AnimalIDs
    for nummice=1:length(mice_list_virus.(virus{v}))
        AnimalIDs(nummice) = {['ID',mice_list_virus.(virus{v})(nummice).name(8:11)]};
    end   
    
    for nummice=1:length(mice_list_virus.(virus{v}))
        % Define the path to the mouse and find the folders to analyze: 
        path2mouse = [PATH2DATA.(virus{v}),'\',mice_list_virus.(virus{v})(nummice).name,'\']; 
        % if there are several sessions inside the mouse folder
        sessions = dir(path2mouse);
        %remove non relevant folders
        for o = length(sessions):-1:1
            if sessions(o).isdir == 0  %remove non-folders
                sessions(o) = [];
            elseif strcmp(sessions(o).name,'data') == 1 || strcmp(sessions(o).name,'figures') == 1 || strcmp(sessions(o).name,'results') == 1 ...
                 || strcmp(sessions(o).name,'.') == 1 || strcmp(sessions(o).name,'..') == 1 || contains(sessions(o).name,'.csv') == 1 ...
                 || contains(sessions(o).name,'figures') == 1 || contains(sessions(o).name,'data') == 1 || contains(sessions(o).name,'other')  || contains(sessions(o).name,'BACKUP')
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
        if exist([PATH2SAVEFOLDER.(virus{v}),'\',mice_list_virus.(virus{v})(nummice).name,'\'],'dir') == 0
            mkdir([PATH2SAVEFOLDER.(virus{v}),'\',mice_list_virus.(virus{v})(nummice).name,'\']);
        end
        path2save_mouse = [PATH2SAVEFOLDER.(virus{v}),'\',mice_list_virus.(virus{v})(nummice).name,'\'];
        
        %% Loop for all the sessions for the mouse: Sessions are experimental days you want to average together (replicates of each other) ............. 
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
            done = exist([PATH2SAVE,'PooledAllMice.mat'],'file'); % .................. 
            if done == 0 || reanalysis == 1 %if not done or if wish to reanalyze
        
                % load the mouse/session matlab space
                load([PATH2SESSION,'\IndividualData.mat']);
       
                % Initialization of pooled structure with individual trials   %................ ..

                if nummice == 1 && s == 1
                    for d=1:length(dFF_names)
                        for pow = 1:length(Opto_Powers)
                            for k = 1:size(datatype,2)
                                % ANIMAL ID .................... add session ID later if needed
                                PooledAnimalID.(virus{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = {};
                                % INDIVIDUAL TRIALS OF EACH MICE FOR ALL MICE
                                PooledINDIV.(virus{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = ones(size(IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})))*nan;  
                                PooledINDIV.(virus{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = ones(length(mice_list_virus.(virus{v})),size(t_trials,2))*nan;  
                                % AVERAGE TRIAL OF EACH MICE (ALL SESSIONS) FOR ALL MICE, EACH MOUSE ON ONE LINE
                                PooledAVE.(virus{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(virus{v})),size(t_trials,2))*nan;  
                                PooledAVE.(virus{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(virus{v})),size(t_trials,2))*nan; 
                                % STREAMS                
                                fullstreams.(virus{v}).(AnimalIDs{nummice}) = {};

                                % STREAMS ALIGNED                
                                fullstreams_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) = {};
                                % TIME VECT ALIGNED
                                timevect_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) = {}; 

     
                            end
                        end
                    end
                end
                
                % Add data one animal at a time
                % ANIMAL ID and SESSION ID
                PooledAnimalID.(virus{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = SessionNames(s);
                PooledAnimalID_only.(virus{v}) = AnimalIDs;
                % STORE INDIVIDUAL TRIALS
                for d=1:length(dFF_names)
                    for pow = 1:length(Opto_Powers)
                        for k = 1:size(datatype,2)
                            PooledINDIV.(virus{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}); 
                            PooledINDIV.(virus{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}); 
                        end
                    end
                end                            
                % STORE AVERAGE AND SELECT TRIALS
                trials_select_VEH = [3:10]; %there are 10 trials, VEH injection is at 0 min, every 3 min
                trials_select_DETQ = [4:11]; % there are 11 trials, DETQ injection is at 0 min, every 3 min

                for d=1:length(dFF_names) %NAc
                    for pow = 1:length(Opto_Powers) %{'VEH'}    {'DETQ'}
                        for k = 1:size(datatype,2)%{'dFF'}    {'dFFwithin'}    {'ZScoredFF_within'}    {'ZScoredFF'}
                            if length(sessions) == 1

                                if pow == 1
                                    PooledAVE.(virus{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = nanmean(IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(trials_select_VEH,:),1); 
                                    PooledAVE.(virus{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = nanmean(IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(trials_select_VEH,:),1); 
                                elseif pow == 2
                                    PooledAVE.(virus{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = nanmean(IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(trials_select_DETQ,:),1); 
                                    PooledAVE.(virus{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(nummice,:) = nanmean(IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(trials_select_DETQ,:),1); 

                                end
                            end
                        end
                    end
                end

                for d=1:length(dFF_names)
                    for pow = 1:length(Opto_Powers)
                        for k = 1; %:size(datatype,2)
                            %STORE STREAMS
                            fullstreams.(virus{v}).(AnimalIDs{nummice}) = streams;                       
            
                            %STORE STREAMS ALIGNED TO OPTO ONSET #1
                            fullstreams_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) = streams_aligned.(datatype{k}).(dFF_names{d}).(Opto_Powers{pow});
                            timevect_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) = time_vect_aligned.(Opto_Powers{pow}); % DETQ/VEH or DRUG levels

                        end
                    end
                end
            end
        end
    end
end





%% Save pooled data for all viruses together
% STORE INDIVIDUAL and AVERAGE TRIALS
save([PATH2DATA_0,'\PooledAllMice.mat'],'PooledINDIV','PooledAVE','virus','PooledAnimalID','PooledAnimalID_only','AnimalIDs','dFF_names','datatype','Opto_Powers',...
    't_trials','stim_duration','TRANGE','BASELINE_WIN','Events2Align2','dt_ds','sampling_rate_ds','length_data','fullstreams','time_vect','pooledtype');




%% Create plot with nanmean and nanstd, Save graph

color2plot_triplet={[0 0 0]/256,[255 128 128]/256};

%plot
TrialType_Number = length(Opto_Powers);

for v=1:length(virus)
    for k = 2; %1:size(datatype,2) %dff within
        for d=1:length(dFF_names)
            for p = 2; %1:length(pooledtype) %raw
                figure; clf; hold on
                if show_plot == 0
                    set(gcf,'visible','off')
                end
                for pow = 1:length(Opto_Powers)
                    tmp_avg_pool = nanmean(PooledAVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1); 
                    tmp_error_pool = nanstd(PooledAVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1,1)./...
                    sqrt(sum(~isnan(PooledAVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(:,1))));
                    error_area(t_trials,tmp_avg_pool,tmp_error_pool,color2plot_triplet{pow},0.25); %t_trials is the time vector
                end

                % individual traces
                % hold on
                % for kui = 1:size(PooledAVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1)
                %     plot(t_trials,PooledAVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(kui,:))
                % end

                % plot properties
                xline(0,'-k');
                xline(stim_duration,'--k');
                yline(0,'-.k');
                %     yline(max(max_avg),'-.r');
                xlabel('Time (s)','FontSize', 20);
                ylabel('\DeltaF/F (%)','FontSize', 20);
                xlim([-30,30]);
                ylim([-3,3])
                xticks([-30 -20 -10 0 10 20 30])
                ax = gca;
                ax.FontSize = 20; 
                % sgtitle([(pooledtype{p}),' ',datatype{k}],'Interpreter','none')

                annotation('textbox', [0.16, 0.90, 0.2, 0.05], 'String', 'Veh', 'Color', 'black', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none');
                annotation('textbox', [0.16, 0.84, 0.2, 0.05], 'String', 'DETQ', 'Color', [255 128 128]/256, 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none');

                % Create line foir opto stim
                annotation('line',[0.52 0.645],...
                    [0.95 0.95],'Color',[0.78 0.02 0.02],'LineWidth',10);

               
                %saveplot or not
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVEFOLDER.(virus{v}),'\pooled figures\AVE all powers optostim',pooledtype{p},' ',datatype{k},'.tif']);
                    saveas(gcf,[PATH2SAVEFOLDER.(virus{v}),'\pooled figures\AVE all powers optostim',pooledtype{p},' ',datatype{k},'.fig']);
                    saveas(gcf,[PATH2SAVEFOLDER.(virus{v}),'\pooled figures\AVE all powers optostim',pooledtype{p},' ',datatype{k},'.pdf']);
                end
            end
        end
    end
    % end
end


%% Quantify

% WINDOW DURING WHICH TO CALCULATE AUC AND OTHER VARIABLES: 
windowcalc_duration = 10; %10 seconds; for the opto calculations.
prestim_epoch_start = -10;
prestim_epoch_stop = -1; 
prestim_duration = -(prestim_epoch_start - prestim_epoch_stop); % 9 seconds, ie -10 to -1 sec
poststim_epoch_start = 10;
poststim_epoch_stop = 20;

%% Compute the AUC and the Peaks
% Determine the indexes where to calculate the AUC and other variables in the stim epoch
temp_t = time_vect - time_vect(1); %values in timevect shifted by one value. We dont use the t_trials cos it doesnt start at 0
dummie = 1:length(temp_t); % indexes of time vector
dummie2 = dummie(temp_t >= windowcalc_duration); % so look for the indexes when the time vector is superior or equal to the stim duration
idx_AUC2 = dummie2(1); %find the first index when the time vector hit the exact value of the end of the stim duration; this is the second AUC index, or the number of indixes you need to go thru to have the stim duration
idx_AUC1 = find(t_trials == 0); % this is the index in t_trials where the stimulation starts
idx_AUC = idx_AUC1:(idx_AUC1+(idx_AUC2-1)); % we will be calculating the AUC from the beginning of the stimulation to the index just before the first index at the end of the stim
% Determine the indexes where to calculate the AUC and other variables in the prestim epoch
dummie3 = dummie(temp_t >= prestim_duration); % so look for the indexes when the time vector is superior or equal to the epoch
idx_AUC3 = dummie3(1); %find the first index when the time vector hit the exact value of the end of the epoch; this is the second AUC index, or the number of indixes you need to go thru to have the epoch duration
dummie4 = find(t_trials >= prestim_epoch_start); 
idx_AUC4 = dummie4(1); % this is the index in t_trials where the prestimulation epoch starts
idx_AUC_pre = idx_AUC4:(idx_AUC4+(idx_AUC3-1)); % we will be calculating the AUC from the beginning of the stimulation to the index just before the first index at the end of the stim
clear dummie dummie2 dummie3 dummie4 idx_AUC1 idx_AUC2 idx_AUC3 idx_AUC4
% Determine the indexes where to calculate the AUC and other variables in the post- stim epoch
temp_t = time_vect - time_vect(1); %values in timevect shifted by one value. We dont use the t_trials cos it doesnt start at 0
dummie5 = 1:length(temp_t); % indexes of time vector
dummie6 = dummie5(temp_t >= poststim_epoch_stop); % so look for the indexes when the time vector is superior or equal to the stim duration
idx_AUC5 = dummie6(1); %find the first index when the time vector hit the exact value of the end of the stim duration; this is the second AUC index, or the number of indixes you need to go thru to have the stim duration
dummie7 = dummie5(temp_t >= poststim_epoch_start); % so look for the indexes when the time vector is superior or equal to the stim duration
idx_AUC6 = dummie7(1); % this is the index in t_trials where the stimulation starts
idx_AUC_post = idx_AUC6:(idx_AUC6+(idx_AUC5-1)); % we will be calculating the AUC from the beginning of the stimulation to the index just before the first index at the end of the stim

%% Initialize the "Measurements" structure where you will be putting the AUC, Peaks and specific values. 
% Here I decided to save all the mice at the same place. That means I need to run this everytime again when I add mice. But also means this can be independent from the initial dataloading
% And accessing the data is easy; since now all the trials come from one matlab space for all groups.

for v=1:length(virus) % which group
    for p = 1:length(pooledtype)   % raw or baselinecorr
        for d=1:length(dFF_names)   % if different recordings eg from 2 different brain regions
            for pow = 1:length(Opto_Powers) % different trial types; eg different opto powers or different drugs
                for k = 1:size(datatype,2) % ZScoredFF and dFF and dFF-within
                    for nummice=1:length(AnimalIDs)    % all mice
                            t_opto_i = PooledAVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k});
                            Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).AUC = ones(size(t_opto_i,1),1)*nan; %same number of rows as trials
                            Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).AbsMaxima = ones(size(t_opto_i,1),1)*nan; %value of maxima
                            Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).AbsMinima = ones(size(t_opto_i,1),1)*nan; %value of minima
                            Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).DegChange = ones(size(t_opto_i,1),1)*nan; %degree of change (onset vs maxima) 
                            Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).Average = ones(size(t_opto_i,1),1)*nan; %value of minima
                            Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).AbsMinima_pre = ones(size(t_opto_i,1),1)*nan; %value of minima
                            Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).Average_pre = ones(size(t_opto_i,1),1)*nan; %value of minima
                            Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).AbsMinima_post = ones(size(t_opto_i,1),1)*nan; %value of minima
                            Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).Average_post = ones(size(t_opto_i,1),1)*nan; %value of minima

                    end
                end
            end
        end
    end
end

%% Determine the AUC for entire section (anything below 0 will count for negative AUC) ; and other variables
for v = 1:length(virus) 
    for p = 1:length(pooledtype)   
        for d = 1:length(dFF_names)   
            for pow = 1:length(Opto_Powers) 
                for k = 1:size(datatype,2)
                    for nummice=1:length(AnimalIDs)
                            %%
                          
                            t_opto_i = PooledAVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k});
                            for w = 1:size(t_opto_i,1)   % = how many rows
                                %AUC, includes positive and negative values (no need to normalized cos fixed duration: 
                                dff_all = t_opto_i(w,:); % put the trace of this row in a new variable cos easier to handle
                                dff_stim = t_opto_i(w,idx_AUC);   %piece of dFF in the stim epoch we are interested in for this row
                                dff_prestim = t_opto_i(w,idx_AUC_pre); %piece of dFF in the pre-stim epoch we are interested in for this row
                                dff_poststim = t_opto_i(w,idx_AUC_post); %piece of dFF in the pre-stim epoch we are interested in for this row
                                tmp_AUC = trapz(dff_stim);
                                Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).AUC(w) = tmp_AUC; %same number of rows as trials

                                % Maxima in stim epoch  
                                AbsMaxima = max(dff_stim);
                                Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).AbsMaxima(w,:) = [AbsMaxima];                                
                                
                                % Minima in stim epoch   
                                AbsMinima = min(dff_stim);
                                Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).AbsMinima(w,:) = [AbsMinima];                                

                                % Average in stim epoch
                                Average = nanmean(dff_stim);
                                Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).Average(w,:) = [Average];                                
                               
                                % Degree change in stim epoch 
                                dFFonset = dff_stim(1);
                                % dFFoffset = t_opto_i(w,idx_AUC(end));
                                DegChange = 100* (AbsMaxima - dFFonset)./dFFonset;
                                Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).DegChange(w,:) = [DegChange]; %value at start, end and degree of change
                             
                                % Minima in pre epoch   
                                AbsMinima_pre = min(dff_prestim);
                                Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).AbsMinima_pre(w,:) = [AbsMinima_pre]; 

                                % Average in pre epoch
                                Average_pre = nanmean(dff_prestim);
                                Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).Average_pre(w,:) = [Average_pre];                                

                                % Minima in post epoch   
                                AbsMinima_post = min(dff_poststim);
                                Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).AbsMinima_post(w,:) = [AbsMinima_post]; 

                                % Average in post epoch
                                Average_post = nanmean(dff_poststim);
                                Measurements_AVE.(virus{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).Average_post(w,:) = [Average_post];                                

                            end
                            %%
                    end
                end
            end
        end
    end
end





%% Streams aligned graph
%% Adjust all to VEH 0 as baseline
for v=1:length(virus)    
    for nummice=1:length(mice_list_virus.(virus{v}))
        for d=1 %:length(dFF_names)
            for k = 1 %:size(datatype,2)
                for pow = 1:length(Opto_Powers)
                    fullstreams_aligned_norm.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) = ...
                        fullstreams_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) ...
                        - nanmean(fullstreams_aligned.(virus{v}).VEH.(AnimalIDs{nummice}))                    
                end
            end
        end
    end
end

%if happy
fullstreams_aligned = fullstreams_aligned_norm;

%% Streams aligned pooled
   
    % create empty pooled array
    nummice = 1; k=2; d=1; pow = 1; s=1; 
    for v=1:length(virus)    
        for d=1:length(dFF_names)
            for pow = 1:length(Opto_Powers)
                for k = 1:size(datatype,2)
                    length_fullstreams_aligned.(Opto_Powers{pow}) = ...
                    size(find(~isnan(fullstreams_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}))),2)       

                    fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}) = ...
                        ones(length(mice_list_virus.(virus{v})),length_fullstreams_aligned.(Opto_Powers{pow}))*nan;
                end
            end
        end

        %STORE STREAMS ALIGNED TO EVENT ONSET #1 with the average of each mouse in each row (eg 3 mice --> 3 rows)
        for nummice=1:length(mice_list_virus.(virus{v}))
            for d=1:length(dFF_names)
                for pow = 1:length(Opto_Powers)
                    for k = 1:size(datatype,2)
                        fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow})(nummice,1:length_fullstreams_aligned.(Opto_Powers{pow})) = ...
                            fullstreams_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice})(1:length_fullstreams_aligned.(Opto_Powers{pow}));
                          
                    end
                end
            end
        end
    end


%% Smooth the streams aligned

for v=1:length(virus)    
    for d=1:length(dFF_names)
        for pow = 1:length(Opto_Powers)
            for k = 1:size(datatype,2)
                for nummice=1:size(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}),1)
                                  
                    fullstreams_aligned_pooled_smooth.(virus{v}).(Opto_Powers{pow})(nummice,:) = smooth(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow})(nummice,:)',1000); %smooth 10 points at 100fps means smooth over 1 sec 
                
                end
                
            end
        end
    end
end

%% Downsample the streams aligned (smoothed)
for v=1:length(virus)    
    for d=1:length(dFF_names)
        for pow = 1:length(Opto_Powers)
            for k = 1:size(datatype,2)
                for nummice=1:size(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}),1)
                
                    data2plot=fullstreams_aligned_pooled_smooth.(virus{v}).(Opto_Powers{pow})(nummice,:);
                    ttt = timetable(data2plot','SampleRate',sampling_rate_ds);
                    tttAvg = retime(ttt,'regular','mean','SampleRate',1); %1fps
                    fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(nummice,:) = tttAvg.Var1;
                end
                
                time_vect_aligned_ds.(Opto_Powers{pow}) = [seconds(tttAvg.Time)] + timevect_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice})(1); %also add the first timepoint from the previous timevector allowing to shift in time                
                time_vect_aligned_ds_minutes.(Opto_Powers{pow}) = time_vect_aligned_ds.(Opto_Powers{pow})/60;
                
            end
        end
    end
end



    
%% Plot 
color2plot_triplet={[0 0 0]/256,[255 128 128]/256};
figure; clf; hold on
if show_plot == 0
    set(gcf,'visible','off')
end

idx_selct = [1,2,3,4,5]; 


for pow = 1:length(Opto_Powers)
    for v=1:length(virus)
        % average plot
        tmp_avg_pool = nanmean(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,:),1); 
        tmp_error_pool = nanstd(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,:),1,1)./...
        sqrt(sum(~isnan(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,1))));
        error_area(time_vect_aligned_ds_minutes.(Opto_Powers{pow}),tmp_avg_pool,tmp_error_pool,color2plot_triplet{pow},0.25); %t_trials is the time vector

    end
end

% Formatting
ylim([-30 80]); 
xlim([0 60]);
xticks([0 10 20 30 40 50 60]);
ax = gca;
ax.FontSize = 20; 
xlabel('Time (min)','FontSize', 20)
ylabel('\DeltaF/F (%)','FontSize', 20)
box off
set(gca,'TickLength',[0 0])

xline(3,':k','Veh','FontSize',20); % grey out the noise in signal due to injection 
xline(32,':r','DETQ','FontSize',20); % grey out the noise in signal due to injection 
% % xline(20,':k'); xline(21,':k','J60 or Sal','FontSize',20); 
error_area_onlyrectangle([0 3],[0 120],[200 200],[165, 167, 168]./255,1); 
error_area_onlyrectangle([29 32],[0 120],[200 200],[165, 167, 168]./255,1); 

% optostim. they start from +60 sec after the data starts; then every 3minutes (that's how the % data was saved)
%NB 
% trials_select_VEH = [3:10]; %there are 10 trials, VEH injection is at 0 min, every 3 min
for oiu = 1:3:10*3
    if oiu == 1
        % plot([oiu+1 oiu+1],[-20 -15],'Color',[0.78 0.02 0.02],'LineWidth',4)
    elseif oiu == 4
        % plot([oiu+1 oiu+1],[-20 -15],'Color',[0.78 0.02 0.02],'LineWidth',4)
    else
        plot([oiu+1 oiu+1],[-20 -15],'Color',[0.78 0.02 0.02],'LineWidth',4)
    end
end
% trials_select_DETQ = [4:11]
for oiu = 1:3:10*3
    if oiu == 1
        % plot([oiu+32 oiu+32],[-20 -15],'Color',[0.78 0.02 0.02],'LineWidth',4)
    elseif oiu == 4
        % plot([oiu+32 oiu+32],[-20 -15],'Color',[0.78 0.02 0.02],'LineWidth',4)
    elseif oiu == 4
        % plot([oiu+32 oiu+32],[-20 -15],'Color',[0.78 0.02 0.02],'LineWidth',4)
    else
        plot([oiu+32 oiu+32],[-20 -15],'Color',[0.78 0.02 0.02],'LineWidth',4)
    end
end

%saveplot or not
if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVEFOLDER.(virus{v}),'\pooled figures\streams_aligned_plot',pooledtype{p},' ',datatype{k},'.tif']);
                    saveas(gcf,[PATH2SAVEFOLDER.(virus{v}),'\pooled figures\streams_aligned_plot',pooledtype{p},' ',datatype{k},'.fig']);
                    saveas(gcf,[PATH2SAVEFOLDER.(virus{v}),'\pooled figures\streams_aligned_plot',pooledtype{p},' ',datatype{k},'.pdf']);
end       