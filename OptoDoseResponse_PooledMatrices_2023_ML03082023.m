% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT
% Name: OptoDoseResponse_PooledMatrices_2023
% Cohort: PFC dLight1.3b with opto 

% May 2023 - Marie Labouesse, marie.labouesse@gmail.com

% 1- POOLED DATA
% loads FP data previously analyzed (with: OptoDoseResponse_DataExtraction) in matlab space "IndividualData.mat" from multiple animals (select folder containing multiple animals) 
% puts individual trials into one big matrix (all mice together): PooledINDIV
% calculates the mean of individual trial trajectories aligned to specific event for each mouse: PooledAVE
% saves a workspace that can then be used for quantifications (OptoDoseResponse_PooledQuantification) or graphs OptoDoseResponse_PooledGraphs

% edit relevant parameters in %% SETUP PARAMETERS

% functions needed: error_area_onlyrectangle


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%% INITIALIZATIONS
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1]);

%% PARAMETERS 
show_plot = 1; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
reanalysis = 1; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed. 
overwrite = 0; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not          
Collapse_sessions = 0; % If 0, analysis only of individual sessions. If 1, analysis of the group data (pooled)   

pooledtype = {'raw','baselinecorr'};

%% Define the path where the data is
PATH2DATA_0 = uigetdir('select folder'); %select the overall folder
TrialDay = {'VEH','DETQ'}; 
PATH2SAVEPOOL = PATH2DATA_0;
mkdir([PATH2SAVEPOOL,'\pooled figures\']);

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
            if  strcmp(mice_list_virus.(TrialDay{v})(o).name,'data') == 1 || strcmp(mice_list_virus.(TrialDay{v})(o).name,'figures') == 1 ...   
                || contains(mice_list_virus.(TrialDay{v})(o).name,'data') || contains(mice_list_virus.(TrialDay{v})(o).name,'figures') ...
                || contains(mice_list_virus.(TrialDay{v})(o).name,'results') || contains(mice_list_virus.(TrialDay{v})(o).name,'other')...
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
        AnimalIDs(nummice) = {['ID',mice_list_virus.(TrialDay{v})(nummice).name(end-4:end)]};
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
                 || strcmp(sessions(o).name,'.') == 1 || strcmp(sessions(o).name,'..') == 1 || contains(sessions(o).name,'.csv') == 1 || contains(sessions(o).name,'other') == 1
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
            done = exist([PATH2SAVE,'PooledAllMice.mat'],'file'); % 
            if done == 0 || reanalysis == 1 %if not done or if wish to reanalyze
        
                % load the mouse/session matlab space
                load([PATH2SESSION,'IndividualData.mat']);
                Opto_Powers = {'DRUG'};
                if exist('IndivStim_data_VEH','var')   % 2022                 
                    IndivStim_data = IndivStim_data_VEH;
                    streams_aligned = streams_aligned_VEH;
                    time_vect_aligned = time_vect_aligned_VEH;
                    time_optotrials = time_optotrials_double.VEH;
                    clear IndivStim_data_VEH streams_aligned_VEH time_vect_aligned_VEH
                    
                elseif exist('IndivStim_data_DETQ','var')  % 2022                  
                    IndivStim_data = IndivStim_data_DETQ;
                    streams_aligned = streams_aligned_DETQ;
                    time_vect_aligned = time_vect_aligned_DETQ;
                    time_optotrials = time_optotrials_double.DETQ;                    
                    clear IndivStim_data_DETQ streams_aligned_DETQ time_vect_aligned_DETQ
                
                else % 2021
                    time_optotrials = time_optotrials_single;
                    clear time_optotrials_single
                end
                
                % Initialization of pooled structure with individual trials   %

                if nummice == 1 && s == 1
                    for d=1:length(dFF_names)
                        for pow = 1:length(Opto_Powers)
                            for k = 1:size(datatype,2)
                                % ANIMAL ID .................... 
                                PooledAnimalID.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = {};
                                % INDIVIDUAL TRIALS OF EACH MICE FOR ALL MICE
                                PooledINDIV.(TrialDay{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = ones(size(IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})))*nan;  
                                PooledINDIV.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan;  
                                % AVERAGE TRIAL OF EACH MICE (ALL SESSIONS) FOR ALL MICE, EACH MOUSE ON ONE LINE
                                PooledAVE.(TrialDay{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan;  
                                PooledAVE.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan; 
                                 % AVERAGE TRIAL OF EACH MICE (ALL SESSIONS) FOR ALL MICE, EACH MOUSE ON ONE LINE VEH NORMALIZED
                                PooledAVE_VEHnormalized.(TrialDay{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan;  
                                PooledAVE_VEHnormalized.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan;                         
                                % SEM
                                PooledSEM.(TrialDay{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan;  
                                PooledSEM.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan; 
                                % SEM, VEH normalized
                                PooledSEM_VEHnormalized.(TrialDay{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan;  
                                PooledSEM_VEHnormalized.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan;                                                               
                                % STREAMS                
                                fullstreams.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = {};
                                % STREAMS ALIGNED                
                                fullstreams_aligned.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = {};
                                % TIME VECT ALIGNED
                                timevect_aligned.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = {};
                                % TIME OPTO TRIALS
                                time_optotrials_aligned.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = {};
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
                            end
                        end
                    end
                end
                %STORE STREAMS
                fullstreams.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = streams;
                
                %STORE STREAMS ALIGNED TO OPTO ONSET #1
                fullstreams_aligned.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = streams_aligned;
                timevect_aligned.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = time_vect_aligned; % DETQ/VEH or DRUG levels
                time_optotrials_aligned.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = time_optotrials;

            end
        end
    end
end

%% Normalization against VEH; i.e. against the average of the maxima in the veh trials + calculate Zscore of this
% initialize
for v=1:length(TrialDay)
    for p=1:length(pooledtype)
        for d=1:length(dFF_names)
            for pow = 1:length(Opto_Powers)
                for k = 1:size(datatype,2)
                    for nummice=1:length(mice_list_virus.(TrialDay{v}))
                        for s = 1:length(sessions)
                            oiu=PooledINDIV.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s});
                            PooledINDIV_VEHnormalized.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = ...
                                ones(size(oiu,1),size(oiu,2))*nan;
                        end
                    end
                end
            end
        end
    end
end

% fill
for v=1:length(TrialDay)
    for p=1:length(pooledtype)
        for d=1:length(dFF_names)
            for pow = 1:length(Opto_Powers)
                for k = 1:size(datatype,2); % 
                    for nummice=1:length(mice_list_virus.(TrialDay{v}))
                        for s = 1:length(sessions)
                            PooledINDIV_VEHnormalized.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = ...
                                100*PooledINDIV.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s})./ ...
                                nanmean(max(PooledINDIV.VEH.(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}),[],2),1);
                        end
                    end
                end
            end
        end
    end
end

% store average
for v=1:length(TrialDay)
    for p=1:length(pooledtype)
        for d=1:length(dFF_names)
            for pow = 1:length(Opto_Powers)
                for k = 1:size(datatype,2)
                    for nummice=1:length(mice_list_virus.(TrialDay{v}))
                        for s = 1:length(sessions)
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

%% Simpler tracking of time for opto trials; granted only 1 session per animal
for nummice=1:length(mice_list_virus.(TrialDay{v}))
    for v=1:length(TrialDay)
        s=1;
        time_optotrials_aligned_simple.(TrialDay{v})(:,nummice) = time_optotrials_aligned.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s});  
    end
end

%% Save pooled data for all viruses together
% STORE INDIVIDUAL and AVERAGE TRIALS
save([PATH2SAVEPOOL,'\PooledAllMice.mat'],'PooledINDIV_VEHnormalized','PooledAVE_VEHnormalized','PooledSEM_VEHnormalized','PooledINDIV','PooledAVE','PooledSEM',...
    'TrialDay','PooledAnimalID','PooledAnimalID_only','dFF_names','datatype','Opto_Powers','t_trials','stim_duration','TRANGE','BASELINE_WIN','Events2Align2','mice_list_virus',...
    'dt_ds','sampling_rate_ds','length_data','fullstreams','time_vect','pooledtype','time_optotrials_aligned_simple','AnimalIDs','SessionIDs','SessionNames','fullstreams_aligned',...
    'time_vect_aligned');


