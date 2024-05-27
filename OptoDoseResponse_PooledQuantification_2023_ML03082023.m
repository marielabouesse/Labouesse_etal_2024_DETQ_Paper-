% OptoDoseResponse_PooledQuantification_2023

% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT
% Name: OptoDoseResponse_PooledQuantification_2023
% Cohort: PFC dLight1.3b with opto 

% May 2023 - Marie Labouesse, marie.labouesse@gmail.com

% 1- POOLED QUANTIFICATION
% loads FP data previously analyzed (with: OptoDoseResponse_DataExtraction then OptoDoseResponse_PooledMatrices) in matlab space "PooledAllMice.mat" from multiple animals (select folder containing multiple animals) 
% or if used as a function the data was loaded previously within OptoDoseResponse_PooledOverall_2023 or OptoDoseResponse_TempWindow_2023
% calculates AUC, Peak Maxima, DecayHalfTime (on individual traces or on the average), Baseline level (before the opto stim)
% quantifications are done on individual trials as well as on the mean of all trials for 1 mouse (in 1 session)
% generates and saves the data into matlab spaces and tables that can be copy-pasted into GraphPad for stats

% edit relevant parameters in %% SETUP PARAMETERS or within the code running above this one

% functions needed: ExponentialFit_DecayHalfTime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('runasindependentscript') % then we don't need the below as the other scripts already initialized everything we need
    %% INITIALIZATIONS  --> run this if you want to use this code independently.
    close all; clear variables; clc;
    set(0,'defaultfigurecolor',[1 1 1])

    % SETUP PARAMETERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD MATLAB SPACES
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
    datatype_2use_for_graphs = 2; % 2 if you want to use dFFwithin within 'datatypes'

    % data4measurements
    data4measurements = 0; % 0 if regular, 1 if normalized;

    % WINDOW DURING WHICH TO CALCULATE AUC AND OTHER VARIABLES: 
    windowcalc_duration = 60; %60 seconds
    prestim_epoch_start = -60;
    prestim_epoch_stop = -1; 
    prestim_duration = -(prestim_epoch_start - prestim_epoch_stop); % 9 seconds, ie -10 to -1 sec

    % CALCULATE OFF DECAY ? MANUALLY PER TRIAL TO LOOK AT THE FIT
    offdecay_todo = 1; % 0 if no, 1 if on individual trials, 2 if on the average

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
    mkdir([PATH2SAVEPOOL_SELECT]);

    % missing
    for v=1:length(TrialDay_2analyze)
        mice_list_virus.(TrialDay_2analyze{v}) = mice_list_virus.(TrialDay{v});
    end
    % pooledtype
    pooledtype={'raw','baselinecorr'}
    animals = AnimalIDs;
end
% ABOVE THIS LINE NOT IN USE WHEN USED AS A INDEPENDENT SCRIPT CALLED OUT FROM OTHER SCRIPTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% info about what will be cropped; we only want 5 trials for VEH and 14 for DETQ (or less if we are focusing on defining the temporal window
if length(TrialDay_2analyze) == 1
    if Trials_2analyze.DETQ(1) == 1; trialcrop_at_start_DETQ = 0; else trialcrop_at_start_DETQ = 1; end % if we are chopping at the start of DETQ: default is start at 1
    if Trials_2analyze.DETQ(end) == 18; trialcrop_at_end_DETQ = 0; else trialcrop_at_end_DETQ = 1; end % if we are chopping at the end of DETQ: default is end at 18
else
    if Trials_2analyze.VEH(1) == 1; trialcrop_at_start_VEH = 0; else trialcrop_at_start_VEH = 1; end % if we are chopping at the start of VEH: default is start at 1
    if Trials_2analyze.VEH(end) == 6; trialcrop_at_end_VEH = 0; else trialcrop_at_end_VEH = 1; end % if we are chopping at the end of VEH: default is end at 6

    if Trials_2analyze.DETQ(1) == 1; trialcrop_at_start_DETQ = 0; else trialcrop_at_start_DETQ = 1; end % if we are chopping at the start of DETQ: default is start at 1
    if Trials_2analyze.DETQ(end) == 18; trialcrop_at_end_DETQ = 0; else trialcrop_at_end_DETQ = 1; end % if we are chopping at the end of DETQ: default is end at 18
end

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

%  % NUM TRIALS
for v=1:length(TrialDay_2analyze) % which group
    numtrials_final.(TrialDay_2analyze{v}) = Trials_2analyze.(TrialDay_2analyze{v})(end) - Trials_2analyze.(TrialDay_2analyze{v})(1) +1;
end

%% Edit PooledINDIV
%% CROP THE DATA TRIALS NOT TO ANALYZE
if alreadycroppeddata == 0
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
                                 PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = ...
                                     PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s})(trials_2_keep,:);
                                 % PooledINDIV_VEHnormalized
                                 PooledINDIV_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}) = ...
                                     PooledINDIV_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s})(trials_2_keep,:);
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
    % NUM TRIALS
    for v=1:length(TrialDay_2analyze) % which group
        numtrials_final.(TrialDay_2analyze{v}) = Trials_2analyze.(TrialDay_2analyze{v})(end) - Trials_2analyze.(TrialDay_2analyze{v})(1) +1;
    end
end

%% Initialize the "Measurements" structure where you will be putting the AUC, Peaks and specific values. 
% Here we decide to save all the mice at the same place. That means we need to run this everytime again when we add mice. But also means this can be independent from the initial dataloading
% And accessing the data is easy; since now all the trials come from one matlab space for all groups.

for v=1:length(TrialDay_2analyze) % which group
    for p = 1:length(pooledtype)   % raw or baselinecorr
        for d=1:length(dFF_names)   % if different recordings eg from 2 different brain regions
            for pow = 1:length(Opto_Powers) % different trial types; eg different opto powers or different drugs
                for k = 1:size(datatype,2) % ZScoredFF and dFF and dFF-within
                    for nummice=1:length(AnimalIDs)    % all mice
                        for s = 1:length(SessionIDs) % only 1 session for now
                            if data4measurements == 0% regular
                                t_opto_i = PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s})(1:numtrials_final.(TrialDay_2analyze{v}),:);
                            elseif data4measurements == 1% VEHnormalized                            
                                t_opto_i = PooledINDIV_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s})(1:numtrials_final.(TrialDay_2analyze{v}),:);
                            end
                            Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AUC = ones(size(t_opto_i,1),1)*nan; %same number of rows as trials
                            Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMaxima = ones(size(t_opto_i,1),1)*nan; %value of maxima
                            Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMinima = ones(size(t_opto_i,1),1)*nan; %value of minima
                            Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).DegChange = ones(size(t_opto_i,1),1)*nan; %degree of change (onset vs maxima)                            
                            Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).DecayHalfTime = ones(size(t_opto_i,1),1)*nan; %off decay
                            Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).Baseline = ones(size(t_opto_i,1),1)*nan; %baseline value (-10 to -1 sec before opto stim)
                        end
                    end
                end
            end
        end
    end
end

%% Determine the AUC for entire section (anything below 0 will count for negative AUC) ; and other variables
for v = 1:length(TrialDay_2analyze) 
    if length(TrialDay_2analyze) == 1 % only DETQ
        v=2;
    end    
    for p = 1:length(pooledtype)   
        for d = 1:length(dFF_names)   
            for pow = 1:length(Opto_Powers) 
                for k = 1:size(datatype,2)
                    for nummice=1:length(AnimalIDs)
                        for s = 1:length(SessionIDs)
                            %%
                            % for off-decay; manual running
                            if offdecay_todo == 1 || offdecay_todo == 2 % 1 if on individual trials, 2 if average
                                p=2; d=1; k=datatype_2use_for_graphs; s=1; pow = 1;
                                v=1;
                                nummice = 1; 
                            end
                            
                            % which data
                            if data4measurements == 0% regular
                                t_opto_i = PooledINDIV.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s})(1:numtrials_final.(TrialDay_2analyze{v}),:);
                            elseif data4measurements == 1% VEHnormalized                            
                                t_opto_i = PooledINDIV_VEHnormalized.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s})(1:numtrials_final.(TrialDay_2analyze{v}),:);
                            end
                            
                            for w = 1:size(t_opto_i,1)   % = how many rows
                                %AUC, includes positive and negative values (no need to normalized cos fixed duration: 
                                dff_all = t_opto_i(w,:); % put the trace of this row in a new variable cos easier to handle
                                dff_stim = t_opto_i(w,idx_AUC);   %piece of dFF in the stim epoch we are interested in for this row
                                dff_prestim = t_opto_i(w,idx_AUC_pre); %piece of dFF in the pre-stim epoch we are interested in for this row
                                tmp_AUC = trapz(dff_stim);
                                Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AUC(w) = tmp_AUC; %same number of rows as trials

                                % Peak Maxima in stim epoch  
                                AbsMaxima = max(findpeaks(dff_stim));                                
%                                 AbsMaxima = max(dff_stim);
                                Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMaxima(w,:) = [AbsMaxima];                                
                                
                                % Minima in stim epoch   
                                AbsMinima = min(dff_stim);
                                Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMinima(w,:) = [AbsMinima];                                

                                % Degree change in stim epoch 
                                dFFonset = dff_stim(1);
%                                 dFFoffset = t_opto_i(w,idx_AUC(end));
                                DegChange = 100* (AbsMaxima - dFFonset)./dFFonset;
                                Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).DegChange(w,:) = [DegChange]; %value at start, end and degree of change
                                            
                                if offdecay_todo == 1 % on individual peaks
                                    % Off decay
                                    % select the data as of the peak maxima
                                    AbsMaximaIndex = find(t_opto_i(w,:) == AbsMaxima);
                                    t_opto_i_postpeak = t_opto_i(w,AbsMaximaIndex:AbsMaximaIndex+round(10./dt_ds));    %AbsMaximaIndex+round(120./dt_ds)
                                    % fit it with an exponential fit and determine off-decay using the function ExponentialFit_DecayHalfTime
                                    % Need to define input parameters: time vector, fitobject, expfit_type: exp1, and start point [x,x]
                                    fit_timevect = 0:dt_ds:dt_ds*(length(t_opto_i_postpeak)-1);
                                    fit_timevect = fit_timevect';
                                    object_2fit = smooth(t_opto_i_postpeak,100);
    %                                 object_2fit = t_opto_i_postpeak';
                                    expfit_type = 'exp2'; %exp1 or 2
                                    StartPoint = [5,0];
                                    DecayHalfTime = ExponentialFit_DecayHalfTime(fit_timevect,object_2fit,expfit_type,StartPoint)
                                    Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).DecayHalfTime(w,:) = [DecayHalfTime]; %value at start, end and degree of change

                                elseif offdecay_todo == 0
                                end
                                    
                                % Baseline
                                BaselinedFF = nanmean(dff_prestim,2);
                                Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).Baseline(w,:) = [BaselinedFF]; %baseline value (-10 to -1 sec before opto stim)
                                
                            end
                            if offdecay_todo == 2 % on the average
                                t_opto_i_mean = nanmean(t_opto_i,1);
                                dff_stim = t_opto_i_mean(1,idx_AUC);
                                AbsMaxima_mean = max(dff_stim);
                                AbsMaximaIndex_mean = find(t_opto_i_mean == AbsMaxima_mean);
                                t_opto_i_postpeak_mean = t_opto_i_mean(1,AbsMaximaIndex_mean:AbsMaximaIndex_mean+round(30./dt_ds));    %AbsMaximaIndex_mean+round(60./dt_ds)
                                % fit it with an exponential fit and determine off-decay using the function ExponentialFit_DecayHalfTime
                                % Need to define input parameters: time vector, fitobject, expfit_type: exp1, and start point [x,x]
                                fit_timevect_mean = 0:dt_ds:dt_ds*(length(t_opto_i_postpeak_mean)-1);
                                fit_timevect_mean = fit_timevect_mean';
                                object_2fit_mean = smooth(t_opto_i_postpeak_mean,100);
    %                                 object_2fit = t_opto_i_postpeak';
                                expfit_type = 'exp2'; %exp1 or 2
                                StartPoint = [5,0];
                                DecayHalfTime = ExponentialFit_DecayHalfTime(fit_timevect_mean,object_2fit_mean,expfit_type,StartPoint)              
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay_2analyze{v}).DecayHalfTime(nummice,pow) = DecayHalfTime;
                            end
                            
                            %%
                        end
                    end
                end
            end
        end
    end
end


%% Generate the AVE Measurements data (average for each mouse, for all its trials) + easy "copy-paste" table for plotting into graph pad (mice are in rows, groups in columns)
for v = 1:length(TrialDay_2analyze) 
    if length(TrialDay_2analyze) == 1 % only DETQ
        v=2;
    end
    for p = 1:length(pooledtype)   
        for d = 1:length(dFF_names)   
            for pow = 1:length(Opto_Powers) 
                for k = 1:size(datatype,2) 
                    for nummice=1:length(AnimalIDs)
                        for s = 1:length(SessionIDs)
                            if s == 1
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay_2analyze{v}).AUC(nummice,pow) = ...
                                    nanmean(Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AUC,1)
                                 Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay_2analyze{v}).AbsMaxima(nummice,pow) = ...
                                    nanmean(Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMaxima,1)
                                 Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay_2analyze{v}).AbsMinima(nummice,pow) = ...
                                    nanmean(Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMinima,1)
                                 Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay_2analyze{v}).DegChange(nummice,pow) = ...
                                    nanmean(Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).DegChange,1)
                                 if offdecay_todo == 1
                                 Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay_2analyze{v}).DecayHalfTime(nummice,pow) = ...
                                    nanmean(Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).DecayHalfTime,1)
                                 else
                                 end
                                 Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay_2analyze{v}).Baseline(nummice,pow) = ...
                                    nanmean(Measurements.(TrialDay_2analyze{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).Baseline,1)                            
                            elseif length(sessions) > 1
                                % N/A
                            end
                        end
                    end
                end
            end
        end
    end
end



%% Save data
save([PATH2SAVEPOOL_SELECT,'\Cohort_analysis.mat'],'Measurements','Measurements_AVE','TrialDay','PooledAnimalID','dFF_names','datatype','Opto_Powers',...
    't_trials','stim_duration','TRANGE','BASELINE_WIN','Events2Align2','dt_ds','sampling_rate_ds','length_data','fullstreams','time_vect','pooledtype');
    
