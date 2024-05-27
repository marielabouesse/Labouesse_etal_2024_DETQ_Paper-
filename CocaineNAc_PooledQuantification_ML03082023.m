% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT
% Name: CocaineNAc_PooledQuantification
% Cohort: NAc AllodLightctrl with opto and drugs (DETQ or cocaine)

% Feb 2023 - Marie Labouesse, marie.labouesse@gmail.com

% 1- POOLED QUANTIFICATION
% loads FP data previously analyzed (with: CocaineNAc_DataExtraction then CocaineNAc_PooledGraphs) in matlab space "PooledAllMice.mat" from multiple animals (select folder containing multiple animals) 
% calculates AUC, Peak Maxima, DecayHalfTime (on individual traces or on the average), Baseline level (before the opto stim)
% quantifications are done on individual trials as well as on the mean of all trials for 1 mouse (in 1 session)
% generates and saves the data into matlab spaces and tables that can be copy-pasted into GraphPad for stats

% edit relevant parameters in %% SETUP PARAMETERS

% functions needed: ExponentialFit_DecayHalfTime


%% INITIALIZATIONS
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% SETUP PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD MATLAB SPACES
PATH2DATA = uigetdir('select folder'); %path 2 folder above folder of different groups
PATH2SAVEPOOL = PATH2DATA;
mkdir([PATH2SAVEPOOL,'\pooled figures\']);
load([PATH2DATA,'\PooledAllMice.mat']);

% SESSIONS ...............(only 1 session)
SessionIDs = {'SessionNum1'}; s=1;

% NUM TRIALS
v=1;d=1;k=1;nummice=1;
for pow = 1:length(Opto_Powers) % 0 and 2mW (for example)
    numtrials.(Opto_Powers{pow}) = size(PooledINDIV.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}),1);
end

% dFF TO USE FOR QUANTIFICATION
datatype_2use_for_graphs = 2; % 2 if you want to use dFFwithin

% WINDOW DURING WHICH TO CALCULATE AUC AND OTHER VARIABLES: 
windowcalc_duration = 60; %10 seconds; for the opto stim calculations.
prestim_epoch_start = -60;
prestim_epoch_stop = -1; 
prestim_duration = -(prestim_epoch_start - prestim_epoch_stop); % 9 seconds, ie -10 to -1 sec

% BASELINE: calculate between what trials
bas_t1=2;
bas_t2=3;

% CALCULATE OFF DECAY ? MANUALLY PER TRIAL TO LOOK AT THE FIT
offdecay_todo = 2; % 0 if no, 1 if on individual trials, 2 if on the average

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

%% Initialize the "Measurements" structure where you will be putting the AUC, Peaks and specific values. 
% Here I decided to save all the mice at the same place. That means I need to run this everytime again when I add mice. But also means this can be independent from the initial dataloading
% And accessing the data is easy; since now all the trials come from one matlab space for all groups.

for v=1:length(TrialDay) % which group
    for p = 1:length(pooledtype)   % raw or baselinecorr
        for d=1:length(dFF_names)   % if different recordings eg from 2 different brain regions
            for pow = 1:length(Opto_Powers) % different trial types; eg different opto powers or different drugs
                for k = 1:size(datatype,2) % ZScoredFF and dFF and dFF-within
                    for nummice=1:length(AnimalIDs)    % all mice
                        for s = 1:length(SessionIDs) % only 1 session for now
                            t_opto_i = PooledINDIV.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s})(1:numtrials.(Opto_Powers{pow}),:);
                            Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AUC = ones(size(t_opto_i,1),1)*nan; %same number of rows as trials
                            Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMaxima = ones(size(t_opto_i,1),1)*nan; %value of maxima
                            Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMinima = ones(size(t_opto_i,1),1)*nan; %value of minima
                            Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).DegChange = ones(size(t_opto_i,1),1)*nan; %degree of change (onset vs maxima)                            
                            Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).DecayHalfTime = ones(size(t_opto_i,1),1)*nan; %off decay
                            Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).Baseline = ones(size(t_opto_i,1),1)*nan; %baseline value (-10 to -1 sec before opto stim)
                        end
                    end
                end
            end
        end
    end
end

%% Determine the AUC for entire section (anything below 0 will count for negative AUC) ; and other variables
for v = 1:length(TrialDay) 
    for p = 1:length(pooledtype)   
        for d = 1:length(dFF_names)   
            for pow = 1:length(Opto_Powers) 
                for k = 1:size(datatype,2)
                    for nummice=1:length(AnimalIDs)
                        for s = 1:length(SessionIDs)
                            %%
                            % for off-decay; manual 
                            if offdecay_todo == 1 || offdecay_todo == 2
                                p=2; d=1; k=2; s=1;
                                v=2; 
                                pow = 1;
                                nummice = 5; 
                            end
                            
                            t_opto_i = PooledINDIV.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s})(1:numtrials.(Opto_Powers{pow}),:);
                            for w = 1:size(t_opto_i,1)   % = how many rows
                                %AUC, includes positive and negative values (no need to normalized cos fixed duration: 
                                dff_all = t_opto_i(w,:); % put the trace of this row in a new variable cos easier to handle
                                dff_stim = t_opto_i(w,idx_AUC);   %piece of dFF in the stim epoch we are interested in for this row
                                dff_prestim = t_opto_i(w,idx_AUC_pre); %piece of dFF in the pre-stim epoch we are interested in for this row
                                tmp_AUC = trapz(dff_stim);
                                Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AUC(w) = tmp_AUC; %same number of rows as trials

                                % Maxima in stim epoch  
                                AbsMaxima = max(dff_stim);
                                Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMaxima(w,:) = [AbsMaxima];                                
                                
                                % Minima in stim epoch   
                                AbsMinima = min(dff_stim);
                                Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMinima(w,:) = [AbsMinima];                                

                                % Degree change in stim epoch 
                                dFFonset = dff_stim(1);
%                                 dFFoffset = t_opto_i(w,idx_AUC(end));
                                DegChange = 100* (AbsMaxima - dFFonset)./dFFonset;
                                Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).DegChange(w,:) = [DegChange]; %value at start, end and degree of change
                                            
                                if offdecay_todo == 1
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
                                    expfit_type = 'exp1';
                                    StartPoint = [5,0];
                                    DecayHalfTime = ExponentialFit_DecayHalfTime(fit_timevect,object_2fit,expfit_type,StartPoint)
                                    Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).DecayHalfTime(w,:) = [DecayHalfTime]; %value at start, end and degree of change

                                elseif offdecay_todo == 0
                                end
                                    
                                % Baseline
                                BaselinedFF = nanmean(dff_prestim,2);
                                Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).Baseline(w,:) = [BaselinedFF]; %baseline value (-10 to -1 sec before opto stim)
                                
                            end
                            if offdecay_todo == 2
                                t_opto_i_mean = nanmean(t_opto_i,1);
                                dff_stim = t_opto_i_mean(1,idx_AUC);
                                AbsMaxima_mean = max(dff_stim);
                                AbsMaximaIndex_mean = find(t_opto_i_mean == AbsMaxima_mean);
                                t_opto_i_postpeak_mean = t_opto_i_mean(1,AbsMaximaIndex_mean:AbsMaximaIndex_mean+round(20./dt_ds));    %AbsMaximaIndex_mean+round(120./dt_ds) AbsMaximaIndex_mean+round(220./dt_ds)
                                % fit it with an exponential fit and determine off-decay using the function ExponentialFit_DecayHalfTime
                                % Need to define input parameters: time vector, fitobject, expfit_type: exp1, and start point [x,x]
                                fit_timevect_mean = 0:dt_ds:dt_ds*(length(t_opto_i_postpeak_mean)-1);
                                fit_timevect_mean = fit_timevect_mean';
                                object_2fit_mean = smooth(t_opto_i_postpeak_mean,100);
    %                                 object_2fit = t_opto_i_postpeak';
                                expfit_type = 'exp2';
                                StartPoint = [5,0];
                                DecayHalfTime = ExponentialFit_DecayHalfTime(fit_timevect_mean,object_2fit_mean,expfit_type,StartPoint)              
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay{v}).DecayHalfTime(nummice,pow) = DecayHalfTime;
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
for v = 1:length(TrialDay) 
    for p = 1:length(pooledtype)   
        for d = 1:length(dFF_names)   
            for pow = 1:length(Opto_Powers) 
                for k = 1:size(datatype,2) 
                    for nummice=1:length(AnimalIDs)
                        for s = 1:length(SessionIDs)
                            if s == 1
                                Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay{v}).AUC(nummice,pow) = ...
                                    nanmean(Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AUC,1)
                                 Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay{v}).AbsMaxima(nummice,pow) = ...
                                    nanmean(Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMaxima,1)
                                 Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay{v}).AbsMinima(nummice,pow) = ...
                                    nanmean(Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).AbsMinima,1)
                                 Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay{v}).DegChange(nummice,pow) = ...
                                    nanmean(Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).DegChange,1)
                                 if offdecay_todo == 1
                                 Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay{v}).DecayHalfTime(nummice,pow) = ...
                                    nanmean(Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).DecayHalfTime,1)
                                 else
                                 end
                                 Measurements_AVE.(pooledtype{p}).(dFF_names{d}).(datatype{k}).(TrialDay{v}).Baseline(nummice,pow) = ...
                                    nanmean(Measurements.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}).Baseline(bas_t1:bas_t2,:),1);                            
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
save([PATH2SAVEPOOL,'Cohort_analysis.mat'],'Measurements','Measurements_AVE','TrialDay','PooledAnimalID','dFF_names','datatype','Opto_Powers',...
    't_trials','stim_duration','TRANGE','BASELINE_WIN','Events2Align2','dt_ds','sampling_rate_ds','length_data','fullstreams','time_vect','pooledtype');
    








