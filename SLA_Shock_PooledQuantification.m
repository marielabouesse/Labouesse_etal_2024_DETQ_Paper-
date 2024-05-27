% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT
% Name: Shock_PooledQuantification
% Cohort: GCAMP or dLight1.3b with shock

% Feb 2023 then March 2024 - Marie Labouesse, marie.labouesse@gmail.com

% 1- POOLED QUANTIFICATION
% loads FP data previously analyzed (with: Shock_DataExtraction then Shock_PooledGraphs) in matlab space "PooledAllMice.mat" from multiple animals (select folder containing multiple animals) 
% calculates AUC, Peak Maxima, DecayHalfTime (on individual traces or on the average)
% quantifications are done on individual trials as well as on the mean of all trials for 1 mouse (in 1 session)
% generates and saves the data into matlab spaces and tables that can be copy-pasted into GraphPad for stats

% edit relevant parameters in %% SETUP PARAMETERS

% functions needed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
for pow = 1:length(Opto_Powers) % 
    numtrials.(Opto_Powers{pow}) = size(PooledINDIV.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s}),1);
end

% WINDOW DURING WHICH TO CALCULATE AUC AND OTHER VARIABLES: 
windowcalc_duration = 10; %10 seconds


%% Compute the AUC and the Peaks
% Determine the indexes where to calculate the AUC and other variables in the stim epoch
temp_t = time_vect - time_vect(1); %values in timevect shifted by one value. We dont use the t_trials cos it doesnt start at 0
dummie = 1:length(temp_t); % indexes of time vector
dummie2 = dummie(temp_t >= windowcalc_duration); % so look for the indexes when the time vector is superior or equal to the stim duration
idx_AUC2 = dummie2(1); %find the first index when the time vector hit the exact value of the end of the stim duration; this is the second AUC index, or the number of indixes you need to go thru to have the stim duration
idx_AUC1 = find(t_trials == 0); % this is the index in t_trials where the stimulation starts
idx_AUC = idx_AUC1:(idx_AUC1+(idx_AUC2-1)); % we will be calculating the AUC from the beginning of the stimulation to the index just before the first index at the end of the stim

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
                            
                            
                            t_opto_i = PooledINDIV.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(AnimalIDs{nummice}).(SessionIDs{s})(1:numtrials.(Opto_Powers{pow}),:);
                            for w = 1:size(t_opto_i,1)   % = how many rows
                                %AUC, includes positive and negative values (no need to normalized cos fixed duration: 
                                dff_all = t_opto_i(w,:); % put the trace of this row in a new variable cos easier to handle
                                dff_stim = t_opto_i(w,idx_AUC);   %piece of dFF in the stim epoch we are interested in for this row
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
    
