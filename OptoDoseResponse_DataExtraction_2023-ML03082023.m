% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT
% Name: OptoDoseResponse_DataExtraction
% Cohort: PFC dLight1.3b with opto-dose response and temporal window

% May 2023 - Marie Labouesse, marie.labouesse@gmail.com

% extracts FP data (select folder from a particular day eg day 1, containing multiple animals inside)
% generates and saves a matrix with individual trials aligned to specific event (opto stim)
% generates and saves graphs aligned to specific event
% keep tracks of timestamp for each opto stim relative to the drug injection
% saves trials separately per drug (veh, drug)

% edit relevant parameters in %% SETUP PARAMETERS

% functions needed: error_area_onlyrectangle, TDTbin2mat, SEV2mat

%% INITIALIZATIONS
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% SETUP PARAMETERS
%%%%%%%%%% TTL EXTRACTION
ttl_extract = 'default'; %'default' or 'special' (needs to be 7 characters)
            
%%%%%%%%%% DATA EXTRACTION
channel_number = 1; % 1 if 1-site recording, 2 if 2-site recording
channel_names = {'PFC'}; % eg brain region
color_number = 2; % 2 if 405 and 465 recording, 3 if also 565
color_names = {'c405','c465'}; % can also be 565
dFF_number = 1; % how many dFF to analyze, can be 1 or 2, based on number of channels and colors (if more than 2, need to edit script)
dFF_names = {'PFC'}; %put them in this order: 1) channel 1, 465 and 2) channel 1, 560 OR channel 2: 465 (if other combinations, need to edit script)
datatype = {'dFF','dFFwithin','ZScoredFF_within','ZScoredFF'}; %data analyses of interest. Note the ZscoredFF_within here is just equal to ZscoredFF; as the initial calc were not conclusive.
keeprawdata = 1; %1 if you want to keep the raw data to calculate dFF within trials. 0 otherwise ; % make sure to add dFFwithin to datatype
raw_or_corr = {'raw','baselinecorr'};

%%%%%%%%%%% PREPROCESSING
% trimming
timetrim_start_set = 180; %time to chop off at start of recording, in seconds: adjust this manually, eg 60 seconds
timetrim_end_set = 1; %time to chop off at end of recording, eg 1 second
% low pass filtering of raw 405 and 470
lowpassfiltering = 1; % 1 if you want to low pass filter- typically yes- 0 for no
lowpassfreq = 1; % typically 1Hz
% detrending dFF
detrending = 0;% %% Write 1 for detrend dFF, normal detrend function of Matlab. Typically no for open field.
% substract 8th percentile as a baseline calculate using a moving window
deletemovingwindow = 0; % 1 to delete moving baseline 8th percentile, 0 if no. Typically no for open field.
% high pass filtering of dFF
highpassfiltering = 0; %'yes' if you want to high pass filter, 0 for no. Typically no for open field.
highpassfreq = 0.005; % typically 0.005Hz

%%%%%%%%%%% TRIAL DEFINITION
% time window for the trials and graphs
TRANGE = [-60 60]; % will create events for a -60 to +60 sec window (you can always plot less later)
% baseline correction window
BASELINE_WIN.OptoStim = [-10 -1];% baseline correction window. Typically -10 to -1 
% variables to align to
Events2Align2 = {'PulseStart_time','PulseStop_time'}; %can be {'optostim_onset','optostim_offset'} for opto, or any other TTL epoc in the "epocs" structure with same organization

%%%%%%%%%%% EXPERIMENT TYPE
stim_interval = 5; % 5 minutes
stim_duration = 2; % 2 seconds
stim_duration2 = [];  % in case their is a ramp or another duration of interest
extratime_prefirstopto = 60; % 60 sec; how much of the data to take for the graphs with the traces over time: how much time to take before the start of first opto stim
extratime_postlastopto = 240; % 240 sec; how much of the data to take for the graphs with the traces over time: how much time to take after end of last opto stim

%%%%%%%%%%% HOW MANY MICE TO ANALYZE
LoopOrNot = 1; % 1 to loop, 0 if you only want to test one mouse and dont want to loop --> if so, edit the mouse number you want below: nummice_set
nummice_set = 1;

%%%%%%%%%%% SHOW, SAVE, OVERWRITE PARAMETERS
show_plot = 1; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
reanalysis = 1; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed. 
overwrite = 1; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not          

%% DEFINE PATHS
path2data = uigetdir('select folder'); % SELECT FOLDER OF GROUP TO ANALYZE - %Location of the data
mice_list = dir(path2data); %all things in this folder
% Define paths and data to analyze
path2savefolder = path2data; %Path to save

%% IDENTIFY MICE TO ANALYZE 
for o = length(mice_list):-1:1
    if mice_list(o).isdir == 0  %remove non-folders
        mice_list(o) = [];
    else
        if  strcmp(mice_list(o).name,'data') == 1 || strcmp(mice_list(o).name,'figures') == 1 ...   %remove folders with data or figures
            || contains(mice_list(o).name,'data') || contains(mice_list(o).name,'figures') || contains(mice_list(o).name,'representative figures')...
            || contains(mice_list(o).name,'results') || contains(mice_list(o).name,'other') ...
            || strcmp(mice_list(o).name,'.') == 1 || strcmp(mice_list(o).name,'..') == 1
            mice_list(o) = [];
        end
    end
end
Nmice = length(mice_list);  %number of mice to analyze

%% LOOP ACROSS MICE
if LoopOrNot == 0  %if you only want to analyze one animal
    Nmice = 1;
end
%%
for nummice = 1:Nmice   
    if LoopOrNot == 0
        nummice = nummice_set;  %need to set the number of the mouse you want to analyze up above in PARAMETERS
    end
    AnimalID = ['ID_',mice_list(nummice).name(end-3:end)];   %last 5 digits of folders should be the animal name
    Dirsplit = strsplit(path2data,'\');
    Virus_cell = Dirsplit(length(Dirsplit));
    Virus = Virus_cell{:}; %just extracting the string out of the cell array   --> (name of the folder you selected)
    BehavData = {};

    %% Setting the paths for the mouse and the sessions
    % Define the path to the mouse and find the folders to analyze:
    path2mouse = [path2data,'\',mice_list(nummice).name,'\'];
    % if there are several sessions inside the mouse folder
    sessions = dir(path2mouse);
    %remove non relevant folders
    for o = length(sessions):-1:1
        if sessions(o).isdir == 0  %remove non-folders
            sessions(o) = [];
        elseif strcmp(sessions(o).name,'data') == 1 || strcmp(sessions(o).name,'figures') == 1 || strcmp(sessions(o).name,'representative figures') == 1 || strcmp(sessions(o).name,'results') == 1 ...
             || strcmp(sessions(o).name,'.') == 1 || strcmp(sessions(o).name,'..') == 1 || contains(sessions(o).name,'.csv') == 1
             sessions(o) = [];
        end
    end
    if isempty(sessions)
        sessions(1).name = []; %to have an existent session
    end
    
    % Create a folder to save the data for this mouse and define the path
    if exist([path2savefolder,'\',mice_list(nummice).name,'\'],'dir') == 0
        mkdir([path2savefolder,'\',mice_list(nummice).name,'\'])
    end
    path2save_mouse = [path2savefolder,'\',mice_list(nummice).name,'\'];       
    
    %% Loop for all the sessions for the mouse
    for s = 1; %:length(sessions)
       
        % Define the path to the session and create folder to save if needed:
        PATH2SESSION = [path2mouse,sessions(s).name];
        if exist([path2save_mouse,sessions(s).name],'dir') == 0
            mkdir([path2save_mouse,sessions(s).name])
        end
        if  length(sessions) == 1
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
        
        % PATH2SAVEWORKSPACE
        Dirsplit = strsplit(path2data,'\');
        Drug_cell = Dirsplit(length(Dirsplit));
        Drug_notcell = Drug_cell{:}; %just extracting the string out of the cell array   --> (name of the folder you selected)
        if contains(path2data,'2022')
            Opto_Powers = {'VEH','DETQ'}; drug_name = Opto_Powers;
            %%%%%%%%%%%% TRIALS BY DRUG
            VEH_trial1 = 1; VEH_trialend = 6; DETQ_trial1 = 7; DETQ_trialend = 24; % how the different trials were organized relative to drug injection
            PATH2SAVEWORKSPACE_DETQ = [fileparts(fileparts(path2data)),'\_ALLDATA\',Drug_notcell,'\','DETQ','\',AnimalID,'\']
            if exist([fileparts(fileparts(path2data)),'\_ALLDATA\',Drug_notcell,'\','DETQ','\',AnimalID,'\'],'dir') == 0
                mkdir([fileparts(fileparts(path2data)),'\_ALLDATA\',Drug_notcell,'\','DETQ','\',AnimalID,'\'])
            end
            PATH2SAVEWORKSPACE_VEH = [fileparts(fileparts(path2data)),'\_ALLDATA\',Drug_notcell,'\','VEH','\',AnimalID,'\'];
            if exist([fileparts(fileparts(path2data)),'\_ALLDATA\',Drug_notcell,'\','VEH','\',AnimalID,'\'],'dir') == 0;
                mkdir([fileparts(fileparts(path2data)),'\_ALLDATA\',Drug_notcell,'\','VEH','\',AnimalID,'\']);
            end
        end
        
        % Check if results are already saved for this session 
        done = 0; % not in use %exist([PATH2SAVE,'IndividualData',AnimalID,'.mat'],'file'); %if you specify 'file' then matlab searches for both files and folders
        %%
        if done == 0 || reanalysis == 1 %if not done or if wish to reanalyze


            %% Select FP data file and extract FP data streams 
            files_tev = dir(fullfile(PATH2SESSION,'*.tev'));
            trace_name = files_tev.name(1:end-4); %remove .mat or.tev
            data = TDTbin2mat(PATH2SESSION, 'TYPE', {'epocs', 'scalars', 'streams'});
            
            %Create 'streams' structure, which will contain extracted raw data from each channel/color
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    streams.rawdata.(channel_names{channel}).(color_names{colore}) = [];
                end
            end
            % first channel
            stream_name_type405A = {'x05A','05A','A05A','x405A'}; %possible names on TDT rigs- can add to this list
            stream_name_type465A = {'x65A','65A','A65A','x465A'}; %stick to same order on these 3 rows (names that go together on same column)
            stream_name_type560A = {'x60A','60A','A60A','x560A'};
            for i=1:length(stream_name_type405A)
                if isfield(data.streams, stream_name_type405A{i}); %isfield checks if 'c' is a field inside structure a.b.c: isfield(a,'c')
    %                     fieldnames(data.streams) %to get directly the field names in data.streams
                    streams.rawdata.(channel_names{1}).c405 = data.streams.(stream_name_type405A{i}).data; %access data in data.stream using the name of the 405 store
                    streams.rawdata.(channel_names{1}).c465 = data.streams.(stream_name_type465A{i}).data; %access data in data.stream using the name of the 465 store
                    % if also recording in red
                    if color_number == 3           
                        streams.rawdata.(channel_names{1}).c560 = data.streams.(stream_name_type560A{i}).data; %access data in data.stream using the name of the 560 store
                    end
                end
            end
            % second channel (if there is)
            if channel_number == 2
                stream_name_type405B = {'x05B','05B','B05B','x405B'};
                stream_name_type465B = {'x65B','65B','B65B','x405B'};
                stream_name_type560B = {'x60A','60A','A60A','x405B'};
                for i=1:length(stream_name_type405B)
                    if isfield(data.streams, stream_name_type405B{i}) %isfield checks if 'c' is a field inside structure a.b.c: isfield(a,'c')
                        streams.rawdata.(channel_names{2}).c405 = data.streams.(stream_name_type405B{i}).data; %access data in data.stream using the name of the 405 store
                        streams.rawdata.(channel_names{2}).c465 = data.streams.(stream_name_type465B{i}).data; %access data in data.stream using the name of the 465 store
                        % if also recording in red
                        if color_number == 3           
                            streams.rawdata.(channel_names{2}).c560 = data.streams.(stream_name_type560B{i}).data; %access data in data.stream using the name of the 560 store
                        end
                    end
                end
            end
            %Calculate min length
            length_stream = [];
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    length_stream(end+1) = length(streams.rawdata.(channel_names{channel}).(color_names{colore}));
                end
            end
            min_length_stream = min(length_stream);
            %Adjust the lengths of the channels
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    if length(streams.rawdata.(channel_names{channel}).(color_names{colore})) > min_length_stream
                       streams.rawdata.(channel_names{channel}).(color_names{colore})...
                        (min_length_stream+1:length(streams.rawdata.(channel_names{channel}).(color_names{colore}))) = []; 
                    end
                end
            end        

            %% Downsampling FP data, create a new structure for it (streams) and define time vector
            % sampling rate
            for h=1:length(stream_name_type405A)
                if isfield(data.streams, stream_name_type405A{h}) 
                    sampling_rate = data.streams.(stream_name_type405A{1}).fs;
                end
            end
            N = 10; %downsample 10 times (from initial sampling rate of 1017.3 Hz)
            sampling_rate_ds = sampling_rate/N;
            % downsample
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    streams.rawdata.(channel_names{channel}).(color_names{colore}) = downsample(streams.rawdata.(channel_names{channel}).(color_names{colore}),10);
                end
            end
            length_data = length(streams.rawdata.(channel_names{channel}).(color_names{colore})); %will just take it from the last data
            % time vector
            max_idx = length_data;
            dt_ds = 1/sampling_rate_ds; % delta t or unit of time (in seconds) per sample --> for 101.73 Hz sampling rate (after downsampling), will be close to 10 ms (9.8 ms)
            time_vect = 0:dt_ds*1:dt_ds*(max_idx-1); %time vector; values in seconds, same number of data points as the data


            %% Extract TTL data (epocs) for opto stim and camera frames
            % Initialize epocs structure, for all possible TTLs we may want to have, depending on the experiment. This is generic for experiments done in various labs, not specific to this particular paper/project.
            if ttl_extract == 'default' 
                epoc_type = {'Drug_time','PulseStart_time','PulseStop_time','StimAmp','StimTime','TTLstart_time','TTLstop_time','Testsync_start_time','Testsync_stop_time','Framesync_start_time',...
                'Possync_start_time','FrontCam_start_time','BackCam_start_time'};
                epocs = struct();
                for i=1:length(epoc_type)
                    epocs.(epoc_type{i}) = [];
                end
                % Overall opto pulse
                if isfield(data.epocs,'Pu1_');
                    epocs.PulseStart_time = data.epocs.Pu1_.onset;
                    epocs.PulseStop_time = data.epocs.Pu1_.offset;
                end
                % Overall opto pulse
                if isfield(data.epocs,'U22_');
                    epocs.Drug_time = data.epocs.U22_.onset;
%                     epocs.Drug_time = data.epocs.U22_.offset;
                end
                %Opto Stim amplitude
                if isfield(data.epocs,'Am1_');
                    epocs.StimAmp = data.epocs.Am1_.data;
                    epocs.StimTime = data.epocs.Am1_.onset;
                end
                % Generic TTL in PCO
                if isfield(data.epocs,'PC0_');
                    epocs.TTLstart_time = data.epocs.PC0_.onset;   
                    epocs.TTLstop_time = data.epocs.PC0_.offset;   
                end
                % Test start/stop sync TTLs arriving from Anymaze
                if isfield(data.epocs,'PC1_');
                    epocs.Testsync_start_time = data.epocs.PC1_.onset;
                    epocs.Testsync_stop_time = data.epocs.PC1_.offset;
                end
                % Cam Frame sync TTLs arriving from Anymaze
                if isfield(data.epocs,'PC2_');
                    epocs.Framesync_start_time = data.epocs.PC2_.onset;    %PC2 in new program
                end 
                % % Position sync TTLs arriving from Anymaze
                if isfield(data.epocs,'PC3_');
                    epocs.Possync_start_time = data.epocs.PC3_.onset;
                end
                % % Cam1 camera frames taken directly in Synapse
                if isfield(data.epocs,'Cam1');
                    epocs.FrontCam_start_time = data.epocs.Cam1.onset;
                end
                % % Cam2 camera frames taken directly in Synapse
                if isfield(data.epocs,'Cam2');
                    epocs.BackCam_start_time = data.epocs.Cam2.onset;
                end
                % if the TTLs are organized differently:
            elseif ttl_extract == 'special'
                epoc_type = epoc_names;
                for i=1:length(epoc_type)
                    epocs.(epoc_type{i}) = data.epocs.(epoc_locations{i}).(epoc_locations2{i});
                end
            end     

            %% TRIMMING
            % Setting up trim indexes
            timetrim_start = timetrim_start_set*1; %in seconds: adjust this manually, eg trim 3 minutes: 60*3 --> set up at the beginning of the code
            
             
            timetrim_end =  timetrim_end_set*1; %in seconds: adjust this manually --> set up at the beginning of the code
            dummie1 = 1:length(time_vect); %indexes, same number as there are samples in time_vect
            dummie2 = dummie1(time_vect > timetrim_start); %only keep the indexes starting from the new start
            idx_start = dummie2(1); %index for the new start
            dummie2 = dummie1(time_vect > time_vect(end) - timetrim_end); %only keep the indexes from the new end to the current end
            idx_end = dummie2(1); %index for the new end
            clear dummie1 dummie2 
                
            % Trimming and readjusting everything
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    streams.rawdata_trim.(channel_names{channel}).(color_names{colore}) = streams.rawdata.(channel_names{channel}).(color_names{colore})(idx_start:idx_end);
                end
            end        
            time_vect_trim = time_vect(idx_start:idx_end)-time_vect(idx_start);
            %pulses
            epocs_trim = struct();
            for i=1:length(epoc_type)
                epocs_trim.(epoc_type{i}) = epocs.(epoc_type{i}) - timetrim_start;
            end
            
            % Trimming: Plot
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    figure; clf; 
                    if show_plot == 0
                       set(gcf,'visible','off')
                    end
                    %Plot after
                    subplot(2,1,1);
                    plot(time_vect,streams.rawdata.(channel_names{channel}).(color_names{colore})); hold on; 
                    for o = 1:size(epocs.PulseStart_time)
                        plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-4 12],'m')
                        plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-4 12],'c')
                    end
                    min_x=-200; max_x=max(time_vect)+200; min_y=mean(streams.rawdata.(channel_names{channel}).(color_names{colore}))-100; 
                    max_y=mean(streams.rawdata.(channel_names{channel}).(color_names{colore}))+100; axis([min_x max_x min_y max_y]); 
                    xlabel('time(sec)'); ylabel('fluorescence'); title('Data before trimming'); 
                    xline(timetrim_start,':k','trim'); xline(time_vect(end)-timetrim_end,':k','trim'); 
                    %Plot after
                    subplot(2,1,2); plot(time_vect_trim,streams.rawdata_trim.(channel_names{channel}).(color_names{colore})); hold on; 
                    for o = 1:size(epocs_trim.PulseStart_time)
                        plot([epocs_trim.PulseStart_time(o) epocs_trim.PulseStart_time(o)],[-4 12],'m')
                        plot([epocs_trim.PulseStop_time(o) epocs_trim.PulseStop_time(o)],[-4 12],'c')
                    end
                    axis([min_x max_x min_y max_y]); %ie same as other subplot
                    xlabel('time(sec)'); ylabel('fluorescence'); title('Data after trimming'); 
                    annotation('textbox',[.25 0 .1 .17],'string','Ru ok with the trimming? if not readjust','FitBoxToText','on','EdgeColor','none');
                    sgtitle([channel_names{channel},' ',color_names{colore}]);
                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\trimming ',channel_names{channel},' ',color_names{colore},'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\trimming ',channel_names{channel},' ',color_names{colore},'.fig'])
                    end
                end 
            end
            
            % plot 405 and 465 together
            figure; clf; 
            if show_plot == 0
            	set(gcf,'visible','off')
            end
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    plot(time_vect,streams.rawdata.(channel_names{channel}).(color_names{colore})); hold on;
                    min_x=-200; max_x=max(time_vect)+200; min_y=mean(streams.rawdata.(channel_names{channel}).(color_names{colore}))-100; 
                    max_y=mean(streams.rawdata.(channel_names{channel}).(color_names{colore}))+100; axis([min_x max_x min_y max_y]); 
                    xline(timetrim_start,':k','trim'); xline(time_vect(end)-timetrim_end,':k','trim'); 
                    axis([min_x max_x min_y max_y]); %ie same as other subplot
                    xlabel('time(sec)'); ylabel('fluorescence'); title('405 and 465'); 
                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\405 and 465 ',channel_names{channel},' ',color_names{colore},'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\405 and 465 ',channel_names{channel},' ',color_names{colore},'.fig'])
                    end
                end
            end
            
            % Reassign data variables if happy with trimming 
            for channel=1:length(channel_names)
                for colore=1:length(color_names)
                    streams.rawdata.(channel_names{channel}).(color_names{colore}) = streams.rawdata_trim.(channel_names{channel}).(color_names{colore});
                end
            end
            time_vect = time_vect_trim;
            for i=1:length(epoc_type)
                epocs.(epoc_type{i}) = epocs_trim.(epoc_type{i});
            end
            
     
            length_data = length(streams.rawdata.(channel_names{channel}).(color_names{colore}));
            clear epocs_trim time_vect_trim newPos_trim new_Pos
            streams = rmfield(streams,'rawdata_trim'); %clear variable within structure.


            %% LOW PASS FILTER OF FP DATA 
            if lowpassfiltering == 1
                %filter option 2: butter
                ftype = 'low';
                n = 2; % 2nd order filter. Will use filtfilt so it's fwd and backward filter which leads to minimal data distorsion
                Wn = lowpassfreq/((sampling_rate_ds)/2); %lowpassfreq defined above.
                % 0.5 Hz = 2 sec ; 1 Hz = 1 sec ; 2 Hz = 0.5 sec ; 3 Hz = 0.33 sec
                [a,b] = butter(n,Wn,ftype);
                for channel=1:length(channel_names)
                    for colore=1:length(color_names)
                        streams.lowfilt.(channel_names{channel}).(color_names{colore}) = filtfilt(a,b,double(streams.rawdata.(channel_names{channel}).(color_names{colore})));
                        % plot
                        figure; clf; 
                        if show_plot == 0
                           set(gcf,'visible','off')
                        end
                        plot(time_vect,streams.rawdata.(channel_names{channel}).(color_names{colore})); hold on; 
                        plot(time_vect,streams.lowfilt.(channel_names{channel}).(color_names{colore}));
                        min_x=-10; max_x=max(time_vect)+10; min_y=min(streams.rawdata.(channel_names{channel}).(color_names{colore}))-1; 
                        max_y=max(streams.rawdata.(channel_names{channel}).(color_names{colore}))+1; axis([min_x max_x min_y max_y]); 
                        xlabel('samples'); ylabel('fluorescence'); title([channel_names{channel},' ', color_names{colore},' data before and after filtering']);
                        L=legend('data -10','data after filter'); L.Location = 'Best';
                        if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                            saveas(gcf,[PATH2SAVE,'figures\lowpass ',channel_names{channel},' ',color_names{colore},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\lowpass ',channel_names{channel},' ',color_names{colore},'.fig'])
                        end
                    end
                end

                %% Save data if happy with low-pass filtering
                for channel=1:length(channel_names)
                    for colore=1:length(color_names)
                        streams.rawdata.(channel_names{channel}).(color_names{colore}) = streams.lowfilt.(channel_names{channel}).(color_names{colore});
                    end
                end
                streams = rmfield(streams,'lowfilt'); %clear variable within structure.
            end %of low pass filtering

            %% Plot 405 and 465 together
           
            % VEH injection and DETQ times
            if contains(path2data,'2022')
                time_drug1.VEH=epocs.Drug_time(1)+60; % estimate since i press on the button right after injection and it takes 60 sec for things to stabilize
                time_drug1.DETQ=epocs.Drug_time(2)+60;
            else
                time_drug1.DRUG=epocs.Drug_time(1)+60; % estimate since i press on the button after injection
            end
            
            figure; clf; 
            if show_plot == 0
            	set(gcf,'visible','off')
            end
            for channel=1:length(channel_names)
%                 for colore=1:length(color_names)
                    plot(time_vect,streams.rawdata.(channel_names{channel}).c405,'Color','k','Linewidth',1); hold on;
                    plot(time_vect,streams.rawdata.(channel_names{channel}).c465,'Color','r','Linewidth',1); hold on;

                    for o = 1:size(epocs.PulseStart_time)
                        if o ~= 1 && o ~= 7
                            plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[55 60],'Color',[15, 128, 242]/255,'Linewidth',1)
                            plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[55 60],'Color',[15, 128, 242]./255,'Linewidth',1)
                        end
                    end
                    text(5630,49,'opto stim','FontSize',16,'Color',[15, 128, 242]/255)

                    min_x=-200; max_x=7200; min_y=40; max_y=120; axis([min_x max_x min_y max_y]); 
                    axis([min_x max_x min_y max_y]); %ie same as other subplot
                    
                    if contains(path2data,'2022')                        
                        xline(time_drug1.VEH-60,':k');  xline(time_drug1.VEH+120,':k','Vehicle','FontSize',16); % grey out the noise in signal due to injection 
                        xline(time_drug1.DETQ-60,':k'); xline(time_drug1.DETQ+120,':k','DETQ','FontSize',16); 
                        %plot the rotarod area
                        error_area_onlyrectangle([time_drug1.VEH-60 time_drug1.VEH+120],[40 120],[80 80],[165, 167, 168]./255,0.8); 
                        error_area_onlyrectangle([time_drug1.DETQ-60 time_drug1.DETQ+120],[40 120],[80 80],[165, 167, 168]./255,0.8); 
                    else
                        xline(time_drug1.DRUG-60,':k');  xline(time_drug1.DRUG+120,':k','Drug','FontSize',16); 
                        %plot the rotarod area
                        error_area_onlyrectangle([time_drug1.DRUG-60 time_drug1.DRUG+120],[40 120],[80 80],[165, 167, 168]./255,0.8); 
                    end
                    
                    % axis
                    set(gca,'FontSize',16,'FontName', 'Arial')  
                    ax = gca;
                    set(gca,'TickLength',[0 0])
                    
                    box off
                    xlabel('Time (s)'); ylabel('Raw fluorescence');  
                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\405 and 465 representative',channel_names{channel},' ',color_names{colore},'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\405 and 465 representative',channel_names{channel},' ',color_names{colore},'.fig'])
%                     end
                end
            end

            %% CALCULATE dFF (synapse and polyfit versions)
            % Initialize
            F0_405 = []; F0_465 = []; F0_560 = []; F0_405_2 = []; F0_465_2 = [];
            dFF_405 = []; dFF_465 = []; dFF_560 = []; dFF_405_2 = []; dFF_465_2 = [];
            % Calculations in each channel and color
            % First channel, 465 color, first dFF
            channel=1;
            % dFF polyfit
            % Step 1: Polyfit on entire data 
            index = 1:length(streams.rawdata.(channel_names{channel}).c465);
            calc_coeff_fit_465 = polyfit(streams.rawdata.(channel_names{channel}).c405,streams.rawdata.(channel_names{channel}).c465,1); %fitting control Dv2 against signal Dv1 using least squares fit of degree 1
            % Step 2: Apply polyfit on 465 
            data405_fitted_465 = calc_coeff_fit_465(1).*streams.rawdata.(channel_names{channel}).c405 + calc_coeff_fit_465(2);
            % Step 3: Normalize
            streams.dFF.(dFF_names{1}) = 100*((streams.rawdata.(channel_names{channel}).c465 - data405_fitted_465)./data405_fitted_465); % of deltaF/F              

            % second dFF, either 2nd color in channel 1 or 2nd channel (1 color)- if ever I had 2 color, 2 channels, would need to edit
            % this script
            if color_number == 3 && channel_number == 1
                channel = 1;
                % % dFF polyfit
                % Step 1: Polyfit on entire data
                index = 1:length(streams.rawdata.(channel_names{channel}).c560);
                calc_coeff_fit_560 = polyfit(streams.rawdata.(channel_names{channel}).c405,streams.rawdata.(channel_names{channel}).c560,1); %fitting control Dv2 against signal Dv1 using least squares fit of degree 1
                % Step 2: Apply polyfit on 560 
                data405_fitted_560 = calc_coeff_fit_560(1).*streams.rawdata.(channel_names{channel}).c405 + calc_coeff_fit_560(2);
                % Step 3: Normalize
                streams.dFF.(dFF_names{2}) = 100*((streams.rawdata.(channel_names{channel}).c560 - data405_fitted_560)./data405_fitted_560); % of deltaF/F 
            elseif color_number == 2 && channel_number == 2
                channel = 2;
                % % dFF polyfit
                % Step 1: Polyfit on entire data 
                index = 1:length(streams.rawdata.(channel_names{channel}).c465);
                calc_coeff_fit_465_2 = polyfit(streams.rawdata.(channel_names{channel}).c405,streams.rawdata.(channel_names{channel}).c465,1); %fitting control Dv2 against signal Dv1 using least squares fit of degree 1
                % Step 2: Apply polyfit on 465 
                data405_fitted_465_2 = calc_coeff_fit_465_2(1).*streams.rawdata.(channel_names{channel}).c405 + calc_coeff_fit_465_2(2);
                % Step 3: Normalize
                streams.dFF.(dFF_names{2}) = 100*((streams.rawdata.(channel_names{channel}).c465 - data405_fitted_465_2)./data405_fitted_465_2); % of deltaF/F  
            end

            % plot all
            for d=1:length(dFF_names)
                figure; clf; 
                if show_plot == 0
                   set(gcf,'visible','off')
                end
                yyaxis left
                plot(time_vect,streams.dFF.(dFF_names{d})); hold on; 
                ylabel('dFF');
                min_y=min(streams.dFF.(dFF_names{d}))-3; max_y=max(streams.dFF.(dFF_names{d}))+3; axis([min_x max_x min_y max_y]);
                yyaxis right
                for o = 1:size(epocs.PulseStart_time)
                    plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-2 2],'m')
                    plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-2 2],'c')
                end
                xlabel('time (sec)'); ylabel('dFF'); title([dFF_names{d},' Normalized data']);
                L=legend('dFF','dFF Synapse','Stim on','Stim off'); L.Location = 'Best';
                min_x=-30; max_x=max(time_vect)+30; min_y=min(streams.dFF.(dFF_names{d}))-3; max_y=max(streams.dFF.(dFF_names{d}))+3; axis([min_x max_x min_y max_y]); 
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVE,'figures\dFF calculations ',dFF_names{d},'.tif'])
                    saveas(gcf,[PATH2SAVE,'figures\dFF calculations ',dFF_names{d},'.fig'])
                end
            end
            
            % plot dFF (normal one) shifted up
            deltashift = 25;
            length_data2 = length(streams.dFF.(dFF_names{d})); %will just take it from the last data
            % time vector
            max_idx2 = length_data2;
            sampling_rate_dmin = sampling_rate_ds * 60; % sampling rate in min = 100*60 frames per min ; sampling rate in sec = 100 frames per sec
            dt_dmin = 1/sampling_rate_dmin; % delta t or unit of time (in min) per sample --> for 101.73 Hz sampling rate (after downsampling), will be close to 10 ms (9.8 ms)
            time_vect_min = 0:dt_dmin*1:dt_dmin*(max_idx2-1); %time vector; values in seconds, same number of data points as the data
            
            for d=1:length(dFF_names)
                figure; clf; 
                if show_plot == 0
                   set(gcf,'visible','off')
                end
                plot(time_vect_min,streams.dFF.(dFF_names{d})+deltashift); hold on; 
                ylabel('dFF');
                min_y=min(streams.dFF.(dFF_names{d}))-3; max_y=max(streams.dFF.(dFF_names{d}))+3; axis([min_x max_x min_y max_y]);
                for o = 1:size(epocs.PulseStart_time)
                    plot([epocs.PulseStart_time(o)/60 epocs.PulseStart_time(o)/60],[-12 -8],'m')
                end
                xlabel('time (min)'); title([dFF_names{d},' Normalized data']);
                L=legend('dFF'); L.Location = 'Best';
                min_x=-5; max_x=max(time_vect_min)+5; min_y=min(streams.dFF.(dFF_names{d}))-3+deltashift-10; max_y=max(streams.dFF.(dFF_names{d}))+3+deltashift; axis([min_x max_x min_y max_y]); 
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVE,'figures\dFF shifted up against baseline ',dFF_names{d},'.tif'])
                    saveas(gcf,[PATH2SAVE,'figures\dFF shifted up against baseline ',dFF_names{d},'.fig'])
                end
            end
 
            %% DETRENDING dFF
            if detrending == 1;% %% Detrend dFF, Option 1: normal detrend (other detrend options not here, but could include them later, perhaps as a function)
                for d=1:length(dFF_names)
                    streams.dFF_dtr.(dFF_names{d}) = detrend(streams.dFF.(dFF_names{d})); % This works but doesn't do much
                    %plot
                    figure; clf
                    if show_plot == 0
                       set(gcf,'visible','off')
                    end
                    h1 = plot(time_vect,streams.dFF.(dFF_names{d})); hold on; h2 = plot(time_vect,streams.dFF_dtr.(dFF_names{d}));
                    for o = 1:size(epocs.PulseStart_time)
                        plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-2 2],'m')
                        plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-2 2],'c')
                    end
                    min_x=-30; max_x=max(time_vect)+30; min_y=min(streams.dFF.(dFF_names{d}))-10; max_y=max(streams.dFF.(dFF_names{d}))+10; axis([min_x max_x min_y max_y]); 
                    xlabel('time (sec)'); ylabel('dFF'); title([dFF_names{d},' Normalized data, Regular detrend']);
                    L=legend([h1 h2],'dFF','dFF detrended'); L.Location='Best'; 
                    if done == 0 && save_plot == 1 || done == 0 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\dFF detrending ',dFF_names{d},'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\dFF detrending ',dFF_names{d},'.fig'])
                    end
                    % if happy with detrend:
                    streams.dFF.(dFF_names{d}) = streams.dFF_dtr.(dFF_names{d});
                    streams = rmfield(streams,'dFF_dtr'); %clear variable within structure.
                end
            else
            end

            %% ADJUSTING dFF BASELINE TO 0 IN A MOVING WINDOW      
            if deletemovingwindow == 1;
                for d=1:length(dFF_names)
                    % Define parameters for moving window % This works but should only be used on specific datasets that need be, with careful inspection of the datasets
                    baseline_window = 10; %Adjust this! %Window to make the calculation from eg if (60 sec): you calculate the baseline for the 60 sec, looking at 50 sec ahead of the current 10sec 
                    baseline_dt = 5; %Moving steps for the threshold --> eg move 10 sec - value has to be < baseline_window
                    percentile = 8; %8th percentile
                    % New time vector
                    temp_t = time_vect - time_vect(1); %creating new time vector set to starting at 0, incase time_vect weas not starting at 0
                    dummie = 1:length(temp_t); %generates index of same amount of samples (better than "find" cos too much computing time)
                    dummie = dummie(temp_t >= baseline_window); %value when the time is beyond the window --> so that we can apply the moving window calculations
                    idx_baseline_window = dummie(1); %index of the timestamp as of which we can start the 60sec moving window
                    dummie = 1:length(temp_t);
                    dummie = dummie(temp_t >= baseline_dt);
                    idx_baseline_dt = dummie(1); %index of the timestamp as of which we start the 10sec moving window, i.e. baseline dt
                    clear dummie temp_t
                    % Index for the moving steps
                    ix1 = 1:idx_baseline_dt:length_data-idx_baseline_window; %take whole duration and then look at index you will look at (eg 1 to 60, 1 to 70 etc). calculate it before doing the loop, makes it easier
                    % create index vector for the first bound of the window, where steps are of the size of 10 sec scaled to the
                    % samples we have, going until the end minus the size of the %window
                    ix2 = idx_baseline_window:idx_baseline_dt:length_data; %index for the second bound of the window
                    ix = [ix1' ix2']; %first and last index you will check in each window (first colum onset, second colum offset) and the rows are the windows
                    if ix(end,2) < length_data
                        [ix] = [ix;[length_data-idx_baseline_window length_data]]; %trick to include the last bit
                    end
                    % Calculate 8th percentile across Pre and Post datasets if you have multiple phases of the data, eg before and after the experiment
                    dFF_perc8.(dFF_names{d}) = ones(1,length_data)*nan; %helps to specify how big your variable will be. Better to fill with Nans ahead of time: preallocation of the space
                    dFF_f0.(dFF_names{d}) = ones(1,length_data)*nan;
                    % calculate the entire 8th percentile (=baseline)
                    for o = 1:size(ix,1) %loop goes thru rows, equivalent of number of windows
                        idx = ix(o,1):ix(o,2); %window for this run of the for loop ; new index, saves space later. allows to go 1 to 60 etc.
                        if o == 1 %tricky to treat the first loop the same as other points, eg opint 0 is represented by index 1 (0 to 10sec),sometimes this doesnt work
                            dFF_perc8.(dFF_names{d})(idx) = prctile(streams.dFF.(dFF_names{d})(idx),percentile); %for this first window (idx), assign 8th percentile value of dFF values during this window, 
                        else
                            dFF_perc8.(dFF_names{d})(idx(end-idx_baseline_dt)+1:end) = prctile(streams.dFF.(dFF_names{d})(idx),percentile); 
                        end
                    end
                    %Substract 8th percentile, or calculated based on average to various epochs
                    %Step 1: option 1: substract moving baseline from all 
                    dFF_f0.(dFF_names{d}) = streams.dFF.(dFF_names{d}) - dFF_perc8.(dFF_names{d});

                    %Plot
                    figure; clf
                    if show_plot == 0
                       set(gcf,'visible','off')
                    end
                    h1=plot(time_vect,streams.dFF.(dFF_names{d})-4); hold on; h2=plot(time_vect,dFF_f0.(dFF_names{d})); yline(0,'--','zero'); yline(-4,'--','zero');
                    min_x=-30; max_x=max(time_vect)+30; min_y=min(streams.dFF.(dFF_names{d}))-10; max_y=max(streams.dFF.(dFF_names{d}))+10; axis([min_x max_x min_y max_y]); 
                    xlabel('time (sec)'); ylabel('dFF'); title([dFF_names{d},' Normalized data, 8t percentile removed']);
                    h3=plot(time_vect,dFF_perc8.(dFF_names{d})-4); 
                    for o = 1:size(epocs.PulseStart_time)
                        plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-5 5],'m','linewidth',1)
                        plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-5 5],'c','linewidth',1)
                    end
                    L=legend([h1,h2,h3],'dFF','dFF, baseline corrected','moving baseline','rotarod on/off'); L.Location = 'Best';
                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\dFF percentilecorr ',dFF_names{d},'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\dFF percentilecorr ',dFF_names{d},'.fig'])
                    end
                    %% CONVERT dFF to dFF_f0 corrected
                    streams.dFF.(dFF_names{d}) = dFF_f0.(dFF_names{d});
                    clear dFF_f0
                end
            end
            
            %% HIGH-PASS FILTERING
            if highpassfiltering == 1; %'yes' if you want to high pass filter- typically off for open field
                % Filter out the slow fluctuation
                ftype = 'high';
                n = 2; % 2nd order filter
                Wn = highpassfreq/((sampling_rate_ds)/2); 
                [a,b] = butter(n,Wn,ftype);

                for d=1:length(dFF_names)
                    streams.dFF_hp.(dFF_names{d}) = filtfilt(a,b,double(streams.dFF.(dFF_names{d})));
                    % plot
                    figure; clf
                    if show_plot == 0
                       set(gcf,'visible','off')
                    end
                    yyaxis left
                    h1=plot(time_vect,streams.dFF.(dFF_names{d})); hold on; 
                    for o = 1:size(epocs.PulseStart_time)
                        plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-2 2],'m')
                        plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-2 2],'c')
                    end
                    yyaxis right
                    h2=plot(time_vect,streams.dFF_hp.(dFF_names{d}));
                    min_x=-10; max_x=max(time_vect)+10; min_y=min(streams.dFF.(dFF_names{d}))-10; max_y=max(streams.dFF.(dFF_names{d}))+10; axis([min_x max_x min_y max_y]); 
                    xlabel('samples'); ylabel('fluorescence'); title([dFF_names{d},' dFF before and after filtering']);
                    L=legend([h1,h2],'dFF','dFF after high-pass filter'); L.Location = 'Best';

                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\dFF highpass ',dFF_names{d},'.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\dFF highpass ',dFF_names{d},'.fig'])
                    end
                    % if happy with detrend:
                    streams.dFF.(dFF_names{d}) = streams.dFF_hp.(dFF_names{d});
                    streams = rmfield(streams,'dFF_hp'); %clear variable within structure.
                end
            else
            end


            %% CALCULATE ROBUST Z-SCORES
            % Calculations
            % Calculate median, Median Absolute Deviation (MAD) and regular "modified" Z-score 
            for d=1:length(dFF_names)
                med_dFF.(dFF_names{d}) = median(streams.dFF.(dFF_names{d}));    
                MAD_dFF.(dFF_names{d}) = mad(streams.dFF.(dFF_names{d}),1);     
                streams.ZScoredFF.(dFF_names{d}) = 0.6745*(streams.dFF.(dFF_names{d})-med_dFF.(dFF_names{d}))./MAD_dFF.(dFF_names{d}); 
                % Plot
                figure; clf; 
                if show_plot == 0
                   set(gcf,'visible','off')
                end
                yyaxis left
                h1 = plot(time_vect,streams.dFF.(dFF_names{d}),'LineWidth',0.5,'Color','k'); hold on; 
                min_y = min(streams.dFF.(dFF_names{d}))-10; max_y = max(streams.dFF.(dFF_names{d}))+10; ylim([min_y max_y]);
                ylabel('dFF');
                yyaxis right
                h2 = plot(time_vect,streams.ZScoredFF.(dFF_names{d}), 'LineWidth',0.5,'Color',[1.0 0.2 0.2]);
                h3 = yline(0,'linewidth',2,'linestyle',':');
                for o = 1:size(epocs.PulseStart_time)
                    plot([epocs.PulseStart_time(o) epocs.PulseStart_time(o)],[-2 2],'m')
                    plot([epocs.PulseStop_time(o) epocs.PulseStop_time(o)],[-2 2],'c')
                end
                min_y = -10; max_y = 10; ylim([min_y max_y]);
                L = legend([h1,h2],'dFF','Z-score');
                L.Location = 'northeast';
                min_y = min(streams.ZScoredFF.(dFF_names{d}))-10; max_y = max(streams.ZScoredFF.(dFF_names{d}))+10;
                min_x = -20; max_x = max(time_vect)+20;  
                xlim([min_x max_x])
                xlabel('time (s)'); ylabel('Z-score dFF'); title([dFF_names{d},' dFF and Z-score']);
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVE,'figures\dFF Zscore ',dFF_names{d},'.tif'])
                    saveas(gcf,[PATH2SAVE,'figures\dFF Zscore ',dFF_names{d},'.fig'])
                end
            end
                  
              %% IF WE CALCULATE "dFF WITHIN" THEN ADD AN EMPTY ARRAY, WHICH WE WILL FILL LATER

            if keeprawdata == 1
                for d=1:length(dFF_names)    
                    streams.dFFwithin.(dFF_names{d}) = streams.ZScoredFF.(dFF_names{d})*nan;
                    streams.ZScoredFF_within.(dFF_names{d}) = streams.ZScoredFF.(dFF_names{d})*nan
                end
            end
            
            %% TRIGGER DATA BASED ON OPTOGENETIC STIMULATION (OR OTHER EVENT)
            % Setup the time vectors for each trial based on how big you want the window to be
            dummie = 1:length(time_vect); % indices 1 to the end of the time vector
            dummie = dummie(time_vect >= abs(TRANGE(1))); %select the indices of the time vector starting from the beginning of the time window
            idx_Init = dummie(1); % this was basically a trick to get the number of time vector indices needed to produce a duration of TRANGE(1) i.e. 15 sec
            dummie = 1:length(time_vect); 
            dummie = dummie(time_vect >= abs(TRANGE(2)));
            idx_End = dummie(1);
            n = idx_Init + idx_End; % Length of each trial in indices
            t_trials = time_vect(1:n) - time_vect(idx_Init); % get the time vector for 1 trial, of duration n, starting at initial index

            % Create 2-column array with the event to align to (start and stop)
            Events2Align2 = {'PulseStart_time','PulseStop_time'}                 
            t_opto = [epocs.(Events2Align2{1}) epocs.(Events2Align2{1})+2]; % t_opto: a double column with all trials: column 1 (all rows) the start of the opto events and column 2 (all rows) the end
            if t_opto(end,2) == Inf %eg if trial was not finished
                t_opto(length(t_opto),:) = [];
            end

            % define the t_opto2 matrix that will contain the individual trials for each trial type (rows) - onset/offset (2 columns)                
            if contains(path2data,'_dLight_2022') 
                 for d=1:length(dFF_names)
                    t_opto2.(dFF_names{d}).VEH = t_opto(VEH_trial1:VEH_trialend,:); % t_opto2: group the ones to average from t_opto 
                    t_opto2.(dFF_names{d}).DETQ = t_opto(DETQ_trial1:DETQ_trialend,:); % t_opto2: group the ones to average from t_opto 
                 end
                 time_drug.VEH=epocs.Drug_time(1)+60;
                 time_drug.DETQ=epocs.Drug_time(2)+60; % this is the approximate time at which we were done injecting the drug
            else
                 for d=1:length(dFF_names)
                     if contains(path2data,'00')
                        t_opto2.(dFF_names{d}).DRUG = t_opto(VEH_trial1:VEH_trialend,:); % t_opto2: group the ones to average from t_opto 
                     else
                        t_opto2.(dFF_names{d}).DRUG = t_opto(DETQ_trial1:DETQ_trialend,:); % t_opto2: group the ones to average from t_opto                          
                     end
                 end
                 time_drug.DRUG=epocs.Drug_time(1)+60;
            end
            
            % keep track of time for opto trials 
            if contains(path2data,'_dLight_2022')
                tmp = round(round(t_opto2.(dFF_names{d}).VEH(:,1) - time_drug.VEH)./60);
                if tmp(1) == 1 || tmp(1) == 2 || tmp(1) == -1 || tmp(1) == - 2 ; tmp = tmp - tmp(1); end;
                if tmp(1) == 0 && tmp(2) == 6 ; tmp = tmp - 1; tmp(1) = 0; end;                
                time_optotrials_double.VEH = tmp;
                tmp = round(round(t_opto2.(dFF_names{d}).DETQ(:,1) - time_drug.DETQ)./60);
                if tmp(1) == 1 || tmp(1) == 2 || tmp(1) == -1 || tmp(1) == - 2; tmp = tmp - tmp(1); end;
                time_optotrials_double.DETQ = tmp;
            else
                tmp = round(round(t_opto2.(dFF_names{d}).DRUG(:,1) - time_drug.DRUG)./60);
                if tmp(1) == 1 || tmp(1) == 2 || tmp(1) == -1 || tmp(1) == - 2; tmp = tmp - tmp(1); end;
                if tmp(1) == 0 && tmp(2) == 2 ; tmp = tmp + 3; tmp(1) = 0; end;                
                time_optotrials_single = tmp;            
            end
 
            %% Generate the split dFF data (and other datatypes) trial matrixes based on trial allocation (t_opto2)
            % datatype = {'ZScoredFF','dFF'} or other, are defined up in the parameters
            for d=1:length(dFF_names)
                % Compute the aligned results of the trials, initialization to have them the number of rows that are in t_opto2 (=trials) and the length of the window: n
                for pow = 1:length(Opto_Powers)
                    for j = 1:length(datatype)
                        IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{j}) = ones(size(t_opto2.(dFF_names{d}).(Opto_Powers{pow}),1),n)*nan;
                        IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{j}) = ones(size(t_opto2.(dFF_names{d}).(Opto_Powers{pow}),1),n)*nan;
                    end
                end
                
                % Compute the aligned results of the trials, assign the real data
                for pow = 1:length(Opto_Powers)
                    for i = 1:length(datatype)
                        for o = 1:size(t_opto2.(dFF_names{d}).(Opto_Powers{pow}),1) % for each row, ie each trial in t_opto2 for this trialtype
                            ix = find(abs(time_vect-t_opto2.(dFF_names{d}).(Opto_Powers{pow})(o,1)) == min(abs(time_vect-t_opto2.(dFF_names{d}).(Opto_Powers{pow})(o,1))));
                            % ix is the closest (just above) index in time_vect at which the trial starts. to find it, you look at the index where time_vect-t_opto2 is minimal
                            if ix+idx_End <= length(time_vect) && ix-(idx_Init-1) >= 1 %if all good, index+end of window or index-beg of window within recording time
                                tmp = ix - (idx_Init-1):ix + idx_End;
                                IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(o,:) = streams.(datatype{i}).(dFF_names{d})(tmp);
                            elseif ix+idx_End > length(time_vect) % if you fall off the end of recording
                                tmp = ix - (idx_Init-1):length(time_vect); %that means the rest will remain NaNs
                                IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(o,1:length(tmp)) = ...
                                        streams.(datatype{i}).(dFF_names{d})(tmp);
                            elseif ix-(idx_Init-1) < 1 % if you fall off at the beginning of recording
                                tmp = 1:(ix + idx_End); %that means the beginning will remain NaNs
                                IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(o,(n-length(tmp)+1):end) = ...
                                        streams.(datatype{i}).(dFF_names{d})(tmp);
                            end
                        end
                        % baseline correction (for each row, (:) normalize against the median of the values of dFF that are in the baseline portion that you decided about
                        % earlier, ie for each value of the time vector bigger than baselinewindow start and smaller than baselinewindow end 
                        IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})= ...
                        IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{i}) - nanmedian...
                        (IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})...
                        (:,t_trials >= BASELINE_WIN.OptoStim(1) & t_trials <= BASELINE_WIN.OptoStim(2)),2); %t_trials is the time vector for each trial
                        % M = median(A,dim) returns the median of elements along dimension dim. 
                        % For example, if A(:,[start end]) is a matrix, then median(A(:,[start end]),2) is a column vector containing the median value of each row.
                        % so nanmedian will be a column with a median for each row, which you substract, row-by-row to a matrix with same amount of rows
                    end
                end                               
            end
            
           %% Generate matrix with the raw data aligned to the first opto stim after the drug injection minus 60 sec in order to generate graphs over time
            if contains(path2data,'_dLight_2022') 
                drug_trials_1 = [VEH_trial1:VEH_trialend];
                drug_trials_2 = [DETQ_trial1:DETQ_trialend];
            else
                drug_trials_1 = [VEH_trial1:VEH_trialend];
                drug_trials_2 = [DETQ_trial1:DETQ_trialend];                
            end
            
            for d=1:length(dFF_names)
                for k=1:length(datatype)
                    for pow = 1:length(Opto_Powers)
                        streams_aligned.(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}) = ones(1,round(((max([length(drug_trials_1),length(drug_trials_2)])-1)*(stim_interval*60)+extratime_postlastopto)./dt_ds))*nan; %size(streams.(datatype{k}).(dFF_names{d}),2)
                        time_vect_aligned.(Opto_Powers{pow}) = ones(1,size(streams_aligned.(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),2))*nan;
                    end
                end
            end
            % 2022
            if contains(path2data,'_dLight_2022')                 
                for d=1:length(dFF_names)
                    for k=1:length(datatype)
                        for pow = 1:length(Opto_Powers)
                            firstopto_after_drug = round((t_opto2.(dFF_names{d}).(Opto_Powers{pow})(1,1)-extratime_prefirstopto)./dt_ds); % this means we only take the first from the trials we took in analysis 
                            
                            if pow==1 % VEH
                                lastdatapoint = round(firstopto_after_drug + round(((length(drug_trials_1)-1)*(stim_interval*60)+extratime_postlastopto)./dt_ds));
                            else % DETQ
                                lastdatapoint = round(firstopto_after_drug + round(((length(drug_trials_2)-1)*(stim_interval*60)+extratime_postlastopto)./dt_ds));
                            end
                            if lastdatapoint - length(streams.(datatype{k}).(dFF_names{d})) > 0
                                delta = lastdatapoint - length(streams.(datatype{k}).(dFF_names{d}));
                                streams_aligned.(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(1,1:lastdatapoint-firstopto_after_drug+1-delta-1) = streams.(datatype{k}).(dFF_names{d})(firstopto_after_drug:lastdatapoint-delta-1);
                                time_vect_aligned.(Opto_Powers{pow})(1,1:lastdatapoint-firstopto_after_drug+1-delta-1) = time_vect(firstopto_after_drug:lastdatapoint-delta-1);                       
                            else
                                streams_aligned.(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(1,1:lastdatapoint-firstopto_after_drug+1) = streams.(datatype{k}).(dFF_names{d})(firstopto_after_drug:lastdatapoint);
                                time_vect_aligned.(Opto_Powers{pow})(1,1:lastdatapoint-firstopto_after_drug+1) = time_vect(firstopto_after_drug:lastdatapoint);                       
                            end
                        end
                    end
                end
                
            % 2021
            else
                for d=1:length(dFF_names)
                    for k=1:length(datatype)
                        for pow = 1:length(Opto_Powers)
                            firstopto_after_drug = round((t_opto2.(dFF_names{d}).(Opto_Powers{pow})(1,1)-extratime_prefirstopto)./dt_ds); % this means we only take the first from the trials we took in analysis
                            if contains(path2data,'00') % VEH
                                lastdatapoint = round(firstopto_after_drug + round(((length(drug_trials_1)-1)*(stim_interval*60)+extratime_postlastopto)./dt_ds));
                            else % DETQ
                                lastdatapoint = round(firstopto_after_drug + round(((length(drug_trials_2)-1)*(stim_interval*60)+extratime_postlastopto)./dt_ds));
                            end
                            if lastdatapoint - length(streams.(datatype{k}).(dFF_names{d})) > 0
                                delta = lastdatapoint - length(streams.(datatype{k}).(dFF_names{d}));
                                streams_aligned.(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(1,1:lastdatapoint-firstopto_after_drug+1-delta-1) = streams.(datatype{k}).(dFF_names{d})(firstopto_after_drug:lastdatapoint-delta-1);
                                time_vect_aligned.(Opto_Powers{pow})(1,1:lastdatapoint-firstopto_after_drug+1-delta-1) = time_vect(firstopto_after_drug:lastdatapoint-delta-1);                       
                            else
                                streams_aligned.(datatype{k}).(dFF_names{d}).(Opto_Powers{pow})(1,1:lastdatapoint-firstopto_after_drug+1) = streams.(datatype{k}).(dFF_names{d})(firstopto_after_drug:lastdatapoint);
                                time_vect_aligned.(Opto_Powers{pow})(1,1:lastdatapoint-firstopto_after_drug+1) = time_vect(firstopto_after_drug:lastdatapoint);                       
                            end
                        end
                    end
                end                
            end
            
            % graph
            color2plot = {'b','g','m','r'};

            for d=1:length(dFF_names)
                for k=1; %:length(datatype)
                    figure;
                    for pow = 1:length(Opto_Powers)
                        plot(time_vect_aligned.(Opto_Powers{pow}),streams_aligned.(datatype{k}).(dFF_names{d}).(Opto_Powers{pow}),'Color',(color2plot{pow}),'Linewidth',1); hold on;
                    end
%                     min_x=-200; max_x=max(time_vect_aligned.(Opto_Powers{2})); min_y=-10; max_y=25; axis([min_x max_x min_y max_y]); 
               % axis
                    set(gca,'FontSize',16,'FontName', 'Arial')  
                    ax = gca;
                    set(gca,'TickLength',[0 0])
                    if contains(path2data,'2022')
                        L=legend(drug_name{1},drug_name{2}); L.Location = 'Best';
                    else
                        L=legend(drug_name{1}); L.Location = 'Best';                        
                    end
                    box off
                    xlabel('Time (s)'); ylabel([datatype{k}]);  
                    if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                        saveas(gcf,[PATH2SAVE,'figures\streams aligned',datatype{k},' ',dFF_names{d},' ','.tif'])
                        saveas(gcf,[PATH2SAVE,'figures\streams aligned',datatype{k},' ',dFF_names{d},' ','.fig'])
        %                     end
                    end
                end
            end
               
                %% Version with calculating dFF within the individual trials
            if keeprawdata == 1
                for d=1:length(dFF_names)
                % Compute the aligned results of the trials, initialization to have them the number of rows that are in t_opto2 (=trials) and the length of the window: n
                    for pow = 1:length(Opto_Powers)
                        IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).dFFwithin = ones(size(t_opto2.(dFF_names{d}).(Opto_Powers{pow}),1),n)*nan;
                        IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).dFFwithin = ones(size(t_opto2.(dFF_names{d}).(Opto_Powers{pow}),1),n)*nan;
                    end
                end
                
             % Compute the aligned results of the trials, assign the real data
                for d=1:length(dFF_names)
                    for pow = 1:length(Opto_Powers)
                        for o = 1:size(t_opto2.(dFF_names{d}).(Opto_Powers{pow}),1) % for each row, ie each trial in t_opto2 for this trialtype
                            ix = find(abs(time_vect-t_opto2.(dFF_names{d}).(Opto_Powers{pow})(o,1)) == min(abs(time_vect-t_opto2.(dFF_names{d}).(Opto_Powers{pow})(o,1))));
                            % ix is the closest (just above) index in time_vect at which the trial starts. to find it, you look at the index where time_vect-t_opto2 is minimal
                            if ix+idx_End <= length(time_vect) && ix-(idx_Init-1) >= 1 %if all good, index+end of window or index-beg of window within recording time
                                tmp = ix - (idx_Init-1):ix + idx_End;
                                % calculate dFF
                                d405 = streams.rawdata.(dFF_names{d}).c405(tmp);
                                d465 = streams.rawdata.(dFF_names{d}).c465(tmp);
                                %fit the signal only during the baseline (before the event)
                                d405_pre=d405(1:idx_Init-1);
                                d465_pre=d465(1:idx_Init-1);  
                                calc_coeff_fit_465 = polyfit(d405_pre,d465_pre,1);
                                data405_fitted_465 = calc_coeff_fit_465(1).*d405 + calc_coeff_fit_465(2);                                
                                IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).dFFwithin(o,:) = 100*(d465 - data405_fitted_465)./data405_fitted_465;
                            elseif ix+idx_End > length(time_vect) % if you fall off the end of recording
                                tmp = ix - (idx_Init-1):length(time_vect); %that means the rest will remain NaNs
                                % calculate dFF
                                d405 = streams.rawdata.(dFF_names{d}).c405(tmp);
                                d465 = streams.rawdata.(dFF_names{d}).c465(tmp);
                                %fit the signal only during the baseline (before the event)
                                d405_pre=d405(1:idx_Init-1);
                                d465_pre=d465(1:idx_Init-1);  
                                calc_coeff_fit_465 = polyfit(d405_pre,Yd465_pre,1);
                                data405_fitted_465 = calc_coeff_fit_465(1).*d405 + calc_coeff_fit_465(2);                                
                                IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).dFFwithin(o,1:length(tmp)) = 100*(d465 - data405_fitted_465)./data405_fitted_465;
                            elseif ix-(idx_Init-1) < 1 % if you fall off at the beginning of recording
                                tmp = 1:(ix + idx_End); %that means the beginning will remain NaNs
                                % calculate dFF
                                d405 = streams.rawdata.(dFF_names{d}).c405(tmp);
                                d465 = streams.rawdata.(dFF_names{d}).c465(tmp);
                                %fit the signal only during the baseline (before the event)
                                d405_pre=d405(1:idx_Init-1);
                                d465_pre=d465(1:idx_Init-1);  
                                calc_coeff_fit_465 = polyfit(d405_pre,Yd465_pre,1);
                                data405_fitted_465 = calc_coeff_fit_465(1).*d405 + calc_coeff_fit_465(2);                                
                                IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).dFFwithin(o,(n-length(tmp)+1):end) = 100*(d465 - data405_fitted_465)./data405_fitted_465;
                            end
                            % baseline correction (for each row, (:) normalize against the median of the values of dFF that are in the baseline portion that you decided about
                            % earlier, ie for each value of the time vector bigger than baselinewindow start and smaller than baselinewindow end 
                            IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).dFFwithin= ...
                            IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).dFFwithin - nanmedian...
                            (IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).dFFwithin...
                            (:,t_trials >= BASELINE_WIN.OptoStim(1) & t_trials <= BASELINE_WIN.OptoStim(2)),2); %t_trials is the time vector for each trial
                        end
                    end
                end
            end
            
               
            
            %% Calculate the Zscore 'dFFwithin'. The code I tried, actually did not work, so artificially we are now assigning this Zscore-within to ZScore, to not induce errors in the code. 
            % Won't be in use.
            pooledtype = {'raw','baselinecorr'};
            for d=1:length(dFF_names)
                for pow = 1:length(Opto_Powers)
                    for p=1:length(pooledtype)
                        IndivStim_data.(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).ZScoredFF_within = IndivStim_data.(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).ZScoredFF;
                    end
                end
            end
           
            
            %% Plot the results per trial types separately
            trial_mode = {'raw','baseline_corrected'}; % raw or baseline_corrected.
            
            for d=1:length(dFF_names)
                % Individual trials
                avg_mode = 1; % 1, use mean and SEM; 2, use median and MAD
                for i = 1:length(datatype)
                    for r = 1:length(raw_or_corr)
                        figure;        
                        if show_plot == 0
                           set(gcf,'visible','off')
                        end
%                         scrsz = get(0,'ScreenSize'); %if you want screen size graphs
%                         set(gcf,'Position',scrsz);
                        for pow = 1:length(Opto_Powers)
                            subplot(length(Opto_Powers),1,pow); 
                            plot(t_trials,IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})); % plot individual trials
                            hold on
                            if avg_mode == 1
                                plot(t_trials,nanmean(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i}),1),'b',...
                                    'LineWidth',1.5) % plot the average on top of the individual trials
                            elseif avg_mode == 2
                                plot(t_trials,nanmedian(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i}),1),'b',...
                                    'LineWidth',1.5)
                            end
                            xline(0,'-k');
                            xline(stim_duration,'-m');
                            yline(0,'-.k');
                            if ~isempty(stim_duration2)
                                xline(stim_duration2,'-c');
                            end
                            xlabel('Time (s)')
                            ylabel(datatype{i});
                            xlim([TRANGE(1) TRANGE(2)])
                            ylim([-5 8])
                            title([Opto_Powers{pow}(1:end)])
                        end
                        sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d}]); % for groups of subplots
                        if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                            saveas(gcf,[PATH2SAVE,'figures\INDIV TRIALS ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\INDIV TRIALS ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.fig'])
                        end
                    end
                end
            end %for d=1:length(dFFnames)
                
            %% Plot the results per trial types separately over time - VERTICAL

            startgraph = TRANGE(1);
            endgraph = TRANGE(2);
            topgraphfactor = 2;
            dottedlines = [1:1:1];
            colormap cool
            for d=1:length(dFF_names)
                avg_mode = 1; % 1, use mean and SEM; 2, use median and MAD
                for i = 2; %1:length(datatype)
                    for r = 2; %1:length(raw_or_corr)
%                         scrsz = get(0,'ScreenSize'); %if you want screen size graphs
%                         set(gcf,'Position',scrsz);
                        for pow = 1:length(Opto_Powers)
                            scalingfactor = -10;
                               
                            data2plot=IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i});
                            figure; set(gcf, 'Position',  [100+(400*pow), 70, 400, 45*size(data2plot,1)])
                            if show_plot == 0
                                set(gcf,'visible','off')
                            end
                            for q=1:size(data2plot,1)
                                plot(t_trials,data2plot(q,:)+q*scalingfactor,'LineWidth',2);   
                                yline(q*scalingfactor,'-.k');
                                time2plot = round((t_opto2.(dFF_names{d}).(Opto_Powers{pow})(q,1)-time_drug.(Opto_Powers{pow}))/60);
                                text(startgraph-12,q*(scalingfactor),{time2plot},'Color','black','FontSize',8);
                                text(startgraph-7,q*(scalingfactor),'min','Color','black','FontSize',8);
                                hold on
                            end
                            xline(0,'-k');
                            xline(stim_duration,'-m');
%                             yline(0,'-.k');
                            if ~isempty(stim_duration2)
                                xline(stim_duration2,'-c');
                            end
                            xlabel('Time (s)')
                            ylabel(datatype{i});
                            set(gca,'YTick',[]);
                            set(gca,'Yticklabel',[]);
                            set(gca,'YLabel',[]);
                            ax=gca; ax.YAxis.Color = 'w';
                            set(gca,'Box','off');
                            xlim([startgraph endgraph])
%                             ylim([size(data2plot,1)*scalingfactor+(scalingfactor*topgraphfactor) 0])
                            yline(0,'-.k');
                            for u=1:length(dottedlines)-1
                                yline(size(data2plot,1)*scalingfactor+(scalingfactor*dottedlines(u)),'-.k');
                            end
                        sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},' ',Opto_Powers{pow}]); % for groups of subplots
                        end
                        if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                            saveas(gcf,[PATH2SAVE,'figures\INDIV TRIALS OVER TIME ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\INDIV TRIALS OVER TIME',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.fig'])
                        end
                    end
                end
            end
 
                
            %% Plot the results per trial types over time - only some trials not all (specify below), overlay on same graph           
            % defining colors
            timesections_min_str = {};
            colormap cool           
            color2plot = {[0 1 0],...
                        [1 0 0],...
                        [1 0 1],...
                        [1 1 0],...
                        [0.4660, 0.9740, 0.9880],[0.6500, 0.4250, 0.0980],[0.9290, 0.6940, 0.1250],...
                        [0.4940, 0.1840, 0.5560],[0, 0.2470, 0.9410],[0.1010, 0.7450, 0.9830],...
                        [0.0350, 0.0780, 0.1840],[0.8, 0.25, 0.75],[0.75, 0.75, 0.2],...
                        [0.25, 0.85, 0.05],[0, 0, 1],[0, 0.5, 0],[1, 0, 0],...
                        [0.4660, 0.8740, 0.1880],[0.9500, 0.1250, 0.0980],[0.9290, 0.6940, 0.1250],...
                        [0.4940, 0.1840, 0.5560],[0, 0.2470, 0.9410],[0.1010, 0.7450, 0.9830],...
                        [0.0350, 0.0780, 0.1840],[0.8, 0.25, 0.75],[0.75, 0.75, 0.2],...
                        [0.4660, 0.8740, 0.1880],[0.9500, 0.1250, 0.0980],[0.9290, 0.6940, 0.1250],...
                        [0.4660, 0.8740, 0.1880],[0.1500, 0.1250, 0.8980],[0.9290, 0.6940, 0.1250],...
                        [0.4660, 0.0740, 0.8880],[0.9500, 0.1250, 0.0980],[0.9290, 0.6940, 0.1250]};
            for d = 1; %1:length(dFF_names)
                for i = 2; %1:length(datatype)
                    for r = 2; %1:length(raw_or_corr)
                        for pow = 1:length(Opto_Powers)                        
                            figure;   
                            if show_plot == 0
                               set(gcf,'visible','off')
                            end
                            % defining time
                            if contains(path2data,'2022')
                                if pow == 1 % VEH
                                    timesections_min = [5 10 15 20 25]; %min
                                    ixt1 = (timesections_min(1)/5)+1-(VEH_trial1-1);
                                    ixt2 = (timesections_min(end)/5)+1-(VEH_trial1-1);
                                else % DETQ
                                    timesections_min = [5 10 15 20 25 30 35 40 45 50 55 60]; %min
                                    ixt1 = (timesections_min(1)/5)+1+VEH_trialend-(DETQ_trial1-1);
                                    ixt2 = (timesections_min(end)/5)+1+VEH_trialend-(DETQ_trial1-1);                                
                                end
                            else
                                if contains(path2data,'00') % VEH
                                    timesections_min = [5 10 15 20 25]; %min
                                    ixt1 = (timesections_min(1)/5)+1-(VEH_trial1-1);
                                    ixt2 = (timesections_min(end)/5)+1-(VEH_trial1-1);
                                else % DETQ
                                    timesections_min = [5 10 15 20 25 30 35 40 45 50 55 60]; %min
                                    ixt1 = (timesections_min(1)/5)+1-(DETQ_trial1-1);
                                    ixt2 = (timesections_min(end)/5)+1-(DETQ_trial1-1);                                
                                end
                            end
                            clear timesections_min_str
                            for o=1:length(timesections_min)
                                timesections_min_str(o) = {[num2str(timesections_min(o)),' min']};
                                timesections_min_str(o) = {[num2str(timesections_min(o)),' min']};
                            end
                            timesections = timesections_min*60;
                            
                            % defining data to plot
                            data2use=IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(ixt1:ixt2,:);
                            for q=1:size(data2use,1)
                                plot(t_trials,data2use(q,:),'Color',color2plot{q}); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM
                                hold on
                            end
                            
                            %Formatting                            
                            xline(0,'-k');
                            xline(stim_duration,'-m');
                            xlabel('Time (s)')
                            ylabel('dFF (%)')
                            ax=gca; 
                            set(gca,'Box','off');
                            xlim([TRANGE(1) TRANGE(2)])
                            ylim([-1 15])

                            %legend
                            hleg = legend(timesections_min_str,'FontSize',8,'TextColor','black','Location','northeast','NumColumns',1);                                
                            legend('boxoff')
                            sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},' ',Opto_Powers{pow}]); % for groups of subplots
                            if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                                saveas(gcf,[PATH2SAVE,'figures\INDIV TRIALS OVER TIME AVE v2',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.tif'])
                                saveas(gcf,[PATH2SAVE,'figures\INDIV TRIALS OVER TIME AVE v2',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.fig'])
                            end
                        end
                    end
                end
            end  
                            
            %% Average plot of all individual trials per type of trials with all types of trials on same graph; selected trials 10-20min postdrug to get an idea
            for d=1:length(dFF_names)
                %RAW
                color2plot = {'b','g','m','r'};
                for i=1:length(datatype)
                    for r=2; %1:length(raw_or_corr)
                        figure; 
                        if show_plot == 0
                            set(gcf,'visible','off')
                        end
                        hold on
                        for pow = 1:length(Opto_Powers)
                            if pow == 1; ix=3:5; else ix=3:5; end % trials 10 to 20 min, i.e. trials 3 to 5
                            tmp_avg = nanmean(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(ix,:),1);
                            tmp_error = nanstd(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(ix,:),1,1)./...
                                sqrt(sum(~isnan(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(:,1))));
                            error_area(t_trials,tmp_avg,tmp_error,color2plot{pow},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM
                        end
                        xline(0,'-k'); 
                        xline(stim_duration,'-m');
                        yline(0,'-.k');
                        xlabel('Time (s)');
                        ylabel(datatype{i});
                        xlim([TRANGE(1) TRANGE(2)]);
                        sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},' select trials 10-20min postdrug'],'Interpreter','none')
                        if contains(path2data,'2022')
                        legend([Opto_Powers{1}(1:end),' error'],[Opto_Powers{1}(1:end)],...
                                [Opto_Powers{2}(1:end),'error'],[Opto_Powers{2}(1:end)],...
                                'Opto On','Opto off','Location','northwest','NumColumns',1)      
                        else
                        legend([Opto_Powers{1}(1:end),' error'],[Opto_Powers{1}(1:end)],...
                                'Opto On','Opto off','Location','northwest','NumColumns',1)                              
                        end
                        if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                            saveas(gcf,[PATH2SAVE,'figures\AVE all powers optostim ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\AVE all powers optostim ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.fig'])
                        end
                    end
                end
            end %for d=1:length(dFFnames)
            
                        
            %% Individual trials- heatmaps
            orderedmap = 0; % 0 to plot in normal order, 1 if you want to plot in the order of deg inhibition by sorting on the min or max value reached within desired timewindow
            for d=1; %1:length(dFF_names)
                for r = 2; %1:length(raw_or_corr)
                    for i=2; %1:length(datatype)                        
                        for pow = 1:length(Opto_Powers)
                            figure;
                            if show_plot == 0
                                set(gcf,'visible','off')
                            end
                            alltrials_array = [];
                            indivdata = [IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})];
                            numrows = size(indivdata,1);
                            Color_scale = [0 10];

                            if isempty(alltrials_array)                              
                                alltrials_array(1:numrows,:) = indivdata;
                            else
                               numcurrentrows = size(alltrials_array,1);
                               alltrials_array(numcurrentrows+1:numcurrentrows+numrows,:) = indivdata;
                            end
                            % if you want each animal to have its trials sorted by order in terms of deg inhibition
                            if orderedmap == 1
                                [~,order] = sort(abs(max(alltrials_array(:,t_trials > 2 & t_trials < stim_duration),[],2))); %if A is a matrix, then min(A,[],2) is a column vector containing the minimum value of each row
                                %B = sort(A) sorts the elements of A in ascending order + [B,I] = sort(___) also returns a collection of index vectors
                            else
                                order = [1:numrows]';
                            end
                            %plot
                            imagesc(t_trials,1,alltrials_array(order,:))
                            colormap('jet')
                            c = colorbar;
                            c.Label.String = 'dF/F'; 
                            c.Label.FontSize = 16
                            if ~isempty(Color_scale)
                                c_limits = Color_scale;
                                caxis(c_limits)
                            end
                            hold on
                            xlim([-10 30])
                            xlabel('Time (s)')
%                             ylabel('                        DETQ    Veh')
                            set(gca,'Yticklabel',[]);
                            set(gca,'YTick',1:1:size(alltrials_array,1));
                            if pow == 1 % VEH
                                for q=1:size(alltrials_array,1)
                                    time2plot = (5 * (q+(VEH_trial1-1)))-5;
                                    text(-9.5,q,{time2plot},'Color','white','FontSize',12);
                                    text(-7.2,q,'min','Color','white','FontSize',12);
                                end
                            else
                                for q=1:size(alltrials_array,1)
                                    time2plot = (5 * (q+(DETQ_trial1-VEH_trialend-1)))-5;
                                    text(-9.5,q,{time2plot},'Color','white','FontSize',12);
                                    text(-7.2,q,'min','Color','white','FontSize',12);
                                end
                            end
                             set(gca,'FontSize',16,'FontName', 'Arial')
                            ax = gca;
                            set(gca,'TickLength',[0 0])
                            box off
                            ylimits = get(gca,'YLim');
                            plot([0 0],ylimits,'k','LineWidth',1)
                            plot([stim_duration stim_duration],ylimits,'k','LineWidth',1,'LineStyle',':')
                            sgtitle([raw_or_corr{r},' ',datatype{i},' ',Opto_Powers{pow}]);
                            if save_plot == 1
                                saveas(gcf,[PATH2SAVE,'figures\Heatmaps of dFF signals ',raw_or_corr{r},' ',datatype{i},' ',Opto_Powers{pow},'.tif'])
                                saveas(gcf,[PATH2SAVE,'figures\Heatmaps of dFF signals ',raw_or_corr{r},' ',datatype{i},'.fig'])
                            end
                        end
                    end
                end
            end
            
            
        %% NOW SEPARATE DATA VEH/DETQ INTO DIFFERENT MATRICES FOR LATER ANALYSES
        % SAVE INDIVDIDUAL SESSION    
        if done == 0 || overwrite == 1        
            if contains(path2data,'2022')
                % VEH 
                Opto_Powers = {'VEH'}; drug_name = Opto_Powers;
                % indiv stim data
                for d=1:length(dFF_names)
                    for pow = 1:length(Opto_Powers)
                        for j = 1:length(datatype)
                            IndivStim_data_VEH.raw.(dFF_names{d}).DRUG.(datatype{j}) = IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{j});
                            IndivStim_data_VEH.baselinecorr.(dFF_names{d}).DRUG.(datatype{j}) = IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{j});
                        end
                    end
                end
                % streams aligned
                for d=1:length(dFF_names)
                    for k=1:length(datatype)
                        for pow = 1:length(Opto_Powers)
                            streams_aligned_VEH.(datatype{k}).(dFF_names{d}).DRUG = streams_aligned.(datatype{k}).(dFF_names{d}).(Opto_Powers{pow});
                            time_vect_aligned_VEH.DRUG = time_vect_aligned.(Opto_Powers{pow});
                        end
                    end
                end
                % save
                save([PATH2SAVEWORKSPACE_VEH,'IndividualData.mat'],'IndivStim_data_VEH','streams','time_vect','streams_aligned_VEH','time_vect_aligned_VEH',...
                                                                'epocs','t_trials','datatype','dFF_names','Opto_Powers','drug_name','Events2Align2',...
                                                                'BASELINE_WIN','stim_duration','TRANGE','dt_ds','sampling_rate_ds','length_data','time_optotrials_double');
                % DETQ
                Opto_Powers = {'DETQ'}; drug_name = Opto_Powers;
                % indiv stim data
                for d=1:length(dFF_names)
                    for pow = 1:length(Opto_Powers)
                        for j = 1:length(datatype)
                            IndivStim_data_DETQ.raw.(dFF_names{d}).DRUG.(datatype{j}) = IndivStim_data.raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{j});
                            IndivStim_data_DETQ.baselinecorr.(dFF_names{d}).DRUG.(datatype{j}) = IndivStim_data.baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{j});
                        end
                    end
                end
                % streams aligned
                for d=1:length(dFF_names)
                    for k=1:length(datatype)
                        for pow = 1:length(Opto_Powers)
                            streams_aligned_DETQ.(datatype{k}).(dFF_names{d}).DRUG = streams_aligned.(datatype{k}).(dFF_names{d}).(Opto_Powers{pow});
                            time_vect_aligned_DETQ.DRUG = time_vect_aligned.(Opto_Powers{pow});
                        end
                    end
                end
                % save
                save([PATH2SAVEWORKSPACE_DETQ,'IndividualData.mat'],'IndivStim_data_DETQ','streams','time_vect','streams_aligned_DETQ','time_vect_aligned_DETQ',...
                                                                'epocs','t_trials','datatype','dFF_names','Opto_Powers','drug_name','Events2Align2',...
                                                                'BASELINE_WIN','stim_duration','TRANGE','dt_ds','sampling_rate_ds','length_data','time_optotrials_double');            

         end    

        %% Ending all the loops
        end % if done
    end %for s=1:length sessions
end % for all mice: i=1:numfiles, see above

    
