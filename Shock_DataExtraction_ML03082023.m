% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT: 
% Name: Shock_DataExtraction
% First creation Sept 2020; Updates for Zurich May 2021, August 2021, March 2023
% Marie Labouesse, marie.labouesse@gmail.com
% Cohort: PFC dLight1.3b with shock

% extracts FP data (select folder from a particular day eg day 1, containing multiple animals inside)
% generates and saves a matrix with individual trials aligned to specific event (shock)
% generates and saves graphs aligned to specific event

% edit relevant parameters in %% SETUP PARAMETERS

% functions needed: error_area_onlyrectangle, TDTbin2mat


%% INITIALIZATIONS
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% SETUP PARAMETERS
%%%%%%%%%% TTL EXTRACTION
ttl_extract = 'default'; %'default' or 'special' (needs to be 7 characters)
          
%%%%%%%%%% DATA EXTRACTION
raw_or_corr = {'raw','baselinecorr'};
channel_number = 1; % 1 if 1-site recording, 2 if 2-site recording
channel_names = {'PFC'}; % eg brain region
color_number = 2; % 2 if 405 and 465 recording, 3 if also 565
color_names = {'c405','c465'}; % can also be 565
dFF_number = 1; % how many dFF to analyze, can be 1 or 2, based on number of channels and colors (if more than 2, need to edit script)
dFF_names = {'PFC'}; %put them in this order: 1) channel 1, 465 and 2) channel 1, 560 OR channel 2: 465 (if other combinations, need to edit script)
datatype = {'ZScoredFF','dFF','dFFwithin'}; %data analyses of interest -- remove speed if not analysing behavior
keeprawdata = 1; %1 if you want to keep the raw data to calculate dFF within trials. 0 otherwise ; % make sure to add dFFwithin to datatype

%%%%%%%%%%% PREPROCESSING
% trimming
timetrim_start_set = 180; %time to chop off at start of recording, in seconds: adjust this manually, eg 60 seconds
timetrim_end_set = 1; %time to chop off at end of recording, eg 1 second
% low pass filtering of raw 405 and 470
lowpassfiltering = 1; % 1 if you want to low pass filter- typically yes- 0 for no
lowpassfreq = 1; % typically 1Hz
% detrending dFF
detrending = 0;% %% Write 1 for detrend dFF, normal detrend function of Matlab- typically no (Write 0) for open field
% substract 8th percentile as a baseline calculate using a moving window
deletemovingwindow = 0; % 1 to delete moving baseline 8th percentile, 0 if no. Typically no 
% high pass filtering of dFF
highpassfiltering = 0; %'yes' if you want to high pass filter, 0 for no. Typically no for open field
highpassfreq = 0.005; % typically 0.005Hz

%%%%%%%%%%% TRIAL DEFINITION
% time window for the trials and graphs
TRANGE = [-90 90]; % will create events for a -5 to +5 sec window (you can always plot less later)
% baseline correction window
BASELINE_WIN.OptoStim = [-10 -1];% baseline correction window% BASELINE_PER = [-5 -1]; 
% variables to align to
Events2Align2 = {'PulseStart_time','PulseStop_time'}; %can be {'optostim_onset','optostim_offset'} for opto, or any other TTL epoc in the "epocs" structure with same organization

%%%%%%%%%%% EXPERIMENT TYPE
ExperimentType = 'other'; % If other, then the trials are not defined and need to go in and change trial allocation.
Opto_Powers = {'STIM'}; % VEH VS DETQ. Just write 'STIM' if not relevant
TrialType_Number = length(Opto_Powers);
stim_duration = 1; % 10 seconds or 3 seconds or 5 seconds
stim_duration2 = [];  % in case their is a ramp or another duration of interest

%%%%%%%%%%% HOW MANY MICE TO ANALYZE
LoopOrNot = 1; % 1 to loop, 0 if you only want to test one mouse and dont want to loop --> if so, edit the mouse number you want below: nummice_set
nummice_set = 6;

%%%%%%%%%%% SHOW, SAVE, OVERWRITE PARAMETERS
show_plot = 0; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
reanalysis = 1; % If 1, the code runs for sessions already analyzed. If 0, session excluded if analysed. This is based on the existence of the matlab space "IndividualData.mat"
overwrite = 1; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not. Code not implementing this.          


%% DEFINE PATHS
path2data = uigetdir('select folder'); % SELECT FOLDER OF GROUP TO ANALYZE,  - %Location of the data
mice_list = dir(path2data); %all things in this folder
% Define paths and data to analyze
path2savefolder = path2data; %Path to save


%% IDENTIFY MICE TO ANALYZE 
for o = length(mice_list):-1:1
    if mice_list(o).isdir == 0  %remove non-folders
        mice_list(o) = [];
    else
        if  strcmp(mice_list(o).name,'data') == 1 || strcmp(mice_list(o).name,'figures') == 1 ...   %remove folders with data or figures or ., ..
            || contains(mice_list(o).name,'data') || contains(mice_list(o).name,'figures') ...
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
for nummice = 1:Nmice   
    if LoopOrNot == 0
        nummice = nummice_set;  %need to set the number of the mouse you want to analyze up above in PARAMETERS
    end
    AnimalID = ['ID_',mice_list(nummice).name(end-3:end)];   %last 4 digits of folders should be the animal name
    Dirsplit = strsplit(path2data,'\');
    Virus_cell = Dirsplit(length(Dirsplit));
    Virus = Virus_cell{:}; %just extracting the string out of the cell array   -->  (name of the folder you selected)
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
        elseif strcmp(sessions(o).name,'data') == 1 || strcmp(sessions(o).name,'figures') == 1 || strcmp(sessions(o).name,'results') == 1 ...
             || strcmp(sessions(o).name,'.') == 1 || strcmp(sessions(o).name,'..') == 1 || contains(sessions(o).name,'.csv') || contains(sessions(o).name,'other') == 1
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
    for s = 1:length(sessions)
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
        % Check if results are already saved for this session 
        done = exist([PATH2SAVE,'IndividualData',AnimalID,'.mat'],'file'); %if you specify 'file' then matlab searches for both files and folders
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
            stream_name_type405A = {'x05A','05A','A05A','x405A'}; %possible names on rigs- can add to this list
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
                    sampling_rate = data.streams.(stream_name_type405A{h}).fs;
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
            dt_ds = 1/sampling_rate_ds; % delta t or unit of time (in seconds) per sample --> for 101.73 Hz sampling rate (after downsampling), should be close to 10 ms (9.8 ms)
            time_vect = 0:dt_ds*1:dt_ds*(max_idx-1); %time vector; values in seconds, same number of data points as the data


            %% Extract TTL data (epocs) for opto stim and camera frames
            % Initialize epocs structure
            if ttl_extract == 'default' 
                epoc_type = {'PulseStart_time','PulseStop_time','StimAmp','StimTime','TTLstart_time','TTLstop_time','Testsync_start_time','Testsync_stop_time','Framesync_start_time',...
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
                    epocs.PulseStart_time = data.epocs.U22_.onset;
                    epocs.PulseStop_time = data.epocs.U22_.offset;
                end
                % Opto Stim amplitude
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

            %% if TTL opto stim didnt work
            
            epocs.PulseStart_time = [];
            epocs.PulseStop_time = [];
            epocs.PulseStart_time = epocs.TTLstart_time;
            epocs.PulseStop_time = epocs.TTLstop_time;
     
            %% TRIMMING
            % Setting up trim indexes
            timetrim_start = timetrim_start_set*1; %in seconds: adjust this manually, eg trim 1.5 minutes: 60*1.5 --> set up at the beginning of the code
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

        
            %% LOW PASS FILTER OF FP DATA TO 1Hz
            if lowpassfiltering == 1
                %filter option 2: butter
                ftype = 'low';
                n = 2; % 2nd order filter, like Arturo and Thomas Akan
                Wn = lowpassfreq/((sampling_rate_ds)/2); %lowpassfreq defined above
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


            %% CALCULATE dFF (synapse and polyfit versions)
            % Initialize
            F0_405 = []; F0_465 = []; F0_560 = []; F0_405_2 = []; F0_465_2 = [];
            dFF_405 = []; dFF_465 = []; dFF_560 = []; dFF_405_2 = []; dFF_465_2 = [];
            % Calculations in each channel and color
            % First channel, 465 color, first dFF
            channel=1;

            % % dFF polyfit
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
                xlabel('time (sec)'); ylabel('dFF '); title([dFF_names{d},' Normalized data']);
                L=legend('dFF','dFF Synapse','Stim on','Stim off'); L.Location = 'Best';
                min_x=-30; max_x=max(time_vect)+30; min_y=min(streams.dFF.(dFF_names{d}))-3; max_y=max(streams.dFF.(dFF_names{d}))+3; axis([min_x max_x min_y max_y]); 
                if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                    saveas(gcf,[PATH2SAVE,'figures\dFF calculations ',dFF_names{d},'.tif'])
                    saveas(gcf,[PATH2SAVE,'figures\dFF calculations ',dFF_names{d},'.fig'])
                end
            end

            %% DETRENDING dFF
            if detrending == 1;% %% Detrend dFF, Option 1: normal detrend (other detrend options not here, but could include them later, perhaps as a function)
                for d=1:length(dFF_names)
                    streams.dFF_dtr.(dFF_names{d}) = detrend(streams.dFF.(dFF_names{d})); % THIS WORKS BUT DOESNT DO MUCH
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
                    % Define parameters for moving window % 
                    baseline_window = 10; %ADJUST! %Window to make the calculation from eg if (60 sec): you calculate the baseline for the 60 sec, looking at 50 sec ahead of the current 10sec 
                    baseline_dt = 5; %Moving steps for the threshold --> eg move 10 sec - value has to be < baseline_window
                    percentile = 8; %8th percentile
                    % New time vector
                    temp_t = time_vect - time_vect(1); %creating new time vector set to starting at 0, incase time_vect weas not starting at 0
                    dummie = 1:length(temp_t); %generates index of same amount of samples (better than find cos too much computing time)
                    dummie = dummie(temp_t >= baseline_window); %value when the time is beyond the window --> so that we can apply the moving window calculations
                    idx_baseline_window = dummie(1); %index of the timestamp as of which we can start the 60sec moving window
                    dummie = 1:length(temp_t);
                    dummie = dummie(temp_t >= baseline_dt);
                    idx_baseline_dt = dummie(1); %index of the timestamp as of which we start the 10sec moving window, i.e. baseline dt
                    clear dummie temp_t
                    % Index for the moving steps
                    ix1 = 1:idx_baseline_dt:length_data-idx_baseline_window; %take whole duration and then look at index you will look at (eg 1 to 60, 1 to 70 etc). calculate it before doing the loop, makes it easier
                    % create index vector for the first bound of the window, where steps are of the size of 10 sec scaled to the
                    ix2 = idx_baseline_window:idx_baseline_dt:length_data; %index for the second bound of the window
                    ix = [ix1' ix2']; %first and last index you will check in each window (first colum onset, second colum offset) and the rows are the windows
                    if ix(end,2) < length_data
                        [ix] = [ix;[length_data-idx_baseline_window length_data]]; %trick to include the last bit
                    end
                    % Calculate 8th percentile across Pre and Post datasets
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
                n = 2; % 2nd order filter, 
                Wn = highpassfreq/((sampling_rate_ds)/2); %highpassfreq: 0.005Hz ? defined above
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

             
            %% IF WE CALCULATE dFF WITHIN THEN ADD AN EMPTY ARRAY

            if keeprawdata == 1
                for d=1:length(dFF_names)    
                    streams.dFFwithin.(dFF_names{d}) = streams.ZScoredFF.(dFF_names{d})*nan
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
%             Events2Align2 = {'TTLstart_time','TTLstop_time'}     
            Events2Align2 = {'PulseStart_time','PulseStop_time'}                
            t_opto = [epocs.(Events2Align2{1}) epocs.(Events2Align2{1})+2]; % t_opto: a double column with all trials: column 1 (all rows) the start of the opto events and column 2 (all rows) the end
            for d=1:length(dFF_names)
                for o = 1:length(Opto_Powers)
                    t_opto2.(dFF_names{d}).(Opto_Powers{o}) = t_opto; % t_opto2: group the ones to average from t_opto 
                end
            end
            timeDETQ=30; % this is approximate, just for plotting on graphs

            end
            
            %% Generate the split dFF data (and other datatypes) trial matrixes based on trial allocation (t_opto2)
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
                            title([Opto_Powers{pow}(2:end)])
                        end
                        sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d}]); % for groups of subplots
                        if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                            saveas(gcf,[PATH2SAVE,'figures\INDIV TRIALS ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\INDIV TRIALS ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.fig'])
                        end
                    end
                end
            end %for d=1:length(dFFnames)

                               
           

            %% Average plot of all individual trials per type of trials with all types of trials on same graph
            for d=1:length(dFF_names)
                %RAW
                color2plot = {'b','g','m','r'};
                for i=1:length(datatype)
                    for r=1:length(raw_or_corr)
                        figure; 
                        if show_plot == 0
                            set(gcf,'visible','off')
                        end
                        hold on
                        for pow = 1:length(Opto_Powers)
                            tmp_avg = nanmean(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i}),1);
                            tmp_error = nanstd(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i}),1,1)./...
                                sqrt(sum(~isnan(IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i})(:,1))));
                            error_area(t_trials,tmp_avg,tmp_error,color2plot{pow},0.25); % error_area(X, Y, barlength, color, alpha, varargin) to plot SEM
                            max_avg(pow) = max(tmp_avg);

                        end
                        xline(0,'-k'); 
                        xline(stim_duration,'-m');
                        yline(0,'-.k');
                        xlabel('Time (s)');
                        ylabel(datatype{i});
                        xlim([TRANGE(1) TRANGE(2)]);
                        sgtitle([raw_or_corr{r},' ',datatype{i},' ',dFF_names{d}],'Interpreter','none')

%                         saveplot or not
                        if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                            saveas(gcf,[PATH2SAVE,'figures\AVE all powers optostim ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\AVE all powers optostim ',raw_or_corr{r},' ',datatype{i},' ',dFF_names{d},'.fig'])
                        end
                    end
                end
            end %for d=1:length(dFFnames)
            
              
                       
            %% Individual trials- heatmaps
            orderedmap = 0; % 0 to plot in normal order, 1 if you want to plot in the order of deg inhibition by sorting on the minimal value between 3 and 10 sec (for 10 sec stim...)
            for d=1:length(dFF_names)
                for r = 1:length(raw_or_corr)
                    for i=1:length(datatype)
                            Color_scale = [-1 3];
                        figure;
                        if show_plot == 0
                            set(gcf,'visible','off')
                        end
                        for pow = 1:length(Opto_Powers)
%                             subplot(length(Opto_Powers),1,pow);
                            alltrials_array = [];
                            indivdata = IndivStim_data.(raw_or_corr{r}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{i});
                            numrows = size(indivdata,1);
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
                            c.Label.String = 'Intensity (A.U.)';
                            if ~isempty(Color_scale)
                                c_limits = Color_scale;
                                caxis(c_limits)
                            end
                            hold on
                            xlim([t_trials(1) t_trials(end)])
            %                 xlim([t_trials(1) 10])
                            xlabel('Time (s)')
                            ylabel('Trial')
%                             set(gca,'YTick',[]);
                            set(gca,'Yticklabel',[]);
                            set(gca,'YLabel',[]);
                            set(gca,'YTick',1:1:size(alltrials_array,1));
                            for q=1:size(alltrials_array,1)
                                time2plot = round((t_opto2.(dFF_names{d}).(Opto_Powers{pow})(q,1)-timeDETQ)/60);
                                text(-13,q,{time2plot},'Color','black','FontSize',8);
                                text(-12,q,'min','Color','black','FontSize',8);
                            end
                            title([raw_or_corr{r},' ',datatype{i}]);
                            box off
                            ylimits = get(gca,'YLim');
                            plot([0 0],ylimits,'k','LineWidth',2)
                            plot([stim_duration stim_duration],ylimits,'k','LineWidth',2,'LineStyle',':')
                            clear alltrials_array
                        end
%                         sgtitle([dFF_names{d}]);
                        if save_plot == 1
                            saveas(gcf,[PATH2SAVE,'figures\Heatmaps of dFF signals ',raw_or_corr{r},' ',datatype{i},'.tif'])
                            saveas(gcf,[PATH2SAVE,'figures\Heatmaps of dFF signals ',raw_or_corr{r},' ',datatype{i},'.fig'])
                        end
                    end
                end
            end

        
        
        %% SAVE INDIVDIDUAL SESSION    

        if done == 0 || overwrite == 1
            save([PATH2SAVE,'IndividualData.mat'],'IndivStim_data','epocs','t_trials','datatype','dFF_names','Opto_Powers','Events2Align2','BASELINE_WIN',...
                'stim_duration','TRANGE','dt_ds','sampling_rate_ds','length_data','streams','time_vect','streams');
        end    


        %% Ending all the loops

    end %for s=1:length sessions
end % for all mice: i=1:numfiles, see above

    
