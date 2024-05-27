% ChemoSLA_PooledGraphs_2023
% Marie Labouesse, marie.labouesse@gmail.com - Sep 2023
% Cohort: NAc dLight1.3b weekly DETQ or food rest

% 1- POOLED DATA
% load matlab spaces generated in: WeeklyDETQ_SLA_DataExtraction
% puts individual trials into one big matrix (all mice together): PooledINDIV
% calculates the mean of individual trial trajectories aligned to specific event for each mouse: PooledAVE
% saves a workspace that can then be used for quantifications

% 2- POOLED GRAPHS
% generates average graphs curves (this used to be in the main script)

% NEED FUNCTIONS
%error_area_only_rectangle

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
PATH2DATA_0 = uigetdir('select folder'); %select the overall folder, where the different groups are saved (as folders)- we will import all
virus = {'wk1','wk2','wk3'}; % weekly
% virus = {'FED','FOODRESTR','FED2'}; 


for v=1:length(virus)
    PATH2DATA.(virus{v}) = [PATH2DATA_0,'\',virus{v}];
    PATH2SAVEFOLDER.(virus{v}) = [PATH2DATA_0,'\',virus{v}];
    mice_list_virus.(virus{v}) = dir(PATH2DATA.(virus{v})); %all things in this folder
    mkdir([PATH2SAVEFOLDER.(virus{v}),'\pooled figures\']);
    mkdir([PATH2SAVEFOLDER.(virus{v}),'\pooled data\']);
    mkdir([PATH2DATA_0,'\pooled figures\']);
    mkdir([PATH2DATA_0,'\pooled data\']);
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
        AnimalIDs(nummice) = {['ID',mice_list_virus.(virus{v})(nummice).name(1:4)]};
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
        
        %% Loop for all the sessions for the mouse: Sessions are experimental days you want to average together (replicates of each other)- we only have 1 session here
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
       
                % Initialization of pooled structure with individual trials   %................  ONLY ONE SESSION PER MOUSE.

                if nummice == 1 && s == 1
                    for d=1:length(dFF_names)
                        for pow = 1:length(Opto_Powers)
                            for k = 1:size(datatype,2)
                                % ANIMAL ID .................... add session ID later if needed
                                PooledAnimalID.(virus{v}) = {};
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
                
                for d=1:length(dFF_names)
                        for pow = 1:length(Opto_Powers)
                            for k = 1:size(datatype,2)
                                % Add data one animal at a time
                                % ANIMAL ID and SESSION ID
                                PooledAnimalID.(virus{v}) = AnimalIDs;

                                %STORE STREAMS
                                fullstreams.(virus{v}).(AnimalIDs{nummice}) = streams.dFF;

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

%% Adjust all to VEH 0 as baseline
for v=1:length(virus)    
    for nummice=1:length(mice_list_virus.(virus{v}))
        for d=1 %:length(dFF_names)
            for k = 1 %:size(datatype,2)
                for pow = 1:length(Opto_Powers)
                    fullstreams_aligned_norm.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) = ...
                        fullstreams_aligned.(virus{v}).(Opto_Powers{pow}).(AnimalIDs{nummice}) ...
                        - mean(fullstreams_aligned.(virus{v}).VEH.(AnimalIDs{nummice}))                    
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
                time_vect_aligned_ds_minutes.(Opto_Powers{pow}) = time_vect_aligned_ds.(Opto_Powers{pow})/60-10;
                
            end
        end
    end
end



    
%% Plot 
color2plot = {'k','r','b','c','g','y','m'}; % food rest
% color2plot = {'r','b','k','c','g','y','m'}; % weeklz DETQ
figure; clf; hold on
if show_plot == 0
    set(gcf,'visible','off')
end

% AVERAGE GRAPH
idx_selct = [1,2,3,4,5,6,7]; % week1,2,3
% idx_selct = [1,2,3,4,5]; %fed vs food-restricted


for pow = 1:length(Opto_Powers)
    for v=1:length(virus)
        % average plot
        tmp_avg_pool = nanmean(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,:),1); 
        tmp_error_pool = nanstd(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,:),1,1)./...
        sqrt(sum(~isnan(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,1))));
        error_area(time_vect_aligned_ds_minutes.(Opto_Powers{pow}),tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector

    end
end

% Formatting
ylim([-10 60]); % food rest
% ylim([-20 80]); % weeklz DETQ

xlim([-10 50]);
xticks([-10 0 10 20 30 40 50]);
ax = gca;
ax.FontSize = 20; 
xlabel('Time (min)','FontSize', 20)
ylabel('\DeltaF/F (%)','FontSize', 20)
box off
set(gca,'TickLength',[0 0])

xline(-10,':k'); xline(-9,':k','Veh','FontSize',20); % grey out the noise in signal due to injection 
xline(0,':k'); xline(1,':k','DETQ','FontSize',20); % grey out the noise in signal due to injection 
% xline(20,':k'); xline(21,':k','J60 or Sal','FontSize',20); 
error_area_onlyrectangle([-10 -9],[40 120],[150 150],[165, 167, 168]./255,0.8); 
error_area_onlyrectangle([0 1],[40 120],[150 150],[165, 167, 168]./255,0.8); 
% error_area_onlyrectangle([20 21],[40 120],[150 150],[165, 167, 168]./255,0.8); 

%WEEKLY DETQ
yline(27.11,'-k'); yline(27.11*0.85,'--k'); yline(27.11*1.15,'--k','Plateau','FontSize',20,'FontAngle','italic'); % plateau
xline(5.238,'--k'); xline(25.1,'--k','FontSize',20); % temporal window
annotation('textbox', [0.348, 0.31, 0.2, 0.05], 'String', 'Temporal window', 'Color', 'black', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none','FontAngle','italic');

% week 1,2,3
annotation('textbox', [0.718, 0.90, 0.2, 0.05], 'String', 'Week 1', 'Color', 'red', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none');
annotation('textbox', [0.718, 0.85, 0.2, 0.05], 'String', 'Week 2', 'Color', 'blue', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none');
annotation('textbox', [0.718, 0.8, 0.2, 0.05], 'String', 'Week 3', 'Color', 'black', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none');

% % FED VS FOOD REST
% annotation('textbox', [0.60, 0.90, 0.2, 0.05], 'String', 'Ad lib fed', 'Color', 'black', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none','FitBoxToText','on');
% annotation('textbox', [0.60, 0.85, 0.2, 0.05], 'String', 'Food restriction', 'Color', 'red', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none','FitBoxToText','on');
% annotation('textbox', [0.60, 0.80, 0.2, 0.05], 'String', 'Ad lib fed', 'Color', 'blue', 'FontSize', 18, 'EdgeColor', 'none','Interpreter', 'none','FitBoxToText','on');

%saveplot or not
if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
    saveas(gcf,[PATH2DATA_0,'\pooled figures\streams_aligned_overlap ',' ',datatype{k},'.tif']);
    saveas(gcf,[PATH2DATA_0,'\pooled figures\streams_aligned_overlap ',' ',datatype{k},'.fig']);
    saveas(gcf,[PATH2DATA_0,'\pooled figures\streams_aligned_overlap ',' ',datatype{k},'.pdf']);
end          

%% SUBPLOT INDIDIVUAL GRAPH
color2plot = {'r','b','k','c','g','y','m'};

idx_selct = [1,2,3,4,5,6,7]; % week1,2,3
% idx_selct = [1,2,3,4,5]; %fed vs. food-restricted

figure; clf; hold on
if show_plot == 0
    set(gcf,'visible','off')
end

figure;

% Get screen size
screenSize = get(0, 'ScreenSize');

% Set the figure position to fill the entire screen
set(gcf, 'Position', [1, 1, screenSize(3)/2-20, screenSize(4)-80]);

for mouse_number = 1:length(idx_selct)
    subplot(4,2,mouse_number)
    for pow = 1:length(Opto_Powers)
        for v=length(virus):-1:1
            % average plot
            tmp_avg_pool = nanmean(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct(mouse_number),:),1); 
            tmp_error_pool = nanstd(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct(mouse_number),:),1,1)./...
            sqrt(sum(~isnan(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct(mouse_number),1))));
            error_area(time_vect_aligned_ds_minutes.(Opto_Powers{pow}),tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector
        end
    end

% Formatting
ylim([-20 80]);  
xlim([-10 50]);
xticks([-10 0 10 20 30 40 50]);
ax = gca;
ax.FontSize = 16; 
xlabel('Time (min)','FontSize', 16)
ylabel('\DeltaF/F','FontSize', 16)
box off
set(gca,'TickLength',[0 0])

xline(-10,':k'); xline(-9,':k','V.','FontSize',16); % grey out the noise in signal due to injection 
xline(0,':k'); xl = xline(1,':k','D.','FontSize',16); % grey out the noise in signal due to injection 
% xl.LabelHorizontalAlignment = 'left';

% xline(20,':k'); xline(21,':k','J60 or Sal','FontSize',20); 
error_area_onlyrectangle([-10 -9],[40 120],[150 150],[165, 167, 168]./255,0.8); 
error_area_onlyrectangle([0 1],[40 120],[150 150],[165, 167, 168]./255,0.8); 
% error_area_onlyrectangle([20 21],[40 120],[150 150],[165, 167, 168]./255,0.8); 

mouse_number_string = string(mouse_number);
% Modify the title line to set the title centered above the subplot
t = title(['mouse ' char(mouse_number_string),':'], 'FontSize', 16, 'FontWeight', 'bold');

% Adjust the position of the title to be centered above the subplot
set(t, 'Position', [31, 68, 0]);

end

% Weekly DETQ
% Add text at the far bottom right of the figure with specified colors
annotation('textbox', [0.56, 0.25-0.015, 0.2, 0.05], 'String', 'Week 1', 'Color', 'red', 'FontSize', 16, 'EdgeColor', 'none');
annotation('textbox', [0.56, 0.225-0.015, 0.2, 0.05], 'String', 'Week 2', 'Color', 'blue', 'FontSize', 16, 'EdgeColor', 'none');
annotation('textbox', [0.56, 0.200-0.015, 0.2, 0.05], 'String', 'Week 3', 'Color', 'black', 'FontSize', 16, 'EdgeColor', 'none');

% Food restriction
% annotation('textbox', [0.56, 0.25-0.015, 0.2, 0.05], 'String', 'Fed', 'Color', 'red', 'FontSize', 16, 'EdgeColor', 'none');
% annotation('textbox', [0.56, 0.225-0.015, 0.2, 0.05], 'String', 'Food-restr', 'Color', 'blue', 'FontSize', 16, 'EdgeColor', 'none');
% annotation('textbox', [0.56, 0.200-0.015, 0.2, 0.05], 'String', 'Fed 2', 'Color', 'black', 'FontSize', 16, 'EdgeColor', 'none');

%saveplot or not
if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
    saveas(gcf,[PATH2DATA_0,'\pooled figures\indiv subplots ',' ',datatype{k},'.tif']);
    saveas(gcf,[PATH2DATA_0,'\pooled figures\indiv subplots ',' ',datatype{k},'.fig']);
    saveas(gcf,[PATH2DATA_0,'\pooled figures\indiv subplots ',' ',datatype{k},'.pdf']);

    % print([PATH2DATA_0,'\pooled figures\indiv subplots ',' ',datatype{k},'.pdf'], '-dpdf', '-fillpage');

end 


%% Quantifications

% idx_selct = [1,2,3,4,5,6,7]; %weekly DETQ
idx_selct = [1,2,3,4,5]; %fed vs. restricted
wanttempwindow = 0; %0 if fed vs restricted, 1 if weekly dETQ

% getting the indexes
Section={'Drug1','Drug2'}
temp_t = time_vect_aligned_ds.VEH - time_vect_aligned_ds.VEH(1); 
dummie = 1:length(temp_t); 
dummie2 = dummie(temp_t >= 599); % 600 sec = 10 min --> first 10 min
idx_Drug1_2 = dummie2(1);  
idx_Drug1_1 = find(temp_t == 0);
idx_Drug1 = idx_Drug1_1:idx_Drug1_2; 

temp_t = time_vect_aligned_ds.DETQ - time_vect_aligned_ds.DETQ(1); 
dummie = 1:length(temp_t); 
dummie3 = dummie(temp_t >= 60*49-1);  % 49 min = end for weekly DETQ and for the plateau calc for food-restr; % 15 min for calculating the AUC food restriction
idx_Drug2_1 = dummie3(1); 
% dummie4 = find(temp_t >= 60); % 1min into DETQ
dummie4 = find(temp_t >= 0); % from time 0, in this case it's 1 min into DETQ (coded earlier in the extraction code)
idx_Drug2_2 = dummie4(1); 
idx_Drug2 = idx_Drug2_2:idx_Drug2_1-1; 

idx_drug_all = {'idx_Drug1','idx_Drug2'};
Day = {'veh','detq'};

% initialize Measurements
for du=1:length(Day)
    for v=1:length(virus)    
        for nummice2=1:length(idx_selct); %size(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}),1)
            nummice = idx_selct(nummice2)

            Measurements_AVE.AUC.(Day{du}) = ones(length(virus),length(mice_list_virus.(virus{v})))*nan; %done like this, for easier copy-paste in GraphPad. day1 (rows=weeks (virus) and columns=animals)
            Measurements_AVE.Baseline.(Day{du}) = ones(length(virus),length(mice_list_virus.(virus{v})))*nan;
            Measurements_AVE.TempWindowStart.(Day{du}) = ones(length(virus),length(mice_list_virus.(virus{v})))*nan;
            Measurements_AVE.TempWindowEnd.(Day{du}) = ones(length(virus),length(mice_list_virus.(virus{v})))*nan;
            Measurements_AVE.TempWindowDur.(Day{du}) = ones(length(virus),length(mice_list_virus.(virus{v})))*nan;
            Measurements_AVE.PlateauVal.(Day{du}) = ones(length(virus),length(mice_list_virus.(virus{v})))*nan;     

        end                
    end
end

for du=1:length(Day) %here Day means drug (veh or detq)
    for v=1:length(virus) 
        for nummice2=1:length(idx_selct); %size(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}),1)
            nummice = idx_selct(nummice2)
            
            if du == 1; idx_select = idx_Drug1;  pow = 1; elseif du == 2; idx_select = idx_Drug2; pow = 2; end
            % find the right dff using the indexes
            dff_select = fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(nummice,idx_select);

            % each row is a new week (here variable virus)
            % and each column is a new mouse
            Measurements_AVE.AUC.(Day{du})(v,nummice) = trapz(dff_select)/length(idx_select);
            Measurements_AVE.Baseline.(Day{du})(v,nummice) = nanmean(dff_select,2);

            if wanttempwindow == 1
                if du ==2
                % calculate the temporal window
                % create vector with 1-min bins 
                dff_select_1min_bins_smoothed = smooth(dff_select',300);  
                ttt = timetable(dff_select_1min_bins_smoothed,'SampleRate',1);%1fps
                tttAvg = retime(ttt,'regular','mean','SampleRate',0.016666666666666666666666666666666666666666666666666666); %0.167fps= every minute 
                dff_select_1min_bins = tttAvg.dff_select_1min_bins_smoothed;
                time_vect_1min_bins = [minutes(tttAvg.Time)];
    
                % Determine the maximal peak maxima, average and average value at the plateau (5-24min). This is to calculate the plateau. (written a bit redundant below, could simplify the code)
                timepoint5min = find(time_vect_1min_bins == 5); timepoint6min = find(time_vect_1min_bins == 6);timepoint7min = find(time_vect_1min_bins == 7);timepoint8min = find(time_vect_1min_bins == 8);timepoint9min = find(time_vect_1min_bins == 9);
                timepoint10min = find(time_vect_1min_bins == 10);timepoint11min = find(time_vect_1min_bins == 11);timepoint12min = find(time_vect_1min_bins == 12);timepoint13min = find(time_vect_1min_bins == 13);timepoint14min = find(time_vect_1min_bins == 14);
                timepoint15min = find(time_vect_1min_bins == 15);timepoint16min = find(time_vect_1min_bins == 16);timepoint17min = find(time_vect_1min_bins == 17);timepoint18min = find(time_vect_1min_bins == 18);timepoint19min = find(time_vect_1min_bins == 19);
                timepoint20min = find(time_vect_1min_bins == 20);timepoint21min = find(time_vect_1min_bins == 21);timepoint22min = find(time_vect_1min_bins == 22);timepoint23min = find(time_vect_1min_bins == 23);timepoint24min = find(time_vect_1min_bins == 24);
                EMPIRIMAX_AbsMaxima_all = mean([dff_select_1min_bins(timepoint5min),dff_select_1min_bins(timepoint6min),dff_select_1min_bins(timepoint7min),dff_select_1min_bins(timepoint8min),dff_select_1min_bins(timepoint9min)...
                                                dff_select_1min_bins(timepoint10min),dff_select_1min_bins(timepoint11min),dff_select_1min_bins(timepoint12min),dff_select_1min_bins(timepoint13min),dff_select_1min_bins(timepoint14min)...
                                                dff_select_1min_bins(timepoint15min),dff_select_1min_bins(timepoint16min),dff_select_1min_bins(timepoint17min),dff_select_1min_bins(timepoint18min),dff_select_1min_bins(timepoint19min)...
                                                dff_select_1min_bins(timepoint20min),dff_select_1min_bins(timepoint21min),dff_select_1min_bins(timepoint22min),dff_select_1min_bins(timepoint23min),dff_select_1min_bins(timepoint24min)]);
          
                % Determine the trials within 15% of the plateau
                WINSTABL_AbsMaxima_all = ones(1,length(dff_select_1min_bins))*nan;
                for ju = 1:size(WINSTABL_AbsMaxima_all,2);
                    if dff_select_1min_bins(ju) > 0.84 * EMPIRIMAX_AbsMaxima_all &&  dff_select_1min_bins(ju) < 1.16 * EMPIRIMAX_AbsMaxima_all
                        WINSTABL_AbsMaxima_all(ju) = 1;
                    else
                        WINSTABL_AbsMaxima_all(ju) = 0;
                    end
                end
                WINSTABL_AbsMaxima_all
    
                % Calculate variables
                % Find indexes (written a bit redundant below, could simplify the code)
                % first ix, its going to be either 5 or 10 min (not before 5 min: the kinetics are too fast, it would not be stable enough or reliable)
                if WINSTABL_AbsMaxima_all(6) == 1;  WINSTABL_IX(1) = 6; % 5 min timepoint               
                    elseif WINSTABL_AbsMaxima_all(7) == 1; WINSTABL_IX(1) = 7; % 6 min timepoint
                    elseif WINSTABL_AbsMaxima_all(8) == 1; WINSTABL_IX(1) = 8; % 7 min timepoint
                    elseif WINSTABL_AbsMaxima_all(9) == 1; WINSTABL_IX(1) = 9; % 8 min timepoint
                    elseif WINSTABL_AbsMaxima_all(10) == 1; WINSTABL_IX(1) =10; % 9 min timepoint
                    elseif WINSTABL_AbsMaxima_all(11) == 1; WINSTABL_IX(1) = 11; % 10 min timepoint                
                end
                % second ix, its going to be anyting as of 15 min        
                if  WINSTABL_AbsMaxima_all(16) == 0; WINSTABL_IX(2) = 16; % 15 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(17) == 0; WINSTABL_IX(2) = 17; % 16 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(18) == 0; WINSTABL_IX(2) = 18; % 17 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(19) == 0; WINSTABL_IX(2) = 19; % 18 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(20) == 0; WINSTABL_IX(2) = 20; % 19 min timepoint            
                    elseif  WINSTABL_AbsMaxima_all(21) == 0; WINSTABL_IX(2) = 21; % 20 min timepoint            
                    elseif  WINSTABL_AbsMaxima_all(22) == 0; WINSTABL_IX(2) = 22; % 21 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(23) == 0; WINSTABL_IX(2) = 23; % 22 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(24) == 0; WINSTABL_IX(2) = 24; % 23 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(25) == 0; WINSTABL_IX(2) = 25; % 24 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(26) == 0; WINSTABL_IX(2) = 26; % 25 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(27) == 0; WINSTABL_IX(2) = 27; % 26 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(28) == 0; WINSTABL_IX(2) = 28; % 27 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(29) == 0; WINSTABL_IX(2) = 29; % 28 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(30) == 0; WINSTABL_IX(2) = 30; % 29 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(31) == 0; WINSTABL_IX(2) = 31; % 30 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(32) == 0; WINSTABL_IX(2) = 32; % 31 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(33) == 0; WINSTABL_IX(2) = 33; % 32 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(34) == 0; WINSTABL_IX(2) = 34; % 33 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(35) == 0; WINSTABL_IX(2) = 35; % 34 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(36) == 0; WINSTABL_IX(2) = 36; % 35 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(37) == 0; WINSTABL_IX(2) = 37; % 36 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(38) == 0; WINSTABL_IX(2) = 38; % 37 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(39) == 0; WINSTABL_IX(2) = 39; % 38 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(40) == 0; WINSTABL_IX(2) = 40; % 39 min timepoint
                    elseif  WINSTABL_AbsMaxima_all(41) == 0; WINSTABL_IX(2) = 41; % 40 min timepoint               
                end                
                WINSTABL_IX(2) = WINSTABL_IX(2) - 1;
    
                % Set variables
                TempWindowStart = (WINSTABL_IX(1)-1); % in minutes, e.g. if it's 6, then it's the 5 min timepoint
                TempWindowEnd = (WINSTABL_IX(2)-1); % in minutes, e.g. if it's 10, then it's the 9 min timepoint
                TempWindowDur = TempWindowEnd-TempWindowStart; % in minutes, e.g. if it's 5 min to 15 min, then it's 10 min
                PlateauVal = EMPIRIMAX_AbsMaxima_all;
    
                % Fill in Measurements
                Measurements_AVE.TempWindowStart.(Day{du})(v,nummice) = TempWindowStart;
                Measurements_AVE.TempWindowEnd.(Day{du})(v,nummice) = TempWindowEnd;
                Measurements_AVE.TempWindowDur.(Day{du})(v,nummice) = TempWindowDur;
                Measurements_AVE.PlateauVal.(Day{du})(v,nummice) = PlateauVal;     
                    
                else
                end
            end
        end                
    end
end  
% end

%% Get table for Graphpad with data in columns in bins of 1 min to calculate window of superiority in graphpad- one new table for each week

% Downsample the streams aligned (smoothed) to 1 frame per minute
for v=1:length(virus)    
    for d=1:length(dFF_names)
        for pow = 1:length(Opto_Powers)
            for k = 1:size(datatype,2)
                for nummice=1:size(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}),1)
                
                    data2plot=fullstreams_aligned_pooled_smooth.(virus{v}).(Opto_Powers{pow})(nummice,:);
                    ttt = timetable(data2plot','SampleRate',sampling_rate_ds);
                    tttAvg = retime(ttt,'regular','mean','SampleRate',0.0167); %1frame per minute
                    fullstreams_aligned_pooled_ds_1frpmin.(virus{v}).(Opto_Powers{pow})(nummice,:) = tttAvg.Var1;
                end
                
                time_vect_aligned_ds_1frpmin.(Opto_Powers{pow}) = [seconds(tttAvg.Time)] ; %also add the first timepoint from the previous timevector allowing to shift in time                
                % time_vect_aligned_ds_1frpmin.(Opto_Powers{pow}) = time_vect_aligned_ds_1frpmin.(Opto_Powers{pow}) - time_vect_aligned_ds.(Opto_Powers{pow})(1)
                time_vect_aligned_ds_1frpmin_minutes.(Opto_Powers{pow}) = round(time_vect_aligned_ds_1frpmin.(Opto_Powers{pow})/60);
                
            end
        end
    end
end


%% Save pooled data for all viruses together
% STORE INDIVIDUAL and AVERAGE TRIALS
save([PATH2DATA_0,'\PooledAllMice.mat'],'virus','PooledAnimalID','dFF_names','datatype','Opto_Powers','fullstreams_aligned_pooled_smooth','Measurements_AVE',...
    'dt_ds','sampling_rate_ds','length_data','fullstreams','time_vect','pooledtype','fullstreams_aligned','timevect_aligned','fullstreams_aligned_pooled_ds','time_vect_aligned_ds','time_vect_aligned_ds_minutes');


    
    
  
        