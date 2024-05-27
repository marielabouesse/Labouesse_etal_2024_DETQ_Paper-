% ChemoSLA_PooledGraphs_2023
% Marie Labouesse, marie.labouesse@gmail.com - Sep 2023
% Cohort: NAc dLight1.3b with chemogenetic inhibition and temporal window

% 1- POOLED DATA
% load matlab spaces generated in: ChemoSLA_DataExtraction
% puts individual trials into one big matrix (all mice together): PooledINDIV
% calculates the mean of individual trial trajectories aligned to specific event for each mouse: PooledAVE
% saves a workspace that can then be used for quantifications

% 2- POOLED GRAPHS
% generates average graphs curves (this used to be in the main script)

% NEED FUNCTIONS
%error_area_only_rectangle, error_area

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
PATH2DATA_0 = uigetdir('select folder'); %select the overall folder, where the diff groups are saved
virus = {'Day1','Day2','Day3','Day4','Day5','Day6'}; 


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
        AnimalIDs(nummice) = {['ID',mice_list_virus.(virus{v})(nummice).name(11:14)]};
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
        
        %% Loop for all the sessions for the mouse: Sessions are experimental days you want to average together (replicates of each other) ........only 1 session
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

        %STORE STREAMS ALIGNED TO  ONSET #1 with the average of each mouse in each row (eg 3 mice --> 3 rows)
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

    % Make new time vectors and also delete the NaNs in the VEH
%     for v=1:length(virus)
%         for d=1:length(dFF_names)
%             for pow = 1:length(Opto_Powers)
%                 for k = 1:size(datatype,2)
%                      length_data_aligned.(Opto_Powers{pow}) = size(find(~isnan(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow})(1,:))),2);
%                      time_vect_aligned_special.(Opto_Powers{pow}) = 0:dt_ds*1:dt_ds*(length_data_aligned.(Opto_Powers{pow})-1);
%                      fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow})(:,length_data_aligned.(Opto_Powers{pow})+1:end) = [] ;                   
%                 end
%             end
%         end
%     end

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



    
%% Plots
% ALL
color2plot = {'k','c','k','c','b','b'};

% order of the drugs: veh,veh,veh; veh,veh,j60-low; veh,detq,veh;
% veh,detq,j60-low; veh,detq,j60-hi; veh,veh,j60-high

% DETQ ONLY
virus_DETQ = [3,4,5];

figure; clf; hold on
if show_plot == 0
    set(gcf,'visible','off')
end

%show all or some of the mice
idx_selct = [1,2,3,4,5,6];

for pow = 1:length(Opto_Powers)
    for jui=1:length(virus_DETQ)
        v=virus_DETQ(jui)
        % average plot
        tmp_avg_pool = nanmean(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,:),1); 
        tmp_error_pool = nanstd(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,:),1,1)./...
        sqrt(sum(~isnan(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,1))));
        error_area(time_vect_aligned_ds_minutes.(Opto_Powers{pow}),tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector

    end
end

% Formatting
ylim([-20 60]);                    
ax = gca;
ax.FontSize = 20; 
xlabel('Time (min)','FontSize', 20)
ylabel('\DeltaF/F (%)','FontSize', 20)
box off
set(gca,'TickLength',[0 0])

xline(0,':k'); xline(1,':k',' ','FontSize',20); % grey out the noise in signal due to injection 
xline(10,':k'); xline(11,':k','  ','FontSize',20); % grey out the noise in signal due to injection 
xline(20,':k'); xline(21,':k',' ','FontSize',20); 
error_area_onlyrectangle([0 1],[40 120],[150 150],[165, 167, 168]./255,0.8); 
error_area_onlyrectangle([10 11],[40 120],[150 150],[165, 167, 168]./255,0.8); 
error_area_onlyrectangle([20 21],[40 120],[150 150],[165, 167, 168]./255,0.8); 

annotation('textbox', [0.18, 0.90, 0.2, 0.05], 'String', 'Veh', 'Color', 'black', 'FontSize', 20, 'EdgeColor', 'none');
annotation('textbox', [0.33, 0.90, 0.2, 0.05], 'String', 'DETQ', 'Color', 'black', 'FontSize', 20, 'EdgeColor', 'none');

annotation('textbox', [0.51, 0.90, 0.2, 0.05], 'String', 'Veh', 'Color', 'black', 'FontSize', 20, 'EdgeColor', 'none');
annotation('textbox', [0.51, 0.84, 0.2, 0.05], 'String', 'J60 low', 'Color', 'cyan', 'FontSize', 20, 'EdgeColor', 'none');
annotation('textbox', [0.51, 0.78, 0.4, 0.05], 'String', 'J60 high', 'Color', 'blue', 'FontSize', 20, 'EdgeColor', 'none');


% L=legend('Veh,Sal',' ','Veh,J60',' ','DETQ,Sal',' ','DETQ,J60-lo',' ','DETQ,J60-hi',' ','Veh,J60-hi','FontSize',10); L.Location = 'NorthEast';     
% L=legend('Veh,Sal',' ','Veh,J60',' ','DETQ,Sal',' ','DETQ,J60-lo',' ','DETQ,J60-hi','FontSize',10); L.Location = 'NorthEast';     
% L=legend('DETQ,Sal',' ','DETQ,J60-hi','FontSize',10); L.Location = 'NorthEast';     

%saveplot or not
if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
    saveas(gcf,[PATH2DATA_0,'\pooled figures\streams_aligned_overlap DETQ ',' ',datatype{k},'.tif']);
    saveas(gcf,[PATH2DATA_0,'\pooled figures\streams_aligned_overlap DETQ ',' ',datatype{k},'.fig']);
    saveas(gcf,[PATH2DATA_0,'\pooled figures\streams_aligned_overlap DETQ ',' ',datatype{k},'.pdf']);
end          

%% VEH ONLY
virus_VEH = [1,2,6];

figure; clf; hold on
if show_plot == 0
    set(gcf,'visible','off')
end

% show all or some mice
idx_selct = [1,2,3,4,5,6];

for pow = 1:length(Opto_Powers)
    for jui=1:length(virus_VEH)
        v=virus_VEH(jui)
        % average plot
        tmp_avg_pool = nanmean(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,:),1); 
        tmp_error_pool = nanstd(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,:),1,1)./...
        sqrt(sum(~isnan(fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(idx_selct,1))));
        error_area(time_vect_aligned_ds_minutes.(Opto_Powers{pow}),tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector

    end
end

% Formatting
ylim([-20 60]);                    
ax = gca;
ax.FontSize = 20; 
xlabel('Time (min)','FontSize', 20)
ylabel('\DeltaF/F (%)','FontSize', 20)
box off
set(gca,'TickLength',[0 0])

xline(0,':k'); xline(1,':k',' ','FontSize',20); % grey out the noise in signal due to injection 
xline(10,':k'); xline(11,':k','  ','FontSize',20); % grey out the noise in signal due to injection 
xline(20,':k'); xline(21,':k',' ','FontSize',20); 
error_area_onlyrectangle([0 1],[40 120],[150 150],[165, 167, 168]./255,0.8); 
error_area_onlyrectangle([10 11],[40 120],[150 150],[165, 167, 168]./255,0.8); 
error_area_onlyrectangle([20 21],[40 120],[150 150],[165, 167, 168]./255,0.8); 

annotation('textbox', [0.18, 0.90, 0.2, 0.05], 'String', 'Veh', 'Color', 'black', 'FontSize', 20, 'EdgeColor', 'none');
annotation('textbox', [0.35, 0.90, 0.2, 0.05], 'String', 'Veh', 'Color', 'black', 'FontSize', 20, 'EdgeColor', 'none');

annotation('textbox', [0.51, 0.90, 0.2, 0.05], 'String', 'Veh', 'Color', 'black', 'FontSize', 20, 'EdgeColor', 'none');
annotation('textbox', [0.51, 0.84, 0.2, 0.05], 'String', 'J60 low', 'Color', 'cyan', 'FontSize', 20, 'EdgeColor', 'none');
annotation('textbox', [0.51, 0.78, 0.4, 0.05], 'String', 'J60 high', 'Color', 'blue', 'FontSize', 20, 'EdgeColor', 'none');

%saveplot or not
if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
    saveas(gcf,[PATH2DATA_0,'\pooled figures\streams_aligned_overlap VEH ',' ',datatype{k},'.tif']);
    saveas(gcf,[PATH2DATA_0,'\pooled figures\streams_aligned_overlap VEH ',' ',datatype{k},'.fig']);
    saveas(gcf,[PATH2DATA_0,'\pooled figures\streams_aligned_overlap VEH ',' ',datatype{k},'.pdf']);
end 


%% Quantifications

% show all or some of the mice
idx_selct = [1,2,3,4,5,6];

Section={'Drug1','Drug2','Drug3'}
temp_t = time_vect_aligned_ds.VEH - time_vect_aligned_ds.VEH(1); 
dummie = 1:length(temp_t); 
dummie2 = dummie(temp_t >= 599); % until 10 min --> first 10 min (but 10 min were cut off for simplicity in the other code)
idx_Drug1_2 = dummie2(1);  
idx_Drug1_1 = find(temp_t == 0); % starting from beginning, so 10 min in total
idx_Drug1 = idx_Drug1_1:idx_Drug1_2; 

temp_t = time_vect_aligned_ds.DETQ - time_vect_aligned_ds.DETQ(1); 
dummie = 1:length(temp_t); 
dummie3 = dummie(temp_t >= 589);  % end = 10 min post DETQ, stop 20 sec before
idx_Drug2_1 = dummie3(1); 
dummie4 = find(temp_t >= 0); % % start at 10 minutes before the 10 min, so 10 min in total
idx_Drug2_2 = dummie4(1); 
idx_Drug2 = idx_Drug2_2:idx_Drug2_1-1; 

dummie5 = dummie(temp_t >= 600 + 120 + 600);  % end = 12 min post J60 (or 22 min post DETQ)
idx_Drug3_1 = dummie5(1); 
dummie6 = find(temp_t >= 600 + 120); % start at 12 min post DETQ, so 10 min in total
idx_Drug3_2 = dummie6(1); 
idx_Drug3 = idx_Drug3_2:idx_Drug3_1-1; 
clear dummie dummie2 dummie3 dummie4 dummie5 dummie6

idx_drug_all = {'idx_Drug1','idx_Drug2','idx_Drug3'};
Day = {'veh','detq'};

for du=1:length(Day)
for v=1:length(virus)    
    for d=1:length(dFF_names)
        for sect = 1:length(Section)
            for k = 1:size(datatype,2)
                for nummice2=1:length(idx_selct); %size(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}),1)
                    nummice = idx_selct(nummice2)         

                    Measurements_AVE.AUC.(Day{du}) = ones(length(Section)*3+2,length(mice_list_virus.(virus{v})))*nan;
                    Measurements_AVE.Baseline.(Day{du}) = ones(length(Section)*3+2,length(mice_list_virus.(virus{v})))*nan;

                end
                
            end
        end
    end
end
end

% for du=1:length(Day)
for v=1:length(virus) 
    if v == 1 || v == 2 || v == 6; du = 1; else du = 2; end;
    for d=1:length(dFF_names)
        for sect = 1:length(Section)
            for k = 1:size(datatype,2)
                for nummice2=1:length(idx_selct); %size(fullstreams_aligned_pooled.(virus{v}).(Opto_Powers{pow}),1)
                    nummice = idx_selct(nummice2)
                    
                    if sect == 1; idx_select = idx_Drug1; elseif sect == 2; idx_select = idx_Drug2; elseif sect == 3; idx_select = idx_Drug3; end
                    if sect == 1; pow = 1; else pow = 2; end;
                    dff_select = fullstreams_aligned_pooled_ds.(virus{v}).(Opto_Powers{pow})(nummice,idx_select);
                    
                    % if v==1, row_num=1; elseif v==2, row_num=2; elseif v==6, row_num=3; elseif v==3, row_num=5, elseif v==4, row_num=6; elseif v==5, row_num=7; end;
                    % order: veh,veh,veh (day 1); veh,veh,j60-low (day 2); veh,veh,j60-hi (day 6); 
                    % veh,detq,veh (day 3); veh,detq,j60-low (day 4); veh,detq,j60-hi (day 5); 
                    
                        if v == 1 % veh,veh,veh
                            if sect == 1, row_num = 1;  elseif sect == 2, row_num = 2;  elseif sect == 3, row_num = 3; end
                        elseif v == 2 % veh,veh,j60-low 
                            if sect == 1, row_num = 5;  elseif sect == 2, row_num = 6;  elseif sect == 3, row_num = 7; end                            
                        elseif v == 6 % veh,veh,j60-hi
                             if sect == 1, row_num = 9;  elseif sect == 2, row_num = 10;  elseif sect == 3, row_num = 11; end
                        elseif v == 3 % veh,detq,veh
                            if sect == 1, row_num = 1;  elseif sect == 2, row_num = 2;  elseif sect == 3, row_num = 3; end
                        elseif v == 4 % veh,detq,j60-low 
                            if sect == 1, row_num = 5;  elseif sect == 2, row_num = 6;  elseif sect == 3, row_num = 7; end                            
                        elseif v == 5 % veh,detq,j60-hi
                            if sect == 1, row_num = 9;  elseif sect == 2, row_num = 10;  elseif sect == 3, row_num = 11; end
                        end
                    

                    % each row is a new section (different drugs; same day)
                    % each chunk in rows is a new day (different drugs injected)
                    % and each column is a new mouse
                    Measurements_AVE.AUC.(Day{du})(row_num,nummice) = trapz(dff_select)/length(idx_select);
                    Measurements_AVE.Baseline.(Day{du})(row_num,nummice) = nanmean(dff_select,2);
                  
                end
                
            end
        end
    end
end  
% end



%% Save pooled data for all viruses together
% STORE INDIVIDUAL and AVERAGE TRIALS
save([PATH2DATA_0,'\PooledAllMice.mat'],'virus','PooledAnimalID','dFF_names','datatype','Opto_Powers','fullstreams_aligned_pooled_smooth','Measurements_AVE',...
    'dt_ds','sampling_rate_ds','length_data','fullstreams','time_vect','pooledtype','fullstreams_aligned','timevect_aligned','fullstreams_aligned_pooled_ds','time_vect_aligned_ds','time_vect_aligned_ds_minutes');


    
  
        