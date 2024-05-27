% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT
% Name: Shock_PooledGraphs
% First creation Sept 2020; Updates for Zurich May 2021, August 2021, March
% 2023, March 2024
% Marie Labouesse, marie.labouesse@gmail.com
% Cohort: GCAMP or dLight1.3b with shock

% 1- POOLED DATA
% loads FP data previously analyzed (with: SLA_Shock_DataExtraction) in matlab space "IndividualData.mat" from multiple animals (select folder containing multiple animals) 
% puts individual trials into one big matrix (all mice together): PooledINDIV
% calculates the mean of individual trial trajectories aligned to specific event for each mouse: PooledAVE
% saves a workspace that can then be used for quantifications (Shock_PooledQuantification)

% 2- POOLED GRAPHS
% generates average graphs curves (this used to be in the main script)
% generates heatmaps (1 average trace/mouse or all trials/all mice in one big heatmap - separated by power)

% edit relevant parameters in %% SETUP PARAMETERS

% functions needed: error_area, error_area_onlyrectangle

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
PATH2DATA_0 = uigetdir('select folder'); %select the overall folder, 1 level above the folder in which all days are (in each day there is multiple animals)
TrialDay = {'Day1','Day2','Day3_SCH'}; % DETQ
TrialDay = {'Day1','Day2'}; % VEH
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
            if  strcmp(mice_list_virus.(TrialDay{v})(o).name,'data') == 1 || strcmp(mice_list_virus.(TrialDay{v})(o).name,'figures') == 1 ...   %remove folders with data or figures or ., ..
                || contains(mice_list_virus.(TrialDay{v})(o).name,'data') || contains(mice_list_virus.(TrialDay{v})(o).name,'figures') ...
                || contains(mice_list_virus.(TrialDay{v})(o).name,'results')  || contains(mice_list_virus.(TrialDay{v})(o).name,'other') ...
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
        AnimalIDs(nummice) = {['ID',mice_list_virus.(TrialDay{v})(nummice).name(1:4)]};
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
                 || strcmp(sessions(o).name,'.') == 1 || strcmp(sessions(o).name,'..') == 1 || contains(sessions(o).name,'.csv') == 1 || contains(sessions(o).name,'Representative')
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
                load([PATH2SESSION,'\IndividualData.mat']);
                
        
                % Initialization of pooled structure with individual trials    %................ ONLY ONE SESSION PER MOUSE
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
                                % STREAMS                
                                fullstreams.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = {};
                                % SEM
                                PooledSEM.(TrialDay{v}).raw.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan;  
                                PooledSEM.(TrialDay{v}).baselinecorr.(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}) = ones(length(mice_list_virus.(TrialDay{v})),size(t_trials,2))*nan; 

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
                                %N/A
                            end
                        end
                    end
                end
                %STORE STREAMS
                fullstreams.(TrialDay{v}).(AnimalIDs{nummice}).(SessionIDs{s}) = streams;
              
            end
        end
    end
end

%% Save pooled data for all viruses together
% STORE INDIVIDUAL and AVERAGE TRIALS
save([PATH2SAVEPOOL,'\PooledAllMice.mat'],'PooledINDIV','PooledAVE','PooledSEM','AnimalIDs','mice_list_virus','TrialDay','PooledAnimalID','PooledAnimalID_only','dFF_names','datatype','Opto_Powers',...
    't_trials','stim_duration','TRANGE','BASELINE_WIN','Events2Align2','dt_ds','sampling_rate_ds','length_data','fullstreams','time_vect','pooledtype');





%% %% Individual trials- curves: ALL GROUPS TOGETHER ON SAME GRAPH
% Create plot with nanmean and nanstd, Save graph

if contains(PATH2SESSION,'DETQ')
color2plot = {'k','r','b','m'};
else
color2plot = {'k','m','b','m'};    
end

TrialType_Number = length(Opto_Powers);

% for k = 1:size(Body_parts,2)
for k = 3; %1:size(datatype,2) %dFF within: n=3
    for d=1:length(dFF_names)
        for p = 2; %:length(pooledtype) % baselinecorr only p=2
            figure; clf; hold on
            subplot(1,2,1)
            if show_plot == 0
                set(gcf,'visible','off')
            end

            for v=length(TrialDay):-1:1 
                for pow = 1:length(Opto_Powers)
                    tmp_avg_pool = nanmean(PooledAVE.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1); 
                    tmp_error_pool = nanstd(PooledAVE.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}),1,1)./...
                    sqrt(sum(~isnan(PooledAVE.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k})(:,1))));
                    error_area(t_trials,tmp_avg_pool,tmp_error_pool,color2plot{v},0.25); %t_trials is the time vector
                end
            end
            % plot properties

            % Formatting
            xlim([-10 20]);
            ylim([-2 6]); %see parameters setting at start                     
            ax = gca;
            ax.FontSize = 20; 
            xlabel('Time (s)','FontSize', 20)
            ylabel('\DeltaF/F (%)','FontSize', 20)
            box off
            set(gca,'TickLength',[0 0])
            box off
            ylimits = get(gca,'YLim');
            xlimits = get(gca,'XLim');
            set(gca,'TickLength',[0 0])
            %             xticks([-30 0 30 60 90]);
            %plot the shock area
            error_area_onlyrectangle([0 1],ylimits,[ylimits(2)-ylimits(1) ylimits(2)-ylimits(1)],[0.1000, 0.0950, 0.0980],0.1); %(X,Y,barlength,color,alpha for shading of SEM). Barlength should be Ylimit_max - Ylimitmin
            plot([0 0],ylimits,'k','LineWidth',0.5,'LineStyle','-.')
            yline(0,'-.k');
            
            if contains(PATH2SESSION,'DETQ')
            TrialDay_Drug = {'Veh','DETQ','SCH-23390'}
            else
            TrialDay_Drug = {'Veh','Veh','SCH-23390'}
            end
            % Legend, self-made - with text only (OPTION 2)
            for v=1:length(TrialDay)
                text4legend = TrialDay_Drug{v};
                % text4legend(end-1:end)=[];
%                 plot([xlimits(1)+1.5 xlimits(1)+3],[(ylimits(2)-0.1)-v*0.8 (ylimits(2)-0.1)-v*0.8],'color',color2plot{v},'LineWidth',2,'LineStyle','-')
                text([xlimits(1)+13 xlimits(1)+13],[(ylimits(2)+0.5)-v*0.6 (ylimits(2)+0.5)-v*0.6],text4legend,'Color',color2plot{v},'FontSize',18)
            end


            %saveplot or not
            if done == 0 && save_plot == 1 || done == 1 && reanalysis == 1 && overwrite == 1
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\AVE_alldays_plots',pooledtype{p},' ',datatype{k},'.tif']);
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\AVE_alldays_plots',pooledtype{p},' ',datatype{k},'.fig']);
                saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\AVE_alldays_plots',pooledtype{p},' ',datatype{k},'.pdf']);

            end
        end
    end
end
% end



%% Individual trials- heatmaps; one line per trial; all mice ; subplots 
if contains(PATH2SESSION,'detq')
color2plot = {'k','r','b','m'};
else
color2plot = {'k','m','b','m'};    
end

figure;

for k = 3; %size(datatype,2) % dFF within
figure;
if show_plot == 0
    set(gcf,'visible','off')
end
hold on
orderedmap = 0; % 0 to plot in normal order, 1 if you want to plot in the order of deg inhibition by sorting on the minimal value between 3 and 10 sec (for 10 sec stim...)
for v=1:length(TrialDay)
    for p = 2; %1:length(pooledtype)
        for d=1:length(dFF_names)
               
               

                subplot(1,length(TrialDay),v);
                for pow = 1 %:length(Opto_Powers)
                    alltrials_array = [];
                    for nummice=1:length(mice_list_virus.(TrialDay{v}))
                        animals = PooledAnimalID_only.(TrialDay{v});
                        for s = 1:length(sessions)
                            indivdata = PooledINDIV.(TrialDay{v}).(pooledtype{p}).(dFF_names{d}).(Opto_Powers{pow}).(datatype{k}).(animals{nummice}).(SessionIDs{s});
                            numrows = size(indivdata,1);
                            if isempty(alltrials_array)                               
                                alltrials_array(1:numrows,:) = indivdata;
                            else
                               numcurrentrows = size(alltrials_array,1);
                               alltrials_array(numcurrentrows+1:numcurrentrows+numrows,:) = indivdata;
                            end
                        end
                    end
                    if orderedmap == 0
                        numrows = size(alltrials_array,1);
                        order = [1:numrows]';
                    % if you want each animal to have its trials sorted by order in terms of deg increase
                    elseif orderedmap == 1
                        [~,order] = sort(abs(max(alltrials_array(:,t_trials > 0 & t_trials < stim_duration+1),[],2))); %if A is a matrix, then min(A,[],2) is a column vector containing the minimum value of each row
                        %B = sort(A) sorts the elements of A in ascending order + [B,I] = sort(___) also returns a collection of index vectors
                    % if you want to put animals in order in order of how good they are (but their trials stay in order)
                    elseif orderedmap == 2
                       A= abs(max(alltrials_array(:,t_trials > 0 & t_trials < stim_duration+3),[],2)); %if A is a matrix, then min(A,[],2) is a column vector containing the minimum value of each row
                        %B = sort(A) sorts the elements of A in ascending order + [B,I] = sort(___) also returns a collection of index vectors
                       numrows = size(alltrials_array,1);
                       numtrials = size(indivdata,1);
                       for i=1:numtrials
                           if i==1
                                mouse_score(i) = mean(A(1:numtrials))
                           elseif i==numtrials
                                mouse_score(i) = mean(A(numrows-numtrials+1:numrows))
                           else
                                mouse_score(i) = mean(A((i-1)*numtrials+1:(i-1)*numtrials+numtrials))
                           end
                       end
                       [~,order_overall] = sort(mouse_score)
                       order = ones(numrows,1)*nan;
                       for nummice=1:length(mice_list_virus.(TrialDay{v}))
                           order((nummice*numtrials)-numtrials+1) = order_overall(nummice)*numtrials-numtrials+1;
                           order((nummice*numtrials)-numtrials+2) = order((nummice*numtrials)-numtrials+1)+1;
                           order((nummice*numtrials)-numtrials+3) = order((nummice*numtrials)-numtrials+1)+2;
                       end

                    end
                    %plot
                    imagesc(t_trials,1,alltrials_array(order,:))
                    colormap(flipud(autumn))
                    colormap(parula)
                    % colormap((autumn))
                    Color_scale = [0 1.6];

                    if ~isempty(Color_scale)
                        c_limits = Color_scale;
                        caxis(c_limits)
                    end
                    xlim([-10 20]);
                    
                    set(gca,'YTick',3:6:size(alltrials_array,1));
%                                 sgtitle([datatype{k},' ',pooledtype{p}],'FontSize', 20);

                    %title
                    if contains(PATH2SESSION,'detq')
                    TrialDay_Drug = {'Veh','DETQ','SCH-23390'}
                    else
                    TrialDay_Drug = {'Veh','Veh','SCH-23390'}
                    end
                    text4legend = TrialDay_Drug{v};
                    box off
                    ylimits = get(gca,'YLim');
%                     
%                     % Formatting
                    ax = gca;
                    ax.FontSize = 30; 
                    title([TrialDay_Drug{v}],'FontSize', 30,'color',color2plot{v});
        
                    xlabel('Time to Shock (s)')
                    box off
                    ylimits = get(gca,'YLim');
                    xlimits = get(gca,'XLim');
                    set(gca,'TickLength',[0 0])
                    %plot the shock area
                    xline(0,':k','Linewidth',2);
                    xline(stim_duration,':k','Linewidth',2);

                    % on for getting the Y axis legends color bar, off
                    % otherwise
                    colorbar('Ticks',[0,3,6])
                    c = colorbar;
                    ylabel('Trials')
                    c.Label.String = '\DeltaF/F (%)';
                    c.Label.FontSize = 30;

                end

%                 clear alltrials_array

                
            end
        end
    end
end
fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1 1];  

                
% if save_plot == 1
%     saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\Heatmaps of dFF signals',pooledtype{p},' ',datatype{k},'.tif']);
%     saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\Heatmaps of dFF signals',pooledtype{p},' ',datatype{k},'.fig']);
%     saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\Heatmaps of dFF signals',pooledtype{p},' ',datatype{k},'.pdf']);
% end

if save_plot == 1
    saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\Heatmaps of dFF signals_withlegend',pooledtype{p},' ',datatype{k},'.tif']);
    saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\Heatmaps of dFF signals_withlegend',pooledtype{p},' ',datatype{k},'.fig']);
    saveas(gcf,[PATH2SAVEPOOL,'\pooled figures\Heatmaps of dFF signals_withlegend',pooledtype{p},' ',datatype{k},'.pdf']);
end




