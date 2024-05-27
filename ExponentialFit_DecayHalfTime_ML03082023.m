function DecayHalfTime = ExponentialFit_DecayHalfTime(fit_timevect,object_2fit,expfit_type,StartPoint)

% DATA ANALYSIS PIPELINE FOR FIBER PHOTOMETRY EXPERIMENT
% Name: ExponentialFit_DecayHalfTime

% Feb 2023 - Marie Labouesse, marie.labouesse@gmail.com

% Exponential Fit a function to get its DecayHalfTime

% These are basic matlab function from here:
% https://ch.mathworks.com/help/curvefit/exponential.html

% TIPS:
% from the matlab function, If you specify start points, choose values appropriate to your data. Set arbitrary start points for coefficients a and b for example purposes. e.g.
    % fitobject= fit(xFit,yFit,'exp1','StartPoint',[0.4,-1]);
% Be careful if you have your data vectors starting with the baseline: you want only to fit the signal after the Laser stim, so introduce an appropriate range for your xfit and yfit
% Plot your fit on top of the data for the inspection. exp1 fit work well on the data which is going to the baseline at some point. 
% If the signal  never goes to baseline, the exponential fit may not be good
% It might well be that exp2 fits your decay data better, but then it is not that standardized to calculate the Decay half time; new code is needed
% Important to choose the start of your fit well(e.g. which range to choose for yfit vector). 
    % For the laser Stim, I think the easy and straightforward way would be to start after the Laser stimulation is over. 
    % In case the signal is still rising afterward, you can decide to start fit after the signal reaches its Maxima

% Double exponential decay model, from https://www.biorxiv.org/content/10.1101/2022.01.09.475513v1.full
% Then the fluorescence profile (dF/F) was normalized to the maximum intensity. 
% The data points following the maximum intensity peak were fit using a double exponential decay model.
% In the fitted curve, the time point where the normalized fluorescence profile passed under 36.8% of the maximum intensity was selected 
    % as the transient decay time (or lifetime) of the phasic dopamine activity in each brain regions.

%% if want to run within the function
% fit_timevect = 0:dt_ds:dt_ds*(length(t_opto_i_postpeak)-1);
% fit_timevect = fit_timevect';
% % object_2fit = smooth(t_opto_i_postpeak,100000);
% object_2fit = t_opto_i_postpeak';
% expfit_type = 'exp1';
% StartPoint = [5,0];
    
%% function
    
if expfit_type == 'exp1'
    %% exponential fit 1
    % yFit -would be your signal to fit, eg dFF after the Laser Stim
    % xFit the respective time vector

    fitobject= fit(fit_timevect,object_2fit,'exp1');
    coefficients = coeffvalues(fitobject);
    FittedCurve = coefficients(1)*exp(coefficients(2)*fit_timevect);
    figure; plot(fit_timevect,object_2fit); hold on; plot(fit_timevect,FittedCurve); L=legend('Plot','Fit'); L.Location = 'Best';
    DecayHalfTime=log(2)/abs(coefficients (2));
    % plot


else
    %% exponential fit 2
    fitobject= fit(fit_timevect,object_2fit,'exp2');
    coefficients = coeffvalues(fitobject);
    FittedCurve = coefficients(1)*exp(coefficients(2)*fit_timevect)+coefficients(3)*exp(coefficients(4)*fit_timevect); 
    figure; plot(fit_timevect,object_2fit); hold on; plot(fit_timevect,FittedCurve); L=legend('Plot','Fit'); L.Location = 'Best';
    DecayHalfTime_ix=find(FittedCurve<(36.8*max(FittedCurve)/100));
    DecayHalfTime=fit_timevect(DecayHalfTime_ix(1));
    % plot

end




