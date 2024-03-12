function plTab = plotDataMaker(resRel,varargin)
% Making plotting data: Finds mean best threshold for each paramSet/method
% and only retains those results in the plotting table
% Inputs: resRel: table of relative results
%         paramIdx = logical index of which rows to include
% Output: plTab: table of only those data to be plotted

if nargin==1 % Angle is default variable to optimize for
    optiVari = 'angleMu';
end
grsum = groupsummary(resRel,{'method','paramSet','thresh'},'mean'); % mean performance over reps
grsum.mean_agg = mean(abs(1-grsum{:,['mean_' optiVari]}),2); % Performance of the threshold aggregated across metrics (but not Std)
optiTab = groupsummary(grsum,{'method','paramSet'},'min'); % Best thresholds
optiIdx = ismember(grsum.mean_agg,optiTab.min_mean_agg); % Idx in grsum of best threshold
optiGrsum = grsum(optiIdx,:);
[~,uni] = unique(optiGrsum(:,1:2)); % Removes equally well performing threshs
optiGrsum = optiGrsum(uni,:);
paramOptiIdx = ismember(resRel(:,{'method' 'paramSet' 'thresh'}),...
                        optiGrsum(:,{'method' 'paramSet' 'thresh'}));
plTab = resRel(paramOptiIdx,:); % table containing only those to be plotted
end