%% Make parameters
% Popp, Codling, Bailey, Dornhaus, Feb 2024
% Requirements: 
%   R2023b or later (table operations & folder selection)
%   Signal Processing Toolbox for DTW (metric)
%   Navigation Toolbox for angdiff (resampling: heading angle to turn angle, could be hand-coded?)

% Glossary/naming conventions
%{
re: resampled
pl: plot
idx: logical index
ind: numerical index
nr...: Number of. Sorry for the German...

param: Ground truth track meta parameters (from which params are made)
params: Parameters from which ground truth tracks are made

For preprocessing empirical tracks for resampling & analysis, run
empiricalMaker.m (in funs folder)
%}

clear

data = "new"; % "new" "load" "empirical" "empiricalNew"
disp("Data: " + data)

% Path stuff
if ispc; slash = '\'; else; slash = '/'; end % OS compatibility
dirName = 'resample'; % Name of the base folder
addpath(genpath(dirName))
if ~license('test','Signal_Toolbox'); error('Requires Signal Processing Toolbox for DTW function'); end
if ~license('test','Navigation_Toolbox'); error('Requires Navigation Toolbox for angdiff function'); end

% Plot appearance
ftsz = {'fontsize',14};
lw = {'linewidth',1};
ps = 2.5; % Point size
clr = lines(14); % Color palette: lines (standard MATLAB discrete coloring)
clrGr = jet; % Color palette: jet (gradient coloring), for last figure
clrGr = clrGr(1:round(256/10):256,:); % Set /10 to number of methods

if data == "new"
    % Make parameters
    param.nrReps = 100; % Nr of simulations per parameter set
    param.length = [1000]; % Length of track (in 'noise steps'), first is default
    
    % Turns
    param.stepDistr = 1; % 1: "norm"; 2: "poiss"; 3: "exp" % Step length distribution
    param.stepMu = 10; % Mean (multiple of noise step-length)
    param.stepSD = 3; % SD for gaussian
    param.angleMu = [50]; % Turn angles
    param.angleSD = [1]; % Std for turn angles: Must !=0 (otherwise Inf results)
    
    % Noise & Smoothing
    param.noiseAngle = [0 5 10 15]; % SD of noise angle [°] at each noise step (0 for smooth)
    param.noiseXY = 0; % SD of x & y displacements from path [noise-steps]
    param.smoothW = [0]; % Smoothing window of LOWESS [noise steps], (def: 0:20)
    
    % Choose Resampling Methods. "corridor" & "cpt" Run super slow, thus not in tests
    methList = ["interval" "local" "nonlocal" "1D" "1D+" "corridor" "DP" "TPA" "CPT" "MRPA" "VW"];
    [methInd,~] = listdlg('PromptString','Select Resampling Methods', ...
                          'ListSize',[150 180], ...
                          'InitialValue',[1 2 3 4 6 7 8 9 10 11],...
                          'ListString',methList);
    nrThreshs = 5; % Number of resampling thresholds tested

    % === End of user input === %
    
    
    % Making parameter matrix
    [g.a, g.b, g.c, g.d, g.e, g.f, g.g, g.h, g.i] = ndgrid(param.length,param.stepDistr,param.stepMu,param.stepSD,param.angleMu,param.angleSD,param.noiseAngle,param.noiseXY,param.smoothW);
    paramMat = [g.a(:), g.b(:), g.c(:), g.d(:), g.e(:), g.f(:), g.g(:), g.h(:), g.i(:)];
    
    paramNames = [
        "length"...
        "stepDistr"...
        "stepMu"...
        "stepSD"...
        "angleMu"...
        "angleSD"...
        "noiseAngle"...
        "noiseXY"...
        "smoothW"];
    params = array2table(paramMat, 'variablenames',paramNames);

    % Resampling
    param.nrThreshs = nrThreshs;
    param.method = methList(methInd);
    param.nrMeths = numel(param.method);
    threshMat = [repmat("interval",nrThreshs,1) round(linspace(10,30,nrThreshs))' zeros(nrThreshs,1)    1*ones(nrThreshs,1);
                 repmat("local",nrThreshs,1)    ones(nrThreshs,1).*2         linspace(10,82,nrThreshs)' 2*ones(nrThreshs,1);
                 repmat("nonlocal",nrThreshs,1) zeros(nrThreshs,1)           linspace(20,76,nrThreshs)' 2*ones(nrThreshs,1);
                 repmat("1D",nrThreshs,1)       linspace(.1,20.1,nrThreshs)' zeros(nrThreshs,1)         1*ones(nrThreshs,1);
                 repmat("1D+",nrThreshs,1)      linspace(1,10,nrThreshs)'    zeros(nrThreshs,1)         1*ones(nrThreshs,1);
                 repmat("corridor",nrThreshs,1) linspace(.1,2.1,nrThreshs)'  zeros(nrThreshs,1)         1*ones(nrThreshs,1);
                 repmat("DP",nrThreshs,1)       linspace(.2,2.2,nrThreshs)'  zeros(nrThreshs,1)         1*ones(nrThreshs,1);
%                 repmat("vDP",nrThreshs,1)      linspace(.1,1,nrThreshs)'    zeros(nrThreshs,1)         1*ones(nrThreshs,1);
                 repmat("TPA",nrThreshs,1)   ones(nrThreshs,1).*param.stepMu linspace(1,41,nrThreshs)'  2*ones(nrThreshs,1);
                 repmat("CPT",nrThreshs,1)      ones(nrThreshs,1).*.05       round(linspace(15,60,nrThreshs))' 2*ones(nrThreshs,1);
                 repmat("MRPA",nrThreshs,1)     linspace(1,5,nrThreshs)'     zeros(nrThreshs,1)         1*ones(nrThreshs,1);
                 repmat("VW",nrThreshs,1)       linspace(10,30,nrThreshs)'   zeros(nrThreshs,1)         1*ones(nrThreshs,1);
                 ];
    threshsAll = array2table(threshMat(:,1:3),'VariableNames',{'method' 'thresh1' 'thresh2'});
    threshsAll.thresh1 = double(threshsAll.thresh1);
    threshsAll.thresh2 = double(threshsAll.thresh2);
    threshsAll.thresh = double(threshsAll.thresh1); % The tested threshold (1)
    threshsAll.thresh(threshMat(:,4)=="2") = threshMat(threshMat(:,4)=="2", 3); % (2)
    threshsAll.threshNr = repmat((1:nrThreshs)',height(threshsAll)/nrThreshs,1);
    threshs = threshsAll(ismember(threshsAll.method, param.method),:); % picks the user-selected methods

    % Nice threshold table format for the paper
    threshTab = table('size',[param.nrMeths,nrThreshs+1],'VariableTypes',["string", repmat("double",1,nrThreshs)]);
    threshTab.Var1 = unique(threshs.method);
    for meth = 1:param.nrMeths
        threshTab{meth,2:end} = threshs.thresh(threshs.method==threshTab{meth,1})';
    end

    % Dictionary between variable names and axis labels 
    % methNames = num2cell(unique(threshs.method)); % Meth names can stay cryptic?
    labelKeys = [fieldnames(param)]; % ; cellfun(@char, methNames, 'un',0)];
    labelChars = {'Number of replications' 'Track length',...
                  'Step Distribution' 'mean(step length)' 'SD(step length)',...
                  'mean(cos(turn angle))' 'SD(cos(turn angle))' 'SD(noise angle)' 'Noise displacement',...
                  'Smooth window' ['Resampling' newline ' method'],...
                  'Number of methods' 'Number of Thresholds'}';
    labelDict = dictionary(labelKeys, labelChars);
    param.labelDict = labelDict;
    
    
    % Save simulation & resampling params
    dataFolderName = [dirName slash char(datetime("today"))]; % Folder w/ date as name
    if ~isfolder(dataFolderName)
        mkdir(dataFolderName)
    else % Already made a folder today
        dataFolderName = [dirName slash char(datetime('now','format','yyyy-MM-dd_HHmm'))];
        mkdir(dataFolderName)
    end
    save([dataFolderName slash 'param'],"param") % Grd trth meta parameters (from which params are made)
    save([dataFolderName slash 'params'],"params") % Parameters from which grd trth tracks are made
    save([dataFolderName slash 'threshs'],"threshs") % Method names & their thresholds included
    mkdir([dataFolderName slash 'plots']) % Where figures are saved
    addpath(genpath(dataFolderName))
elseif data == "load"
    dataFolderName = [dirName slash '2023-11-29_bigNoLength']; % [Change to input dlg w/ 2023]
    load([dataFolderName slash 'param'])
    load([dataFolderName slash 'params'])
    load([dataFolderName slash 'threshs'])
    disp(['Params of ' dataFolderName ' loaded.'])
elseif data == "empirical" || data == "empiricalNew"
    dataFolderName = [dirName slash 'empirical'];
    load([dataFolderName slash 'param'])
    load([dataFolderName slash 'threshs'])
end

clearvars threshMat threshsAll nrThreshs g paramNames paramMat methInd labelKeys labelChars % [Ideally, delete all but 3 param files?]

%% Make tracks from parameters
% Cell array of size params x reps

if data  == "new"
    tracksWide = trackMaker(param,params);
    save([dataFolderName slash 'tracksGTWide'],"tracksWide")
    disp('Ground truth tracks made')
elseif data  == "load"
    load([dataFolderName slash 'tracksGTWide'],"tracksWide")
    disp('Ground truth tracks loaded')
elseif data == "empirical" || data == "empiricalNew"
    load([dataFolderName slash 'tracks']);
    disp('Empirical tracks loaded')
end

%% Resampling
% reTracks is a table w/ 4 columns for metadata + 1 col for the resampled
% tracks (which are a cell array)
% Grouped by param>method>thresh

if data == "new"
    nrParams = height(params);
    % Make meta data columns of resampled tracks
    tracksReWide = [table(repelem((1:nrParams)',height(threshs)),'variablenames',{'paramSet'}) repmat(threshs,nrParams,1)];
    
    % Make resampled tracks
    reNrMax = height(tracksReWide);
    wBar = waitbar(0,'Resampling');
    for reNr = 1:reNrMax
        for rep = 1:width(tracksWide)
            waitbar(reNr/reNrMax,wBar,['tracksRe set ' num2str(reNr) ' of ' num2str(reNrMax) newline,...
                                       'param set ' num2str(tracksReWide.paramSet(reNr)) ', method ' char(tracksReWide.method(reNr))])
            tracksReWide.track{reNr,rep} = resampler(tracksWide{tracksReWide.paramSet(reNr),rep}, ... % Track of this paramSet & rep
                                                     tracksReWide.thresh1(reNr), tracksReWide.thresh2(reNr), tracksReWide.method{reNr});
            % Adding s & alpha columns
            xyDist = diff(tracksReWide.track{reNr,rep}{:,{'x' 'y'}},1); % For s calc
            tracksReWide.track{reNr,rep}.s = [NaN; sqrt(sum(xyDist.^2, 2))];
            theta = [NaN; atan2(xyDist(:,1),xyDist(:,2))]; % angle in space [rad]
            tracksReWide.track{reNr,rep}.alpha = rad2deg([angdiff(theta); NaN]); % turn angle [deg]
        end
    end
    close(wBar)
    save([dataFolderName slash 'tracksReWide' num2str(reNrMax)],"tracksReWide")
    disp('Resampled tracks made')
elseif data == "load"
    load([dataFolderName slash 'tracksReWide'],"tracksReWide")
    disp('Resampled tracks loaded')
elseif data == "empiricalNew"
    tracksRe = [table(repelem(param.animal,height(threshs),1),'variablenames',{'animal'}) repmat(threshs,height(param.animal),1)];
    reNrMax = height(tracksRe);
    wBar = waitbar(0,'Resampling');
    for reNr = 1:reNrMax
        waitbar(reNr/reNrMax,wBar,['trackRe ' num2str(reNr) ' of ' num2str(reNrMax) ', animal ' char(tracksRe.animal(reNr))])
        tracksRe.track{reNr,1} = resampler(tracks{param.animal==tracksRe.animal(reNr)},tracksRe.thresh1(reNr),tracksRe.thresh2(reNr),tracksRe.method{reNr});
        
        % Adding s & alpha columns
        xyDist = diff(tracksRe.track{reNr}{:,{'x' 'y'}},1); % For s calc
        tracksRe.track{reNr,1}.s = [NaN; sqrt(sum(xyDist.^2, 2))];
        theta = [NaN; atan2d(xyDist(:,1),xyDist(:,2))]; % angle in space
        tracksRe.track{reNr,1}.alpha = [angdiff(deg2rad(theta)); NaN];
    end
    close(wBar)
    save([dataFolderName slash 'tracksRe'],'tracksRe')
elseif data == "empirical"
    load([dataFolderName slash 'tracksRe'],"tracksRe")
    disp('Empirical tracks loaded')
end

clearvars xyDist theta nrParams

%% Analysis
% resGt: results ground truth tracks
% resRe: results resampled tracks
% resRel: resampled results relative to grd truth (Re/Gt; not for curvy)
if data == "new"
    % Reshaping wide (params*meths*threshs x reps) results tracks tables to tall format
    tracksGt = repelem(tracksReWide(:,{'paramSet' 'method' 'thresh' 'threshNr'}),param.nrReps,1);
    tracksRe = tracksGt;
    rep = repmat((1:param.nrReps)',height(tracksGt)/param.nrReps,1);
    tracksGt.rep = rep;
    tracksRe.rep = rep;
    tracks = repelem(tracksWide,param.nrMeths*param.nrThreshs,1); % Imitating reTracksWide table
    tracksGt.track = reshape(tracks',[],1);
    tracksRe.track = reshape(tracksReWide.track',[],1);
    clearvars rep tracks

    disp('|-   | Making resRe & resGT')
    % ======================================================================= %
    % Metrics of ground truth (Gt) and resampled (Re) tracks (absolute)
    paramsCos = param;
    resGt = [repelem([table((1:height(params))','variablenames',{'paramSet'}) params],...
             param.nrMeths*param.nrThreshs*param.nrReps,1)...
             tracksGt(:,{'rep' 'method' 'thresh' 'threshNr'})];
    
    resRe = tracksRe(:,{'paramSet' 'rep' 'method' 'thresh' 'threshNr'});
    
    % #of turns
    resGt.nrTurns = cellfun(@(n) sum(n.turn), tracksGt.track);
    resRe.nrTurns = cellfun(@(n) height(n)-2, tracksRe.track); % -2 to account for 1st & last point
    
    % Track Length
    resGt.length = cellfun(@(n) sum(sqrt(diff([0; n.x]).^2 + diff([0; n.y]).^2)), tracksGt.track); % Smooth are mostly not 1000...
    resGt.length(repelem(params.smoothW==0,param.nrReps,1)) = 1000; % ...but in smooth, noiseXY is erroneously also counted
    resRe.length = cellfun(@(n) sum(sqrt(diff([0; n.x]).^2 + diff([0; n.y]).^2)), tracksRe.track);
    
    % Mean turn angles
    resRe.angleMu = cellfun(@(n) mean(abs(n.alpha(2:end-1))), tracksRe.track); % Excluding 1st & end point
    resRe.angleSD = cellfun(@(n) std(abs(n.alpha(2:end-1))), tracksRe.track);
    
    % Step-lengths (mean & SD)
    resRe.stepMu = resRe.length./resRe.nrTurns;
    resRe.stepSD = cellfun(@(n) std(n.s(2:end-1)), tracksRe.track);
    
    % Reshaping & adding params for ease of use
    paramsRep = repelem(params,height(threshs)*param.nrReps,1); % For easy calling (corresponds to rows in all res)
    paramsRep.Properties.VariableNames = cellfun(@(n) ['param_' n], paramsRep.Properties.VariableNames, 'un',0);
    resRe = [paramsRep resRe]; % For curvy track analysis & plotting
    
    disp('|--  | Making resRel')
    % ======================================================================= %
    % Relative values for step-and-turn
    resRel = [paramsRep resRe(:,{'paramSet' 'rep' 'method' 'thresh' 'threshNr'})];
    
    resRel.nrTurns = resRe.nrTurns./resGt.nrTurns;
    resRel.length = resRe.length./resGt.length;
    resRel.angleMu = resRe.angleMu./resGt.angleMu;
    resRel.angleSD = resRe.angleSD./resGt.angleSD;
    resRel.stepMu = resRe.stepMu./resGt.stepMu;
    resRel.stepSD = resRe.stepSD./resGt.stepSD;
    % Smooth doesn't have ground truth
    % resRel(resRel.param_smoothW>1,14:19) = resRe(resRel.param_smoothW>1,14:19); % WARNING: hard coded, would change if column numbers change
    resRel(resRel.param_smoothW>1,15:20) = resRe(resRel.param_smoothW>1,15:20); % WARNING: hard coded, would change if column numbers change
    
    % Curve similarity measures
    resRel.dtw = zeros(size(resRe.paramSet)); % Dynamic Time Warping distance
    paramSetList = find(params.smoothW==0); % List of indices of step&turn
    for paramSet = 1:numel(paramSetList) % Current step&turn parameter set
        for rep = 1:param.nrReps
            grdAll = tracksWide{paramSetList(paramSet),rep};
            grd = grdAll(grdAll.turn == 1,:);
            for method = 1:param.nrMeths
                for threshNr = 1:param.nrThreshs
                    reIdx = tracksRe.paramSet == paramSet &...
                            tracksRe.rep == rep &...
                            tracksRe.method == param.method(method) &...
                            tracksRe.threshNr == threshNr;
                    re = tracksRe.track{reIdx};
                    resRel.dtw(reIdx,1) = dtw(grd.t,re.t); % Difference in time of turns
                end
            end
        end
    end
    
    
    % Adding metrics to the basic param meta data
    param.metric = convertCharsToStrings(resRel.Properties.VariableNames(15:end));
    param.labelDict({'nrTurns';    'dtw'}) =...
                    {'No of turns';'DTW distance (lower=better)'}; % 'CDTW (lower=better)' may be added

    disp('|--- | Making bestRel & resCI')
    % ======================================================================= %
    % Getting best threshold, median & 90% interval for each param set (for plotting)
    % Table of only those resampled tracks with the best thresholds
    bestRel = plotDataMaker(resRel); % Finding optimal thresholds of step&turn GT tracks
    
    resCI5 = grpstats(resRe,["paramSet" "method" "thresh"], @(resRe) prctile(resRe,5), "DataVars",15:20);
    resCI95 = grpstats(resRe,["paramSet" "method" "thresh"], @(resRe) prctile(resRe,95), "DataVars",15:20);
    resCIrange = resCI95(:,5:10) - resCI5(:,5:10);
    % If I need the upper and lower bounds separately:
    % resCI5.Properties.VariableNames = strrep(resCI5.Properties.VariableNames,'Fun1','5prc');
    % resCI95.Properties.VariableNames = strrep(resCI95.Properties.VariableNames,'Fun1','95prc');
    % resCI = [resRe(1:param.nrReps:end,[1:8 13]) resCI5(:,[1:3 5:end]) resCI95(:,5:end)];
    resCI = [resRe(1:param.nrReps:end,[1:10 12:14]) resCIrange];
    resCI.Properties.VariableNames = strrep(resCI.Properties.VariableNames,'Fun1','rnge95'); % Replace in colnames 'Fun1' w/ rnge95
    
    disp('|----| Saving alysis data')
    % Saving output
    save([dataFolderName slash 'resRe'],"resRe")
    save([dataFolderName slash 'resRel'],"resRel")
    save([dataFolderName slash 'param'],"param")
    save([dataFolderName slash 'bestRel'],"bestRel")
    save([dataFolderName slash 'resCI'],"resCI")
    disp('Analysis done')
    
    clearvars paramsRep resGt grd grdAll re newData resCI5 resCI95 resCIrange
elseif data == "load"
    load([dataFolderName slash 'resRe'],"resRe")
    load([dataFolderName slash 'resRel'],"resRel")
    load([dataFolderName slash 'param'],"param")
    load([dataFolderName slash 'bestRel'],"bestRel")
    load([dataFolderName slash 'resCI'],"resCI")
    disp('Analysis data loaded')
elseif data == "empirical" || data == "empiricalNew"
    resRe = tracksRe(:,{'animal' 'method' 'thresh' 'threshNr'});
    
    resRe.nrTurns = cellfun(@(n) height(n)-2, tracksRe.track); % -2 to account for 1st & last point
    resRe.length = cellfun(@(n) sum(sqrt(diff([0; n.x]).^2 + diff([0; n.y]).^2)), tracksRe.track);
    resRe.angleMu = cellfun(@(n) mean(abs(n.alpha(2:end-1))), tracksRe.track); % Excluding 1st & end point
    resRe.angleSD = cellfun(@(n)  std(abs(n.alpha(2:end-1))), tracksRe.track);
    resRe.stepMu = resRe.length./resRe.nrTurns;
    resRe.stepSD = cellfun(@(n) std(n.s(2:end-1)), tracksRe.track);
    disp('Empirical data loaded')
end

%% Fig 2 'tracks': Tracks of example simulations
% Plots 2 example tracks of opposite parameters with zoom inset

zoomIdx = 1:50; % Which points to be displayed in zoomed version

% Selecting tracks
plSets = zeros(1,4);
% Step-and-turn tracks
plSets(1,1) = find(params.angleMu ==    param.angleMu(1) &...
                   params.angleSD ==    param.angleSD(1) &...
                   params.noiseAngle == param.noiseAngle(1) &...
                   params.smoothW ==    0);
plSets(1,2) = find(params.angleMu ==    param.angleMu(end) &...
                   params.angleSD ==    param.angleSD(1) &...
                   params.noiseAngle == param.noiseAngle(end-3) &...
                   params.smoothW ==    0);
% Curvy tracks
plSets(1,3) = find(params.angleMu ==    param.angleMu(1) &...
                   params.angleSD ==    param.angleSD(1) &...
                   params.smoothW ==    param.smoothW(2) &...
                   params.noiseAngle == 0);
plSets(1,4) = find(params.angleMu ==    param.angleMu(end) &...
                   params.angleSD ==    param.angleSD(1) &...
                   params.smoothW ==    param.smoothW(end) &...
                   params.noiseAngle == 0);

offs.x = [0   400 -400    0];
offs.y = [400   0    0 -400];
plTrackC = cell(4,1);
plTrZoomC = cell(4,1);
boxDim = zeros(4); % Corners of the zoom boxes over the original tracks

f2 = figure('units','centimeters','OuterPosition',[.2 1 19.05 18]);
hold on
for tr = 1:4 % Track = panel
    plTrackC{tr} = tracksWide{plSets(tr),1};
    plTrackC{tr}.x = plTrackC{tr}.x+offs.x(tr);
    plTrackC{tr}.y = plTrackC{tr}.y+offs.y(tr);
    plot(plTrackC{tr}.x,plTrackC{tr}.y,'b',lw{:})
    turnIdx = plTrackC{tr}.turn; % True turns marked w/ dot
    scatter(plTrackC{tr}.x(turnIdx),plTrackC{tr}.y(turnIdx),10,'blue','filled')
    % Zoom
    plTrZoomC{tr} = plTrackC{tr}(zoomIdx,:); % Zoomed-in parts
    boxDim(tr,:) = [min(plTrZoomC{tr}.x)   min(plTrZoomC{tr}.y),...
                        range(plTrZoomC{tr}.x) range(plTrZoomC{tr}.y)];
    rectangle('Position',boxDim(tr,:))
end
hold off
% axis equal
xlim([-515 515])
ylim([-515 515])
xlabel('x',ftsz{:})
ylabel('y',ftsz{:})

ax = gca;

% Insets of zoomed track (1 2 3 4 in 'Z' layout)
zoomPos = [.18 .7 .2 .2; .65 .7 .2 .2; .18 .2 .2 .2; .65 .2 .2 .2];
for tr = 1:4
    axes('Position',zoomPos(tr,:))
    box on
    plot(plTrZoomC{tr}.x,plTrZoomC{tr}.y,'b',lw{:})
    hold on
    zoomTurnIdx = plTrZoomC{tr}.turn; % True turns marked w/ dot
    scatter(plTrZoomC{tr}.x(zoomTurnIdx),plTrZoomC{tr}.y(zoomTurnIdx),10,'b','filled')
    hold off
    axis equal
    title(['Parameter set ' num2str(plSets(tr))])
end

% Lines for zoom illustration
% Corners of box around unzoomed track (rows: boxes, cols: corners clockwise, starting top left)
boxCornersAx.x = [boxDim(:,1) boxDim(:,1)+boxDim(:,3) boxDim(:,1)+boxDim(:,3) boxDim(:,1)];
boxCornersAx.y = [boxDim(:,2)+boxDim(:,4) boxDim(:,2)+boxDim(:,4) boxDim(:,2) boxDim(:,2)];
[boxCorners.x, boxCorners.y] = ds2nfu(ax, boxCornersAx.x, boxCornersAx.y);

% Corners of boxs around zoomed track (rows: boxes, cols: corners clockwise, starting top left)
zoomBoxCorners.x = [zoomPos(:,1) zoomPos(:,1)+zoomPos(:,3) zoomPos(:,1)+zoomPos(:,3) zoomPos(:,1)];
zoomBoxCorners.y = [zoomPos(:,2)+zoomPos(:,4) zoomPos(:,2)+zoomPos(:,4) zoomPos(:,2) zoomPos(:,2)];

% annotation('line',[0.372 0.500],[0.900 0.925]); % top left
annotation('line',[zoomBoxCorners.x(1,2) boxCorners.x(1,1)],[zoomBoxCorners.y(1,2) boxCorners.y(1,1)]); % top left
annotation('line',[zoomBoxCorners.x(1,3) boxCorners.x(1,1)],[zoomBoxCorners.y(1,3) boxCorners.y(1,3)]);
annotation('line',[zoomBoxCorners.x(2,3) boxCorners.x(2,2)],[zoomBoxCorners.y(2,3) boxCorners.y(2,2)]); % top right
annotation('line',[zoomBoxCorners.x(2,4) boxCorners.x(2,1)],[zoomBoxCorners.y(2,4) boxCorners.y(2,1)]);
annotation('line',[zoomBoxCorners.x(3,1) boxCorners.x(3,4)],[zoomBoxCorners.y(3,1) boxCorners.y(3,4)]); % bot Left
annotation('line',[zoomBoxCorners.x(3,2) boxCorners.x(3,3)],[zoomBoxCorners.y(3,2) boxCorners.y(3,3)]);
annotation('line',[zoomBoxCorners.x(4,1) boxCorners.x(4,2)],[zoomBoxCorners.y(4,1) boxCorners.y(4,2)]); % bot right
annotation('line',[zoomBoxCorners.x(4,4) boxCorners.x(4,3)],[zoomBoxCorners.y(4,4) boxCorners.y(4,3)]);

% exportgraphics(f2,[dataFolderName slash 'plots' slash 'f2_gtTracks.png'],'Resolution',300)
% exportgraphics(f2,[dataFolderName slash 'plots' slash 'f2_gtTracks.pdf'],'Resolution',300,'ContentType','vector')

%% Fig 3: Track & all its resampled tracks (best thresholds)
% Pick which GT track, methods, & thresholds
paramSet = 35; % High angle, high SD; High angle, high smooth (20, 35)
plSmoothThresh = 1; % Which thresh to show if curvy tracks (1st, 2nd, in threshs)
zoom = false; % For now: set to 'true' to run for zoomed-in on track start version
zoomW = 40; % Window length in steps around the origin

paramIdx = bestRel.paramSet==paramSet;
plTab = bestRel(paramIdx,:);
panelPre = methodSplitter(param); % Those resampling methods will be plotted (same layout as following plots)
panels = [panelPre{1}' panelPre{2}'];

% =============================== Plotting ============================== %
% Tracks separated into panels
ax = gobjects(size(panels));
f3 = figure('units','centimeters','OuterPosition',[.2 1 8 17.6]);
ti = tiledlayout(6,2,'TileSpacing','none', Padding='none');

% Ground Truth
ax(1) = nexttile(1);
plot(tracksWide{paramSet,1}.x,tracksWide{paramSet,1}.y,'b',lw{:})
hold on
turnInd = find(tracksWide{paramSet,1}.turn);
scatter(tracksWide{paramSet,1}.x(turnInd),tracksWide{paramSet,1}.y(turnInd),5,'b')
hold off
plLims = [min(tracksWide{paramSet,1}.x)-10 max(tracksWide{paramSet,1}.x)+10,...
          min(tracksWide{paramSet,1}.y)-10 max(tracksWide{paramSet,1}.y)+10];
axis([plLims(1:2) plLims(3)-30 plLims(4)]) % Shifting grd trth up to accommodate 'Method set 1' title
yline(plLims(3)+5,'LineWidth',1)
ax(1).XTickLabel = '';
ylabel('Ground Truth')
if zoom; axis([-zoomW zoomW -zoomW zoomW]); end

for pnlRow = 1:height(panels)
    for pnlCol = 1:width(panels)
        ax(pnlRow+1,pnlCol) = nexttile((1+pnlRow)*2 - 2 + pnlCol);
        % Threshold to be plotted: optimal determined in section above or rnd
        if params.smoothW(paramSet) == 0 % Step-and-turn GT tracks
            thrOpti = plTab.thresh(plTab.paramSet == paramSet & plTab.method == panels(pnlRow,pnlCol));
        else % Curvy tracks
            thrOpti = threshs.thresh(threshs.method == panels(pnlRow,pnlCol) & threshs.threshNr == plSmoothThresh);
        end
        reInd = find(tracksReWide.paramSet == paramSet & ...
                     tracksReWide.method == panels(pnlRow,pnlCol) & ...
                     tracksReWide.thresh == thrOpti(1));
        re = tracksReWide.track{reInd, 1};
        plot(re.x,re.y,'b',lw{:})
        if pnlRow<height(panels); ax(pnlRow+1,pnlCol).XTickLabel = ''; end % Removes x ticks of all but the last row
        ylabel(panels(pnlRow,pnlCol),'FontWeight','bold')
        axis(plLims)
        % axis('equal') % Doesn't fit that snuggly...
        if zoom; axis([-zoomW zoomW -zoomW zoomW]); end
    end
    ax(pnlRow+1,end).YAxisLocation = 'right';
end
% xlabel(ti,'x',ftsz{:}); ylabel(ti,'y',ftsz{:})
title(ax(2,1),'Method set 1'); title(ax(2,2),'Method set 2')

disp(['Turns '      num2str(params.angleMu(paramSet)) '°', newline,...
      'angle SD '   num2str(params.angleSD(paramSet)) '°', newline,...
      'noise max '  num2str(params.noiseAngle(paramSet)) '°', newline, ...
      'smooth window ' num2str(params.smoothW(paramSet)) ' steps'])

if plTab.param_smoothW == 0
    % exportgraphics(f3,[dataFolderName slash 'plots' slash 'f3_resTr.png'],'Resolution',300)
    % exportgraphics(f3,[dataFolderName slash 'plots' slash 'f3_resTr.pdf'],'Resolution',300,'ContentType','vector')
else
    % exportgraphics(f3,[dataFolderName slash 'plots' slash 'f3_resTrSmooth.png'],'Resolution',300)
    % exportgraphics(f3,[dataFolderName slash 'plots' slash 'f3_resTrSmooth.pdf'],'Resolution',300,'ContentType','vector')
end

%% Fig 3 for empirical tracks
if data ~= "empirical"; error('Run 1st sections w/ data = "empirical"'); end
zoom = true; % For now: set to 'true' and run for zoomed-in version
zoomW = 100; % Window length in steps around the origin
paramSet = 1; % Actually number of animal
animal = param.animal(paramSet); % A string of the name of the animal

panelPre = methodSplitter(param); % Those resampling methods will be plotted (same layout as following plots)
panels = [panelPre{1}' panelPre{2}'];
thrOpti = [5 4 2 4 3; 3 1 4 5 1]'; % Threshold to be plotted: where the angles are ~similar

% =============================== Plotting ============================== %
% Tracks separated into panels
ax = gobjects(size(panels));
f4 = figure('units','centimeters','OuterPosition',[.2 1 8 17.6]);
ti = tiledlayout(6,2,'TileSpacing','none','Padding','tight');

% Ground Truth
ax(1) = nexttile;
plot(tracks{paramSet}.x,tracks{paramSet}.y,'b',lw{:})
plLims = [min(tracks{paramSet,1}.x)-10 max(tracks{paramSet,1}.x)+10,...
          min(tracks{paramSet,1}.y)-10 max(tracks{paramSet,1}.y)+10];
axis([plLims(1:2) plLims(3)-(max(tracks{paramSet,1}.x)+10 - min(tracks{paramSet,1}.x)-10)/10 plLims(4)]) % Shifting grd trth up to accommodate 'Method set 1' title
yline(plLims(3)+5,'LineWidth',1)
ylabel('Raw animal track','FontWeight','bold')
if zoom
    xlim([tracks{paramSet}.x(1)-zoomW tracks{paramSet}.x(1)+zoomW])
    ylim([tracks{paramSet}.y(1)-zoomW tracks{paramSet}.y(1)+zoomW])
end

for pnlRow = 1:height(panels)
    for pnlCol = 1:width(panels)
        ax(1+pnlRow,pnlCol) = nexttile((1+pnlRow)*2 - 2 + pnlCol);
        reInd = find(tracksRe.animal == animal & ...
                     tracksRe.method == panels(pnlRow,pnlCol) & ...
                     tracksRe.threshNr == thrOpti(pnlRow,pnlCol));
        re = tracksRe.track{reInd};
        plot(re.x,re.y,'b',lw{:})
        if pnlRow<height(panels); ax(1+pnlRow,pnlCol).XTickLabel = ''; end
        axis(plLims)
        if zoom
            xlim([tracks{paramSet}.x(1)-zoomW tracks{paramSet}.x(1)+zoomW])
            ylim([tracks{paramSet}.y(1)-zoomW tracks{paramSet}.y(1)+zoomW]);
        end
        ylabel(panels(pnlRow,pnlCol),'FontWeight','bold')
    end
    ax(1+pnlRow,end).YAxisLocation = 'right';
end
% xlabel(ti,'x',ftsz{:}); ylabel(ti,'y',ftsz{:})
title(ax(2,1),'Method set 1'); title(ax(2,2),'Method set 2')

% exportgraphics(f4,[dataFolderName slash 'plots' slash 'f4_resampTrEmp.png'],'Resolution',300)
% exportgraphics(f4,[dataFolderName slash 'plots' slash 'f4_resampTrEmp.pdf'],'Resolution',300,'ContentType','vector')

%% Fig 4&5 & S SD,noise,length: Low accuracy & high variation (2 panels each, for method sets)
% Variables on the axes/panels/colors
xValName = 'noiseAngle'; % angleMu (fig 4), angleSD, noiseAngle (fig S 'noise'), or length (fig S 'length')
yValName = 'angleMu';

% Default fixed values that are plotted
length =     1000; % default for angleMu
angleMu =    50; % default for length
noiseXY =    0;
smoothW =    0;
angleSD =    1;
noiseAngle = 0;

% Selecting results of the selected parameter sets
switch xValName
    case 'angleMu'
        paramIdx = bestRel.param_noiseXY == noiseXY & bestRel.param_smoothW == smoothW &...
                   bestRel.param_angleSD == angleSD & bestRel.param_noiseAngle == noiseAngle & bestRel.param_length == length;
    case 'angleSD'
        paramIdx = ismember(bestRel.paramSet, find(params.angleSD>1 & params.angleMu == angleMu));
    case 'noiseAngle'
        paramIdx = ismember(bestRel.paramSet, find(params.noiseAngle>0 & params.angleMu == angleMu));
    case 'length'
        paramIdx = bestRel.param_noiseXY == noiseXY & bestRel.param_smoothW == smoothW &...
                   bestRel.param_angleMu == angleMu & bestRel.param_angleSD == angleSD & bestRel.param_noiseAngle == noiseAngle;
end

% Plot data
plTab = bestRel(paramIdx,{'method' 'thresh' 'paramSet' ['param_' xValName] yValName});

% Median line data
medTab = groupsummary(plTab,{'method' 'paramSet'},'median');
medTab.acc = abs(1 - medTab.(['median_' yValName]));

% Data for 5-95% interval of the data plotted in the first row
ciTab = resCI(ismember(resCI(:,["paramSet" "method" "thresh"]),plTab(:,["paramSet" "method" "thresh"]),"rows"),:);
ciTab.rnge95_dtw = nan(height(ciTab),1);

% Panels w/ >7 lines are hard to read. We split them up.
series = methodSplitter(param);

% =============================== Plotting ============================== %
jttrWdth = range(plTab.(['param_' xValName]))/numel(param.(xValName))/3; % Jitter width adjusted to x-axis scale
ax = gobjects(2, 2);
f4 = figure('units','centimeters','OuterPosition',[.2 1 19.05 20]);
ti = tiledlayout(2,2,'TileSpacing','tight','Padding','compact');
for panelRow = 1:2
    for methSet = 1:2
        ax(panelRow,methSet) = nexttile;
        hold on
        if panelRow == 1
            for s = 1:numel(series{methSet})
                sIdx = plTab.method==series{methSet}(s);
                jitter = jttrWdth.*rand(sum(sIdx),1) - jttrWdth/2; % Centering
                scatter(plTab{sIdx, ['param_' xValName]}+jitter, plTab{sIdx,yValName},ps,clr(s,:),'MarkerEdgeAlpha',.5)
                medIdx = medTab.method==series{methSet}(s);
                mp(s) = plot(medTab{medIdx, ['median_param_' xValName]}, medTab{medIdx,['median_' yValName]}, 'color',clr(s,:),lw{:});
            end
            if smoothW == 0
                yline(1,'--',lw{:})
            end
            title(['Method set ' num2str(methSet)])

            % Legends
            leg = legend(mp, series{methSet});
            title(leg,'Method')
            leg.ItemTokenSize = [10;0.1]; % Shortens lines in legend
        else
            for s = 1:numel(series{methSet})
                ciIdx = ciTab.method==series{methSet}(s);
                plot(ciTab{ciIdx, ['param_' xValName]}, ciTab{ciIdx,['rnge95_' yValName]}, 'color',clr(s,:),lw{:});
            end
        end
        hold off
    end
end

% Labels
if smoothW == 0
    ylabel(ax(1,1),'Turn angle accuracy (1=best)',ftsz{:})
    xlabel(ti,['Ground truth ' lblDct(xValName,param)])
else
    ylabel(ax(1,1),['Resampled ' lblDct(yValName,param)],ftsz{:})
    xlabel(ti,['Ground truth ' lblDct(xValName,param) ' before smoothing'])
end
ylabel(ax(2,1),['5-95% range of resampled ' newline lblDct(yValName,param)],ftsz{:})
ylabel(ax(1,2),[])

% Only if I wanted to link y axes
% linkaxes(ax(1,1:2),'y')
% linkaxes(ax(2,3:4),'y')
% set(ax(:,2),'YAxisLocation','left')

xticklabels(ax(1,1:2),[])
linkaxes(ax(:,1),'x'); linkaxes(ax(:,2),'x')


% Save
switch xValName
    case 'angleMu'
        % exportgraphics(f4,[dataFolderName slash 'plots' slash 'r1_accVar.png'],'Resolution',300)
        % exportgraphics(f4,[dataFolderName slash 'plots' slash 'r1_accVar.pdf'],'Resolution',300,'ContentType','vector')
        if smoothW > 0
            % exportgraphics(f4,[dataFolderName slash 'plots' slash 'r2_varCurvy.png'],'Resolution',300)
            % exportgraphics(f4,[dataFolderName slash 'plots' slash 'r2_accVar.pdf'],'Resolution',300,'ContentType','vector')
        end
    case 'angleSD'
        % exportgraphics(f4,[dataFolderName slash 'plots' slash 'S_angleSD.png'],'Resolution',300)
        % exportgraphics(f4,[dataFolderName slash 'plots' slash 'S_angleSD.pdf'],'Resolution',300,'ContentType','vector')
    case 'noiseAngle'
        % exportgraphics(f4,[dataFolderName slash 'plots' slash 'S_noise.png'],'Resolution',300)
        % exportgraphics(f4,[dataFolderName slash 'plots' slash 'S_noise.pdf'],'Resolution',300,'ContentType','vector')
    case 'length'
        % exportgraphics(f4,[dataFolderName slash 'plots' slash 'S_length.png'],'Resolution',300)
        % exportgraphics(f4,[dataFolderName slash 'plots' slash 'S_length.pdf'],'Resolution',300,'ContentType','vector')
end

%% Fig 6: Threshold susceptibility along noise/SD/smoothW
% Thresholds for each method across different params

yValName = 'angleMu'; % Metric on y axis (nrTurns, angleMu, angleSD, dtw)
xValName = 'thresh';
seriesName = 'noiseAngle'; % noiseAngle or angleSD or smoothW

fixedNames = ["noiseAngle" "angleSD" "smoothW"];
fixed = fixedNames(~ismember(fixedNames,string(seriesName)));

% Selecting correct tracks (of fixed values)
paramSetIdx = resRel.param_angleMu      == param.angleMu(1) &...
              resRel.param_stepMu       == param.stepMu(1) &...
              resRel.("param_" + fixed(1)) == param.(fixed(1))(1) &...
              resRel.("param_" + fixed(2)) == param.(fixed(2))(1);

plTab = resRel(paramSetIdx,:);
series = param.(seriesName); if strcmp(seriesName,'smoothW'); series = series(2:end); end
medTab = groupsummary(resRel(paramSetIdx,:),{['param_' seriesName] xValName 'method'},'median');

panelPre = methodSplitter(param);
panels = [panelPre{1}' panelPre{2}'];

% =============================== Plotting ============================== %
ax = gobjects(size(panels));
mp = gobjects(size(series));

f6 = figure('units','centimeters','OuterPosition',[.2 1 19.05 20]);
ti = tiledlayout(height(panels),width(panels),'TileSpacing','tight','Padding','compact');
for pnlRow = 1:height(panels)
    for pnlCol = 1:2
        ax(pnlRow,pnlCol) = nexttile;
        hold on
        for s = 1:numel(series)
            plIdx = plTab.method == panels(pnlRow,pnlCol) &...
                    plTab.(['param_' seriesName])==series(s);
            jttrWdth = repmat(max(diff(plTab.thresh(plIdx)))/10, sum(plIdx),1); % Adjusting to x-scale
            jitter = jttrWdth.*rand(sum(plIdx),1) - jttrWdth/2;
            scatter(plTab{plIdx,xValName}+jitter, plTab{plIdx,yValName}, ps,clrGr(s,:))
            medIdx = medTab.(['param_' seriesName])==series(s) &...
                     medTab.method==panels(pnlRow,pnlCol);
            mp(s) = plot(medTab{medIdx,xValName}, medTab{medIdx, ['median_' yValName]}, 'color',clrGr(s,:),lw{:});
        end
        if ~strcmp(yValName, 'dtw') && ~strcmp(seriesName,'smoothW') % DTW: lower is better; curvy: not relative
            yline(1,'--',lw{:})
        end
        hold off
        ylabel(ax(pnlRow,pnlCol),panels(pnlRow,pnlCol),"FontWeight","bold") % 'title' of that row
    end
    ax(pnlRow,end).YAxisLocation = 'right';
    linkaxes(ax(pnlRow,:),'y')
end

% Labels n stuff
title(ax(1,1),'Method set 1')
title(ax(1,2),'Method set 2')
xlabel(ti,'Threshold')
if strcmp(seriesName, 'smoothW')
    ylabel(ti,['Resampled ' lblDct(yValName,param)])
else % Step-and-turn
    ylabel(ti,'Turn angle accuracy (1=best)')
end

linkaxes(ax(panels=="local" | panels=="nonlocal"),'x') % For local method, which has some NaNs

leg = legend(mp,num2str(series'));
leg.ItemTokenSize = [10;18]; % Shortens lines in legend
title(leg,seriesName)

% exportgraphics(f6,[dataFolderName slash 'plots' slash 'r3_thresh_' seriesName '_' yValName '.png'],'Resolution',300)
% exportgraphics(f6,[dataFolderName slash 'plots' slash 'r3_thresh_' seriesName '_' yValName '.pdf'],'Resolution',300,'ContentType','vector')

%% Fig 7: Metrics within & between methods differ
% Between and within methods, there is a tradeoff between accuracies of
% different metrics

% methods in panels, x=metrics, series=threholds
paramSet = 7;
panels = param.method;

paramIdx = resRel.paramSet == paramSet;
plTabLong = resRel(paramIdx,10:end-1); % Currently excluding DTW bc of different y scale

plTab = stack(plTabLong, param.metric(1:end-1), "IndexVariableName","metricName", "NewDataVariableName","metric");
plTab.metricNameNum = double(plTab.metricName); % So that jitter can be applied
medTab = groupsummary(plTab,{'method' 'threshNr' 'metricNameNum'},'median',{'threshNr' 'metric'});

% =============================== Plotting ============================== %
ax = gobjects(1,numel(panels));
mp = gobjects(1,param.nrThreshs);
jttrWdth = .2; % Width of jitter

r4 = figure('units','centimeters','OuterPosition',[.2 1 19.05 20]);
ti = tiledlayout('flow','TileSpacing','tight','Padding','compact');
for panel = 1:param.nrMeths
    ax(panel) = nexttile;
    hold on
    for s = 1:param.nrThreshs
        sIdx = plTab.method == panels(panel) & plTab.threshNr == s;
        jitter = jttrWdth.*rand(sum(sIdx),1) - jttrWdth/2; % Centering
        scatter(plTab.metricNameNum(sIdx)+jitter, plTab.metric(sIdx),ps,clrGr(s,:))
        medIdx = medTab.method == panels(panel) & medTab.median_threshNr == s;
        mp(s) = plot(medTab{medIdx, 'metricNameNum'}, medTab{medIdx,'median_metric'}, 'color',clrGr(s,:),lw{:});
    end
    
    yline(1,'--',lw{:})
    hold off
    title(param.method(panel))
end
ylabel(ti,'Accuracy of the respective metric (1=best)')
xticklabels(ax(1:8),'')
xticks(ax(9:10),1:6)
xticklabels(ax(9:10),param.metric(1:6))

% Legends
leg = legend(mp,num2str(unique(plTab.threshNr)));
title(leg,'Threshold No.')
leg.ItemTokenSize = [10;0.1]; % Shortens lines in legend

% exportgraphics(r4,[dataFolderName slash 'plots' slash 'r4_metrics.png'],'Resolution',300)
% exportgraphics(r4,[dataFolderName slash 'plots' slash 'r4_metrics.pdf'],'Resolution',300,'ContentType','vector')

%% Fig 8: Empirical data

% Plot: panels = animals, x = thresholds, colors = methods
yValName = 'angleMu'; % angleMu or other metric
xValName = 'thresh';
plTab = resRe;

% Panels w/ >7 lines are hard to read. We split them up.
series = methodSplitter(param);

nrAnimals = numel(param.animal);

ax = gobjects(nrAnimals, numel(series));
f8 = figure('units','centimeters','OuterPosition',[.2 1 19.05 20]);
ti = tiledlayout(nrAnimals,2,'TileSpacing','tight','Padding','compact');
for pnlRow = 1:nrAnimals
    for methSet = 1:2
        ax(pnlRow,methSet) = nexttile;
        hold on
        plIdx = resRe.animal == param.animal(pnlRow);
        for s = 1:numel(series{methSet})
            sIdx = plIdx & plTab.method==series{methSet}(s);
            mp(s) = plot(resRe.threshNr(sIdx),resRe{sIdx,yValName},'-o',lw{:});
        end
        hold off
        if pnlRow == 1
            % Legends
            leg = legend(mp, series{methSet});
            title(leg,'Method')
            leg.ItemTokenSize = [10;0.1]; % Shortens lines in legend
        end
    end
    ax(pnlRow,end).YAxisLocation = 'right';
    ylabel(ax(pnlRow,end),param.animal(pnlRow),"FontWeight","bold") % 'title' of that row
    linkaxes(ax(pnlRow,:),'y')
end

title(ax(1,1), 'Method set 1'); title(ax(1,2), 'Method set 2')
linkaxes(ax(1,:))
linkaxes(ax(2,:))
xticklabels(ax(1,1:2),[])
xticks(ax(2,:),1:5)
ylabel(ti,['Resampled ' lblDct(yValName,param)],ftsz{:})
xlabel(ti,'Threshold Number',ftsz{:})

% exportgraphics(f8,[dataFolderName slash 'plots' slash 'r5_empirical.png'],'Resolution',300)
% exportgraphics(f8,[dataFolderName slash 'plots' slash 'r5_empirical.pdf'],'Resolution',300,'ContentType','vector')

%% SUPP Fig S_tracks: Plotting all simulated tracks (rep 1) in panels
% Currently plotted parameters
xPar = "angleMu";
yPar = "angleSD"; % angleSD
zPar = "smoothW"; % what's plotted w/in panels ("noiseAngle" or smoothW") 
zParPlot = [1; 2]; % [number]: those parameter numbers are plotted w/in panels

% Parameter list
xPars = unique(params{:,xPar});
yPars = unique(params{:,yPar});
zPars = unique(params{:,zPar});
nrxPars = numel(xPars); nryPars = numel(yPars);
% List of all possible combinations of x & y parameters
pars = array2table([repmat(xPars,numel(yPars),1) repelem(yPars, numel(xPars),1)],'variablenames',{'xPar' 'yPar'});

if any([xPar yPar zPar] == "smoothW"); curveInd = params.smoothW>0; % Only plot smoothed tracks
else; curveInd = params.smoothW==0; end % Only plot step-end-turn tracks

f2 = figure('units','centimeters','OuterPosition',[.2 1 19.05 10]);
ti = tiledlayout(nryPars,nrxPars,'tilespacing','tight', 'Padding','compact');
ax = gobjects(height(pars),1); % Tile handles to add labels
for i = 1:height(pars)
    ax(i) = nexttile(i);
    % Get indices of tracks in 'tracks' w/ current x, y, z params
    trIdx = params{:,xPar} == pars.xPar(i) &...
            params{:,yPar} == pars.yPar(i) & curveInd;
    trInd = find(trIdx & ismember(params{:,zPar}, zPars(zParPlot)));
    
    hold on % Plot all those in the same tile which have the specified zParam
    for z = 1:numel(trInd)
        plot(tracksWide{trInd(z),1}.x,tracksWide{trInd(z),1}.y)
    end
    hold off
    axis equal
end
leg = legend(num2str(zPars(zParPlot)));
title(leg,zPar)
leg.ItemTokenSize = [10;18]; % Shortens lines in legend
leg.Layout.Tile = 'east';

% Axis labels
for i = nrxPars*nryPars-nrxPars+1:nrxPars*nryPars
    ax(i).XLabel.String = num2str(pars.xPar(i));
end
for i = 1:nryPars+1:nryPars*nrxPars
    ax(i).YLabel.String = num2str(pars.yPar(i));
end
xlabel(ti,lblDct(char(xPar)))
ylabel(ti,lblDct(char(yPar)))

% Save
% exportgraphics(S1,[dataFolderName slash 'plots' slash 'S1_gtTracks.png'],'Resolution',300)
% exportgraphics(S1,[dataFolderName slash 'plots' slash 'S1_gtTracks.pdf'],'Resolution',300,'ContentType','vector')
%% BIG alternative to Fig 6: with 3 GT angles
% Thresholds for each method across different params

yValName = 'angleMu'; % Metric on y axis (nrTurns, angleMu, angleSD, dtw)
xValName = 'thresh';
seriesName = 'angleSD'; % noiseAngle or angleSD or smoothW
panelRowName = 'method';
panelColName = 'angleMu'; % angleMu or angleSD/noiseAngle (what isn't already seriesName)

% Fixed values
namesStr = [string(seriesName) string(xValName) string(panelRowName) string(panelColName)]; % All params somehow plotted
namesStrAll = ["angleMu" "angleSD" "noiseAngle" "method" "thresh" "noiseXY" "smoothW"]; % All params which could be changed
namesFixed = namesStrAll(~contains(namesStrAll, namesStr));

% Selecting correct tracks (of fixed values)
paramSetIdx = resRel.(strcat("param_", namesFixed(1))) == param.(namesFixed(1))(1) &...
              resRel.(strcat("param_", namesFixed(2))) == param.(namesFixed(2))(1) &...
              resRel.(strcat("param_", namesFixed(3))) == param.(namesFixed(3))(1);
              % (resRel.(['param_' panelColName]) == param.(panelColName)(1) |...
              %  resRel.(['param_' panelColName]) == param.(panelColName)(2) |...
              %  resRel.(['param_' panelColName]) == param.(panelColName)(3)); % Only selecting lowest, medium, and highest angle
plTab = resRel(paramSetIdx,:);
medTab = groupsummary(resRel(paramSetIdx,:),{['param_' seriesName] xValName panelRowName ['param_' panelColName]},'median');

series = param.(seriesName); if strcmp(seriesName,'smoothW'); series = series(2:end); end
panelRows = param.(panelRowName);
panelCols = param.(panelColName);

% =============================== Plotting ============================== %
ax = gobjects(numel(panelRows),numel(panelCols));
mp = gobjects(1,numel(series));

f6 = figure('units','centimeters','OuterPosition',[.2 1 19.05 20]);
ti = tiledlayout(numel(panelRows),3,'TileSpacing','tight','Padding','compact');
for pnlRow = 1:numel(panelRows)
    for pnlCol = 1:2:5%numel(panelCols)
        ax(pnlRow,pnlCol) = nexttile;
        hold on
        for s = 1:numel(series)
            plIdx = plTab.method == param.method(pnlRow) &...
                    plTab.(['param_' seriesName])==series(s) &...
                    plTab.(['param_' panelColName])==panelCols(pnlCol);
            jttrWdth = repmat(max(diff(plTab.thresh(plIdx)))/10, sum(plIdx),1); % Adjusting to x-scale
            jitter = jttrWdth.*rand(sum(plIdx),1) - jttrWdth/2;
            scatter(plTab{plIdx,xValName}+jitter, plTab{plIdx,yValName}, ps,clrGr(s,:))
            medIdx = medTab.(['param_' seriesName])==series(s) &...
                     medTab.(panelRowName)==panelRows(pnlRow) &...
                     medTab.(['param_' panelColName])==panelCols(pnlCol);
            mp(s) = plot(medTab{medIdx,xValName}, medTab{medIdx, ['median_' yValName]}, 'color',clrGr(s,:),lw{:});
        end
        if ~strcmp(yValName, 'dtw') && ~strcmp(seriesName,'smoothW') % DTW: lower is better; curvy: not relative
            yline(1,'--',lw{:})
        end
        hold off
    end
    ax(pnlRow,end).YAxisLocation = 'right';
    ylabel(ax(pnlRow,end),panelRows(pnlRow),"FontWeight","bold") % 'title' of that row
    linkaxes(ax(pnlRow,:),'y')
end

% Labels n stuff
xlabel(ti,'Threshold')
if strcmp(seriesName, 'smoothW')
    ylabel(ti,['Resampled ' lblDct(yValName,param)])
    title(ti,['Ground truth ' lblDct(panelColName,param) ' before smoothing'])
else % Step-and-turn
    ylabel(ti,['Relative ' lblDct(yValName,param)])
    title(ti,['Ground truth ' lblDct(panelColName,param)])
end

for col = 1:2:5; title(ax(1,col),[num2str(panelCols(col))]); end

yticks(ax(:,3),'')
linkaxes(ax(2,:),'x') % For local method, which has some NaNs

leg = legend(mp,num2str(series'));
title(leg,seriesName)
leg.ItemTokenSize = [10;18]; % Shortens lines in legend

disp([namesFixed', num2str([param.(namesFixed(1))(1);
                            param.(namesFixed(2))(1);
                            param.(namesFixed(3))(1)])]) % For caption

% exportgraphics(f6,[dataFolderName slash 'plots' slash 'f6_thresh_' seriesName '_' yValName '.png'],'Resolution',300)
% exportgraphics(f6,[dataFolderName slash 'plots' slash 'f6_thresh_' seriesName '_' yValName '.pdf'],'Resolution',300,'ContentType','vector')
