function reTrack = resampler(track, thresh1, thresh2, method)
% Resampling a track according to 1 of 4 threshold-based methods
% Inputs: track (individual): table with columns x,y,t,id,v,theta,alpha
%         thresh(old): scalar or vector [spatial thresh, angle thresh]
%         method: to be used ('interval, 'local', 'non-local', 'circvar')
% Output: resampled track

switch method
    case "corridor" % Turchin et al. 1991
        % Points are eliminated when they are <d distance away from the straight
        % line connecting the first and last point of the current segment

        p1 = 1; % Start point of current sleeve
        p2 = 3; % End point of current sleeve
        i = 1; % Counter of points in the resampled track
        reTrack(1,:) = track(1,1:4);
        while p2 < height(track)
            p = p2-1;
            d = point2LineDist(track{p1+1:p2-1,1:2},track{p1,1:2},track{p2,1:2});
            if any(d > thresh1)
                reTrack(i,:) = track(p,1:4);
                p1 = p;
                p2 = p2+1;
                i = i+1;
            end
            p2 = p2+1;
        end
        
    case "TPA" % Resample w/ Turning Points Algorithm (Circular Variance) (Potts 2018)
        window_size = thresh1; % Set window size
        minThresh = thresh2; % Set threshold angle in degrees
        reTrack = circVarResamplingFun(track, window_size, minThresh);

    case "DP"
        [~, reInd] = dpsimplify(track,thresh1);
        reTrack = track(reInd,1:4);

    case "vDP" % Speed-modified Douglas-Peucker (Thiebault 2013)
        track.v = ones(size(track.x));
        theta360 = wrapTo360(track.theta);
        x = track.v .* cos(theta360);
        y = track.v .* sin(theta360);
        d1 = diff(x);
        d2 = diff(y);
        [~, reIndx] = dpsimplify([track.t(2:end),d1],thresh1);
        [~, reIndy] = dpsimplify([track.t(2:end),d2],thresh1);
        reInd = sort([reIndx; reIndy]);
        reInd([false; diff(reInd)<4]) = [];
        reTrack = track(reInd,1:4);

    case "MRPA"
        reTrack = track(TDMRPA_SED(track.x,track.y,track.t,thresh1, 2,2),1:4);

    case "interval" % Resample @ a fixed point interval
        reTrack = track(1:thresh1:end,1:4);

    case "local" % Resample w/ Local method/STTrace
        a = [NaN;
            diff([NaN; sqrt(sum(diff(track{:,{'x' 'y'}},1).^2,2))] ./... % Distance
                 [NaN; diff(track.t)]) ./ diff(track.t) % time interval
            ];
        reTrack = track(a>=thresh1 | abs(track.alpha)>=thresh2,1:4);
        
    case "nonlocal" % Resample w/ Non-local method (Reynolds 2007)
        % Keep only points where cumsum turn angle hits threshold
        a = track.alpha(~isnan(track.alpha));
        ptInd = 1;
        currPtInd = 1;
        i = 1; % Loop counter
        stop = 0; % Flips 1 when the loop should end
        while stop == 0
            i = i+1;
            cuma = cumsum(a(currPtInd+ptInd(i-1):end));
            currPtInd = find(abs(cuma)>thresh2,1);
            if ~isempty(currPtInd)
                ptInd(i,1) = currPtInd+ptInd(i-1);
            else
                stop = 1;
                ptInd(i,1) = height(track);
            end
        end
        reTrack = track(ptInd,1:4);
        
    case "CPT" % Change Point Test (Byrne et al. 2009), isosceles triangles
        q = thresh1; % Step-length, I guess
        sigLvl = thresh2; % Significance level (.05 or .01 in source)
        reTrack = cptFun(track, q, sigLvl);

    case "1D" % Humphreys 2013 "3D 1D symmetry". Only aims to get s dists
        reInd = find(islocalmax(track.x) | islocalmin(track.x));
        reInd([false; diff(reInd) < thresh1]) = []; % thresh = min s. In source determined through GOF
        reTrack = track(reInd,1:4);

    case "1D+" % Improved 1D, by Tromer et al 2015: 1. resample w/ lamda, 2. run 1D on resampled.
        presamp = (1:thresh1:height(track))'; % thresh = lamda
        reInd = presamp(islocalmax(track.x(presamp)) | islocalmin(track.x(presamp)));
        reTrack = track(reInd,1:4);

    case "VW" % Visvalingam-Whyatt (Time-Sensitive) (Hunnik 2017, Utrecht University)
        % Make triangles for each consecutive 3 points
        reTrack = track(:,1:4);
        area = [inf; triangleArea(reTrack); inf]; % Custom function, padding 1st & last point
        areaMin = min(area);
        % Eliminate the mid point of the triangle w/ the smalles area
        % Repeat until the smallest triangle has an area >threshold
        while areaMin < thresh1
            [areaMin, minIdx] = min(area);
            reTrack(minIdx,:) = [];
            area(minIdx) = [];
            if minIdx < numel(area)-2
                area(minIdx-1:minIdx) = triangleArea(reTrack(minIdx-1:minIdx+2,1:3));
                area(1) = inf; % In case the 2nd triangle is smallest
            elseif minIdx == numel(area)-2
                area(minIdx-1) = triangleArea(reTrack(minIdx-1:minIdx+1,1:3));
            end
        end

    case "inflections" % Resample @ turn direction change points
        % ToDo: add thresholds of length/angle or combination of both to avoid
        % biol. non-significant ultra short wiggles
        warning('Threshold not implemented yet. Wiggles will be broken up')
        theta = [NaN; atan2d(track.y(2:end)-track.y(1:end-1), track.x(2:end)-track.x(1:end-1))]; % heading angle
        alpha = [rad2deg(angdiff(deg2rad(theta))); NaN]; % Turning angle
        r = alpha > 0; % Right turn
%             l = alpha <= -thresh2; % Left turn
%             s = alpha<thresh & alpha>-thresh2; % Straight
        inflection = diff(r)~=0; % | diff(l)~=0 | diff(s)~=0;
        reTrack = track(inflection,:);
end

if isempty(reTrack) % No turns detected
    reTrack = [track(1,1:4); track(end,1:4)];
end
if any(reTrack{1,:} ~= track{1,1:4}) % First point was not carried over
    reTrack = [track(1,1:4); reTrack];
end
if any(reTrack{end,:} ~= track{end,1:4}) % Last point was not carried over
    reTrack = [reTrack; track(end,1:4)];
end
