function tracksWide = trackMaker(param,params)
% Inputs: meta parameters, GT track parameters
% Output: track w/ x,y,t,id,turn angle,heading angle, turn(bool)
nrReps = param.nrReps;
tracksWide = cell(height(params),nrReps);
for setNr = 1:height(params)
    % === Unpacking parameters === %
    l = params.length(setNr); % Length of tracks (in 'noise steps')
    % Step length
    stepDistr = params.stepDistr(setNr); % turn-step-length distribution: norm, exp
    stepMu = params.stepMu(setNr); % Mean (multiple of noise step-length)
    stepSD = params.stepSD(setNr); % SD for gaussian

    % Turn angle
    angleMu = params.angleMu(setNr); % 
    angleSD = params.angleSD(setNr);
    
    % Noise steps & Smoothing
    smoothW = params.smoothW(setNr); % Smoothing window
    noiseAngle = params.noiseAngle(setNr); % Max abs noise angle
    noiseXY = params.noiseXY(setNr); % x, y displacement noise SD
    
    % ========================== %
    % === Making simulations === %
    
    for rep = 1:nrReps
        % === Turns === %
        % Turns (w/ or w/o noise)
        a = zeros(l,1); % a = α = turn angles at each point
        if angleMu == 0 % No turns → only noise angle
            turnInd = false(l,1);
        else
            turnInd = sMaker(stepMu, stepSD, l, stepDistr); % Where the α must be changed
            turnInd(turnInd>l) = []; % Way more turn inds are made above than fit in the track
            nrOfTurns = length(turnInd);
            % turn = vmrand(angleMu,angleKappa,nrOfTurns,1);
            turn = normrnd(angleMu,angleSD,nrOfTurns,1); % Creates turn angles to be assigned
            turn(1:round(nrOfTurns/2)) = -turn(1:round(nrOfTurns/2)); % Half turn right, half left
            turn = turn(randperm(length(turn))); % Above, the first half turned right. Now random succession
            a(turnInd) = turn; % Replaces angle at point turnInd in track w/ respective angle made in normrnd
        end
            
        % === Heading angles & noise/smoothing === %
        % Smoothing w/ LOWESS
        if smoothW > 0
            a = smooth((1:l), a, smoothW, 'loess'); % x vs t
            turnInd = [2;l+1]; % Ain't no turns in a smooth track
        end
        theta = cumsum(a); % Initial heading theta=0 (=right (south))
        
        % Angle noise (from animal). (Note: positional noise comes below)
        if noiseAngle > 0 % Adds heading jitter to each step
            % % noise drawn from vanMises (K would need to be inverted)
            % theta = wrapTo180(theta + rad2deg(vmrand(0,noiseAngleMax,l,1)));
            theta  = wrapTo180(theta + normrnd(0,noiseAngle,l,1)); % Gaussian noise
        end
    
        % === Making x,y from heading angles === %
        s = ones(l,1); % Unit step length
        % For future speed expansion: v & α correlation: s around turnInd shall be squished (lin or exp)
        dx = s.*cosd(theta); % Initial position [0,0]
        dy = s.*sind(theta);
    
            % Positional noise (e.g. from GPS)
            dx = dx + normrnd(0,noiseXY,l,1);
            dy = dy + normrnd(0,noiseXY,l,1);
    
        sim = struct();
        sim.x = cumsum(dx);
        sim.y = cumsum(dy);
        sim.t = (1:l)';
        sim.id = repmat(setNr,l,1);
        sim.alpha = a;
        sim.theta = theta;
        turnIdx = false(l,1); turnIdx(turnInd-1) = true;
        sim.turn = turnIdx;
        sim = struct2table(sim);
    
        tracksWide{setNr,rep} = sim;
    end
end