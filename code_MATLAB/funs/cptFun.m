function out = cptFun(track,alpha,q)
% R code for performing the "circular change point test" (CPT)
% as described in "How did they get here from there?
% Detecting changes of direction in terrestrial ranging."
% by R W Byrne, R G Noser, L A Bates & P E Jupp. 
% Animal Behaviour 77, 619-631, 2009.

% Translated to MATLAB code by Stefan Popp 2023, for resampling paper
% Code for the original version is available at
% http://www.mcs.st-andrews.ac.uk/~pej/CPT_Rcode (needs waybackwhenmachine)

% This is an "automated" version in which 
% (a) q  and  alpha  are set by the user;
% (b) the parameters N (number of permutations) and
% tol (maximum distance between indistinguishable 
% positions) may be changed by the user; 
% (c) the output is a list (in cps) of (row numbers and coordinates of)
% waypoints that are detected as  significant

% A manual for this code can be found at 
% http://www.mcs.st-andrews.ac.uk/~pej/CPTauto_Rcode_Manual.pdf.

% P. E. Jupp 26 June 2012
% Code supplied without guarantee
x1 = track{:,1};
x2 = track{:,2};
xy = [x1, x2];

% INPUTS

% alpha = significance level
% alpha = 0.05;

% Enter value of q
% q = 4 is used here just as an example.
% The user should replace this by a value of q that is appropriate for the
% species and data considered, e.g. that obtained by running CPT_Rcode
% on part of the data.
% q = 20;

% input N
% N = total number of permutations (1 observed and N-1 simulated)
% N = 10000 is a convenient number
N = 99;

% input tol
% tol = tolerance = maximum distance between indistinguishable positions
% tol = 0 is a convenient default
% The user may like to replace this by twice the average GPS error in research area
tol = 0;

% PRELIMINARIES

% Reverse the time-ordering,
% so that (bx1[1], bx2[1]) refers to (final)
% putative goal
bx1 = zeros(size(x1));
bx2 = zeros(size(x2));
n = height(x1);
for j = 1:n
    bx1(j) = x1(n - j + 1);
    bx2(j) = x2(n - j + 1);
end

% Calculate the steps (bxdiff1, bxdiff2)
bx1diff = diff(bx1);
bx2diff = diff(bx2);

% REMOVE POINTS AT WHICH ANIMAL STAYS STILL
% ind is (reverse) time ordering of points
% newp  > 0 if point differs from previous point
tolsq = tol^2;
ind = 1:n;
newp = 1:n;
for j = 2:n
    newp(j) = ind(j) * (bx1diff(j - 1)^2 + bx2diff(j - 1)^2 > tolsq);
end

% bz1, bz2 are coordinates of points(in reverse time order) at which there is movement
bz1 = bx1(newp > 0);
bz2 = bx2(newp > 0);
nz = length(bz1);

% We shall apply the CPT to points with coordinates (bz1,bz2)

% Calculate the steps (bzdiff1, bzdiff2)
bz1diff = diff(bz1);
bz2diff = diff(bz2);

% SOME DECLARATIONS

% goal_no = number (from end) of current putative goal (with goal_no = 1 for end position
% Thus goal_no = goal.t + 1
% where goal.t was used in original code to refer to
% time of current putative goal (with goal.t = 0 for # end position
% start with goal_no = 1
goal_no = 1;

% last_no = number (backwards in time) of last position of interest
% last_no = length(bz1) is a convenient default
% (and includes all the positions)
last_no = numel(bz1);

no_of_nos = numel(bz1);

Rsumrand = zeros(1, N);
Rsumrandr = zeros(1, N);

% Pr will store observed p-values in a run of r
Pr = ones(1, no_of_nos);

% sig is a vector indicating whether or not 
% waypoint k is detected as a possible change point:
% sig(k) = 1 if change point detected at waypoint k
% sig(k) = 0 otherwise
sig = zeros(1, nz);

% k is number of steps in "k-leg"
k = 0;

% LOOK (SEQUENTIALLY) FOR NEXT POSSIBLE CHANGE POINT
while(goal_no < last_no - q) % A
    k = 0;
    
    % INCREASE k UNTIL NEXT POSSIBLE CHANGE POINT IS FOUND
    
    % Re-initialise Pr
    Pr = ones(1, no_of_nos);
    P = 1;
    
    while (P > alpha) && (goal_no + q + k < last_no) % B
        k = k + 1;
        
        R1 = sqrt((bz1(goal_no+k  ) - bz1(goal_no  ))^2 + (bz2(goal_no+k  ) - bz2(goal_no  ))^2);
        R2 = sqrt((bz1(goal_no+k+q) - bz1(goal_no+k))^2 + (bz2(goal_no+k+q) - bz2(goal_no+k))^2);
        
        Rsum = R1 + R2;
        
        % Rsumrand(1) = observed value of statistic R1 + R2
        Rsumrand(1) = Rsum;
        
        % Now calculate statistic R1 + R2 for a further N-1 random permutations
        % and store in Rsumrand
        for it = 2:N % C
            perm = randperm(nz-goal_no-1)';

            bz1r = bz1(goal_no);
            bz2r = bz2(goal_no);

            for j = 1:k % D
                bz1r = bz1r + bz1diff(goal_no-1+perm(j));
                bz2r = bz2r + bz2diff(goal_no-1+perm(j));
            end

            Rsumrand(it) = sqrt((bz1r - bz1(goal_no))^2 + (bz2r - bz2(goal_no))^2) +...
                           sqrt((bz1(goal_no+k+q) - bz1r)^2 + (bz2(goal_no+k+q) - bz2r)^2);
        end
        
        P = sum(Rsumrand >= Rsum) / N;
        
        Pr(k) = P;
        
    end
    
    % f is the first value of k in the current run that is significant
    f = k;
    
    while (P <= alpha) && (goal_no + q + k < last_no) % E
        k = k + 1;
        R1 = sqrt((bz1(goal_no+k) - bz1(goal_no    ))^2 + (bz2(goal_no+k) - bz2(goal_no    ))^2);
        R2 = sqrt((bz1(goal_no+k+q) - bz1(goal_no+k))^2 + (bz2(goal_no+k+q) - bz2(goal_no+k))^2);

        Rsum = R1 + R2;
        
        % Rsumrand(1) = observed value of statistic R1 + R2
        Rsumrand(1) = Rsum;
        
        % Now calculate statistic R1 + R2 for a further N-1 random permutations
        % and store in Rsumrand
        for it = 2:N % F
            perm = randperm(nz-goal_no-1)';
            bz1r = bz1(goal_no);
            bz2r = bz2(goal_no);
            
            for j = 1:k % G
                bz1r = bz1r + bz1diff(goal_no-1+perm(j));
                bz2r = bz2r + bz2diff(goal_no-1+perm(j));
            end
            
            Rsumrand(it) = sqrt((bz1r - bz1(goal_no    ))^2 + (bz2r - bz2(goal_no    ))^2) +...
                           sqrt((bz1(goal_no+k+q) - bz1r)^2 + (bz2(goal_no+k+q) - bz2r)^2);
        end
        
        P = sum(Rsumrand >= Rsum) / N;
        
        Pr(k) = P;
    end
    % l is the last value of k in the current run that is significant
    l = k - 1;

    % apply "peak rule"
    rmin = find(Pr(f:l) == min(Pr(f:l)),1);
    if l < f
        rmin = 0;
    end
    
    if rmin > 0
        goal_no = goal_no + f + rmin - 1;
    else
        goal_no = last_no;
    end

    if goal_no == last_no
        sig(goal_no) = 0;
    else
        sig(goal_no)= 1;
    end
end

sig(1) = 1;
sig(end) = 1;

% OUTPUTS

% bsig is sig in reverse time-order
bsig = flip(sig);

% (first, last) (row nos. of first and last times at "significant" change points) 
% (east, north) (their coordinates) 
cp.ind = find(bsig == 1)';

% "times" (row nos.) and coordinates of change points

% newpp = newp(newp > tol);
% cp.time = newpp(bsig == 1)';
cp.time = track{:,3}(sig == 1);
cp.bx1 = flip(bz1(bsig == 1));
cp.bx2 = flip(bz2(bsig == 1));
cp.no = n + 1 - cp.time;


out = array2table([cp.bx1 cp.bx2 cp.time ones(numel(cp.bx1),1).*track.id(1)],VariableNames=track.Properties.VariableNames(1:4));
