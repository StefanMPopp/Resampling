function reTrack = circVarResamplingFun(track, window_size, thresh)
% Resample w/ Circrular Variance of Turn Angle method (Potts et al 2018)
% Keep only points whose heading 'spikes' in a sliding window
% Inputs: tracks in one table (must contain x,y,t,id,theta)
%         sliding window size (W in the paper)
%         threshold below which points are discarded (theta_{thresh} in the
%         paper)
% Output: Resampled track (x,y,t,id)

%{
Potts JR, Börger L, Scantlebury DM, Bennett NC, Alagaili A, Wilson RP.
Finding turning-points in ultra-high-resolution animal movement data.
Methods Ecol Evol. 2018;9:2091–2101
%}
% %% Code translated from R code %%

headings = track.theta(~isnan(track.theta)); % Heading in degrees

% Calculate sin and cosine of headings, together with x- and y-values 
% (assuming a constant speed)
cosines = cos(headings(:,1)*pi/180);
sines = sin(headings(:,1)*pi/180);
x_vals = cumsum(cosines);
y_vals = cumsum(sines);

% Find average cos and sin over sliding window, as well as 
% squared circular standard deviations (SCSD)
sd_length = length(cosines) - window_size;
ave_cos_array = cosines(1:sd_length)/window_size;
ave_sin_array = sines(1:sd_length)/window_size;
for counter = 1:(window_size-1)
  % Shift cosine and sin arrays back by 1
  cosines = cosines(2:end);
  sines = sines(2:end);
  ave_cos_array = ave_cos_array + cosines(1:sd_length)/window_size;
  ave_sin_array = ave_sin_array + sines(1:sd_length)/window_size;
end
circ_sd = (-2)*log(sqrt((ave_sin_array).^2 + (ave_cos_array).^2));
ave_circ_sd = sum(circ_sd) / length(circ_sd);


% Find candidate turning points by looking for spikes in SCSD
turning = 0;
chgpt_array = 1;
for counter = 1:length(circ_sd)
  if (circ_sd(counter) > ave_circ_sd) && turning == 0
    % Started turning.  Note turning point
    start_chgpt = counter;
    turning = 1;
  end
  
  if (circ_sd(counter) < ave_circ_sd) && (turning == 1)
    % Turning point is the mean of the start turning-point and the end shifted to the right by resolution/2 
    % to account for averaging being done forwards in time
    chgpt = floor((start_chgpt + counter)/2 + window_size/2);
    chgpt_array = [chgpt_array,chgpt];
    turning = 0;
  end
end

chgpt_array = [chgpt_array, length(cosines)];

% Post processing of changepoints to remove any where the switch in direction is less than thresh_angle
new_chgpt_array = 1;
prev_heading = atan2(y_vals(chgpt_array(2))-y_vals(chgpt_array(1)),...
                     x_vals(chgpt_array(2))-x_vals(chgpt_array(1)));
for counter = 2:(length(chgpt_array)-1)
    next_heading = atan2(y_vals(chgpt_array(counter+1))-y_vals(chgpt_array(counter)),...
                         x_vals(chgpt_array(counter+1))-x_vals(chgpt_array(counter)));
    if abs(next_heading - prev_heading) > pi
        turn_angle = pi*2-abs(next_heading - prev_heading);
    else
        turn_angle = abs(next_heading - prev_heading);
    end
    if turn_angle*180/pi > thresh
        new_chgpt_array = [new_chgpt_array, chgpt_array(counter)];
        prev_heading = next_heading; 
    end
end
new_chgpt_array(end+1) = height(track); % Last point also included

reTrack = track(new_chgpt_array,1:4);
