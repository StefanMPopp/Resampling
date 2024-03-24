function label = lblDct(in,param)
% Outputs the label for a plot axis, legend, title, etc. from the labelDict
label = cell2mat(param.labelDict({in}));
