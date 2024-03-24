function series = methodSplitter(param)
% Panels w/ >7 lines are hard to read. We split them up into 2 panels
% Input: param (nrMeths, method)
% Output: 1x2 cell w/ names of ~1/2 of the methods in each cell

methHalf = ceil(param.nrMeths/2);
if param.nrMeths == 10 % Final set?
    param.method(5) = "CPT"; % Switching the two around b/c CPT inaccurate, should be 1st panel
    param.method(8) = "corridor"; % corridor is better
    param.method(2) = "VW";
    param.method(10) = "local";
end
series = {param.method(1:methHalf) param.method(methHalf+1:end)}; % Lines of different colors. Named after Excel