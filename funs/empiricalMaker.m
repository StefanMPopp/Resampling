clear
if ispc; slash = '\'; else; slash = '/'; end % OS compatibility
dirName = 'resample'; % Name of the base folder
addpath(genpath(dirName))
dirEmp = [dirName slash 'empirical'];

% ======================================================================= %
% Ants in arena, equidistantly resampled (each step 2mm long)
raw = readtable([dirName slash 'empirical' slash 'HRM_W1_ants.txt']);
raw.alpha = deg2rad(raw.alpha);
tracks{1,1} = raw(raw.id ==  2,{'x' 'y' 't' 'id' 'v' 'alpha' 'theta'}); tracks{1,1} = tracks{1,1}(1:1000,:);
tracks{2,1} = raw(raw.id == 18,{'x' 'y' 't' 'id' 'v' 'alpha' 'theta'}); tracks{2,1} = tracks{2,1}(1:1000,:);
save([dirEmp slash 'tracks'],'tracks')

dataFolderName = [dirName slash '2023-11-29_bigNoLength']; % Where the param file to be copied is located
% param
load([dataFolderName slash 'param'])
param.animal = ["temno1"; "temno2"];
param.length = 1500;
nrTrheshs = 5;
param = rmfield(param,["nrReps" "stepDistr" "stepMu" "stepSD" "angleMu" "angleSD" "noiseAngle" "noiseXY" "smoothW"]);
save([dirEmp slash 'param'],"param")

% threshs
threshs = load([dataFolderName slash 'threshs']);
threshs{threshs.method=="nonlocal","thresh2"} = (5:4:21)';
save([dirEmp slash 'threshs'],"threshs");