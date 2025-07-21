
% script for embeding

% First, I have created an environment in python (in anaconda terminal, run)
conda create -n matlab_env python=3.10.12
conda activate matlab_env
pip install numpy==1.23.5 scipy==1.15.2 scikit-learn==1.2.2 joblib==1.4.2 matplotlib==3.10.1
% then, I copy the path obtained with (runned in the terminal) here
where python 
% which in my case was
py_path = 'C:\Users\mvalero\.conda\envs\matlab_env\python.exe';

% now, we can run it
embeding_path = 'C:\Users\mvalero\OneDrive - imim.es\Documents\Code\HippoCookBook\utilities\python\compute_lfp_embeding.py';
command = ['"' py_path '" "' embeding_path '"'];
[status, output] = system(command);
disp(output)

% it saves a basepath.lfp_embedding.channelinfo.mat file that can be load
% into matlab
load('python.lfp_embedding.channelinfo.mat'); % (the name of my folder is python)

% --- 1. Channel to trajectory
traj = lfp_embeding.model_trajectory;
projs = lfp_embeding.data_projection_2d;
precision = lfp_embeding.precision;


n_channels = size(projs, 1);
trajis = zeros(n_channels, 1);

for i = 1:n_channels
    dists = sum((traj - projs(i,:)).^2, 2);
    [~, idx] = min(dists);
    trajis(i) = idx;
end

% --- 2. Reference points
ctrl_pts = lfp_embeding.all_training_data_points;
ctrl_labels = cellstr(strtrim(string(lfp_embeding.all_training_data_labels)));

ctrl_traj = zeros(size(ctrl_pts,1), 1);
for i = 1:size(ctrl_pts,1)
    dists = sum((traj - ctrl_pts(i,:)).^2, 2);
    [~, idx] = min(dists);
    ctrl_traj(i) = idx;
end