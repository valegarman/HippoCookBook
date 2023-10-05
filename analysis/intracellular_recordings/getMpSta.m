

function [mpSta] = getMpSta(varargin)
% [mpSta] = getMpSta(varargin)
%
% Perform Spike trigger averages (STA) over intracell data
%
% INPUTS
% <optional>
% basepath      Folder containing a spikes.cellinfo.mat file and a
%                   timeSeries.intracellular.mat file and a spikes.cellInfo.mat
%                   Default pwd.
% intracellular Intracellular time series with at least:
%                   - .data containing intracellular signal 
%                       [nSamples x nChannels].
%                   - .timestamps [nSamples x 1] vector 
%                       with timestamps. 
%                   - (Optionally) .cellInts (start stop times)
%                       indicating different cells segments [start_cell_1, 
%                       end_cell_1; start_cell2, end_cell2; ...].
% spikes         buzcode ripple structure (from bz_GetSpikes). If not provided, 
%                   it loads it from 'basepath' (if provided), or from current 
%                   folder (if not)
% sr            (scalar) Sampling frequency in HZ (by default intracellular.sr
%                   or 20000).
% doPlot        (logical) Creates plot to inspect STA results (default
%                   false).
% STAInts       2 x 1 vector with STA time interval, default [-100 100]
% surr          (logical/scalar) If true or  scalar (N, default N=500) create
%                   surrogate template with N repetitions.
% highPass      Passband frequency (in Hz) for the high pass filter (default 20).
% jitterWin     Jitter windows size, in ms (default 15)
% exclusion     (start stop times) interval matrix indicating segment of d 
%                   (in s) where STA shouldn't be computed, ie SPW-r.
% force         By default (false), load any .mpSta.intracellular.mat on the
%                   folder.
% saveMat       True or false save a buzcode event file with results
%                   (default true)
% saveSurrogates
%               Save complete surrogates matrix (N x StaInts size). False
%               by default.
%
% OUTPUTS
% mpSta         Structure containing all STA results:
% .data         MATLAB cell [1 x N, being N # of different intracellular 
%                   cell signals] containing STA from extracellular
%                   spikes [timestamps x # of extracellular clusters].
% .timestamps   [nSamples x 1] vector with timestamps. 
% .staJittered Same to .data but STA after surrogates substraction.
% .staJitteredStd 
%              Same to .data but containing the dispersion [std] of the
%                   surrogates.
% .allTraces   Structure containing all data for STA computing, including
%                   all STA traces (.data), STA sweeps membrane potential
%                   (.membranePotential), STA sweeps curren (.current),
%                   spikes empited on the sweep (true or false, on .subThre)
%                   and extracellular spikes used on each STA (.spk_times).
% .processingInfo
%              Structure containing information about the STA computing,
%              including sampling rate (.sr), high pass band (.highPass),
%              jitter window for surrogates (.jitterWin) and intervals
%              exclusiongs (.exclusions)
%
%   Manu Valero 2018
%   Updated on 2020 for buzcode compatibility.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'intracellular',[],@isstruct);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'sr',20000,@isscalar);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'STAInts',[]);
addParameter(p,'surr',true);
addParameter(p,'highPass',20,@isscalar);
addParameter(p,'jitterWin',15,@isscalar);
addParameter(p,'exclusion',[]);
addParameter(p,'force',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveSurrogates',true,@islogical)

parse(p,varargin{:});
basepath = p.Results.basepath;
intracellular = p.Results.intracellular;
spikes = p.Results.spikes;
sr = p.Results.sr;
doPlot = p.Results.doPlot;
STAInts = p.Results.STAInts;
surr = p.Results.surr;
highPass = p.Results.highPass;
jitterWin = p.Results.jitterWin;
exclusion = p.Results.exclusion;
force = p.Results.force;
saveMat = p.Results.saveMat;
saveSurrogates = p.Results.saveSurrogates;

% Dealing with inputs
prevBasepath = pwd;
cd(basepath);

staFile = dir('*.mpSta.intracellular.mat');
if ~isempty(staFile) && ~force
    disp('Membrane potential STA already computed! Loading file.');
    load(staFile.name);
    return
end

if isempty(intracellular)
    fileIntracellular = dir('*timeSeries.intracellular.mat');
    try load(fileIntracellular.name,'intracellular');
    catch
        error('Intracellular file not found!');
    end
end
d = intracellular.data;

if isempty(spikes)
    spikes = loadSpikes;
end

if isfield(intracellular, 'sr')
    sr = intracellular.sr;
end

if isempty(STAInts)
    STAInts = [-100 100];
end
STAInts = STAInts * sr/1000;

if islogical(surr) || surr == 1
    surr = 100;
end

cellsInts = intracellular.cellsInts;

% Hight pass filtering
fprintf('%iHz high pass filtering... \n',highPass);
hpFilt = designfilt('highpassiir','FilterOrder',3, 'PassbandFrequency',highPass,'PassbandRipple',0.1, 'SampleRate',sr);
d_filt = filtfilt(hpFilt,d(:,1));
timestamps_intracell = intracellular.timestamps;

% RUN STA
disp('Running STA (could take some minutes)...');
N = STAInts(2) - STAInts(1) + 1; % STA samples size
for ii = 1 : size(cellsInts,1)
    fprintf('Cell %i out of %i... \n',ii, size(cellsInts,1));    
    for jj = 1 : size(spikes.times,2)
        fprintf('  Cluster %i out of %i... \n',jj, size(spikes.times,2));
        tsCTemp = spikes.times{jj}(spikes.times{jj} >= cellsInts(ii,1) & spikes.times{jj} <= cellsInts(ii,2));
        if numel(tsCTemp) > 1 & isempty(find(tsCTemp(1) == timestamps_intracell))
            temp = []; count = 0;
            disp('  Cell timestamps to samples...');
            for kk = 1:length(tsCTemp)
                lineLength = fprintf('    %i out of %i... \n',kk, length(tsCTemp));
                [~,temp(kk)] =  min(abs(timestamps_intracell - tsCTemp(kk)));
                fprintf(repmat('\b',1,lineLength));
            end
            tsCTemp = timestamps_intracell(temp);
        end
        tsCTemp_times = tsCTemp;
        [~,tsCTemp] = intersect(timestamps_intracell,tsCTemp);
        
        if numel(exclusion) > 0
            tmpExc  = zeros(size(tsCTemp));
            for kk = 1 : size(exclusion,1)
                tmpExc = tmpExc + (tsCTemp > exc(kk,1) & tsCTemp < exc(kk,2));
            end
            tsCTemp(tmpExc > 0) = [];
        end
        staTemp = zeros(N,numel(tsCTemp));
        mpTemp = zeros(1,numel(tsCTemp));
        cTemp = zeros(1,numel(tsCTemp));
        subThrTemp = zeros(1,numel(tsCTemp));
        spk_time = zeros(1,numel(tsCTemp));        
        tsCTemp(tsCTemp < abs(STAInts(1)) | tsCTemp > size(d,1)-STAInts(2)) = [];
        for kk = 1 : numel(tsCTemp)
            staTemp(:,kk) = d_filt(int32(tsCTemp(kk) + STAInts(1) : tsCTemp(kk) + STAInts(2)));
            mpTemp(kk) = median(d(int32(tsCTemp(kk) + STAInts(1) : tsCTemp(kk) + STAInts(2)),1));
            if size(d,2) >1
                cTemp(kk) = median(d(int32(tsCTemp(kk) + STAInts(1) : tsCTemp(kk) + STAInts(2)),2));
            else
                cTemp(kk) = NaN;
            end
            subThrTemp(kk) = ~any(d_filt(int32(tsCTemp(kk) + STAInts(1) : tsCTemp(kk) + STAInts(2)))>10 | ...
                d_filt(int32(tsCTemp(kk) + STAInts(1) : tsCTemp(kk) + STAInts(2)))<-10);
            spk_time(kk) = double(tsCTemp(kk))/sr;
        end
        
        if ~isempty(tsCTemp) && surr % create surrogate
            disp('    Computing surrogates...');
            surrTemp = zeros(N,surr); 
            parfor mm = 1 : surr
                u_s = int32(jitter_surrogate(double(tsCTemp(find(subThrTemp>0))),jitterWin,'fs',sr));
                u_s(u_s<abs(STAInts(1))) = []; u_s(u_s > length(d_filt) - STAInts(2)) = [];
                s_temp = zeros(N,numel(tsCTemp(find(subThrTemp>0))));
                for nn = 1 : numel(tsCTemp(find(subThrTemp>0)))
                    try s_temp(:,nn) = d_filt(u_s(nn) + STAInts(1) : u_s(nn) + STAInts(2)); 
                    end
                end
                surrTemp(:,mm) = mean(s_temp,2);
            end
        else
            surrTemp = [];
        end
        
        staData.dfilt{ii}{jj} = staTemp;
        staData.mp{ii}{jj} = mpTemp;
        staData.c{ii}{jj} = cTemp;
        staData.subThre{ii}{jj} = subThrTemp;
        staData.surr{ii}{jj} = surrTemp;
        staData.spk_times{ii}{jj} = tsCTemp_times;
    end
end

timestamps = (STAInts(1)/sr:1/sr:STAInts(2)/sr)';
% computing traces
disp('Getting traces...');
mpSta = [];
mpSta.timestamps = timestamps;
for ii = 1 : size(cellsInts,1)
    for jj = 1 : size(spikes.times,2)
        if ~isempty(staData.dfilt{ii}{jj})
            idx = find(staData.subThre{ii}{jj});
            mpSta.data{ii}(:,jj) = mean(staData.dfilt{ii}{jj}(:,idx),2);
            mpSta.dataJittered{ii}(:,jj) = mean(staData.dfilt{ii}{jj}(:,idx),2) - mean(staData.surr{ii}{jj},2);
            mpSta.jitterStd{ii}(:,jj) = std(staData.surr{ii}{jj},[],2);
        else
            mpSta.data{ii}(:,jj) = nan(1,length(timestamps));
            mpSta.dataJittered{ii}(:,jj) = nan(1,length(timestamps));
            mpSta.jitterStd{ii}(:,jj) = nan(1,length(timestamps));
        end
    end
end
mpSta.allTraces.data = staData.dfilt;
mpSta.allTraces.membranePotential = staData.mp;
mpSta.allTraces.current = staData.c;
mpSta.allTraces.subThre = staData.subThre;
mpSta.allTraces.spk_times = staData.spk_times;
mpSta.processinginfo.sr = sr;
mpSta.processinginfo.highPass = highPass;
mpSta.processinginfo.jitterWin = jitterWin;
mpSta.processinginfo.exclusion = exclusion;
if saveSurrogates
    mpSta.allTraces.surr = staData.surr;
end

if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.mpSta.intracellular.mat'],'mpSta','-v7.3');
end

if doPlot
    inspectMpSta;
end
cd(prevBasepath);

end

function ts_jitter = jitter_surrogate(ts,delta,varargin)
% ts_jitter = jitter_surrogate(ts,delta,'fs',_,'nb_jitter',_)
%
% Load kilosort/phy clusters
%
% INPUTS
% ts                (R x 1 vector) Presynaptic spike train in samples.
% delta             (scalar) jitter interval.
%                   
% <options> Optional list of property-value pairs (see table below)
%   fs              Sampling frequency, default 20000.
%   nb_jitter       Number of surrogates (better to give a perfect square),
%                   default 1
% 
% OUTPUTS
% ts_jitter         R x nb_jitter matrix with jittered spike train in
%                   samples.
%
%   Manu Valero 2018

% Default values
fs = 20000;
nb_jitter = 1;

% Parse options
for ii = 1:2:length(varargin)
    if ~ischar(varargin{ii})
		error(['Parameter ' num2str(ii+3) ' is not a property (type ''help <a href="matlab:help jitter_surrogate">jitter_surrogate</a>'' for details).']);
    end
    switch(lower(varargin{ii}))
        case 'fs'
            fs = varargin{ii+1};
        case 'nb_jitter'
            nb_jitter = varargin{ii+1};
    end
end

if size(ts,1) < size(ts,1)
    ts = ts';
end

% Prepare jitter surrogates
nb_jitt = sqrt(nb_jitter);
ts_jitter = repmat(ts, 1, nb_jitt);
delta = delta * fs/1000;

% Jitter surrogates
% ts_jitter = delta * floor(ts_jitter/delta) + delta * rand(nb_jitt,length(ts))'; % Code from Sam, jitter distribution is a pyramid, why?
ts_jitter = ts_jitter + delta * (rand(length(ts), nb_jitt) - 0.5); % Flat jitter distribution

end