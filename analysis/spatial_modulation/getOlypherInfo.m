function [olypherInfo] = getOlypherInfo(varargin)
% USAGE
%        [olypherInfo] = getOlypherInfo(vararging)
%
% INPUTS
%         <OPTIONALS>
%         firingTrialsMap - Output structure of getfFiringTrialsMap
%         round_to - integer value that data is discritized to. A value of
%                    2 means all data will be rounded to nearest 2 (i.e.
%                    2,4,8...)
%         smoothing - 0 no smoothing, else smooth with N bins
%         basepath - Default, pwd
%         saveMat
%         force
% OUTPUTS
%         olypherInfo
%
% Calculates the information carried in the firing rate of single
% neurons per spatial/temporal bin
%   
% Based on bz_olypherInfo.m
% MV-NCL 2024
%%

% Parse options
p = inputParser;
addParameter(p,'firingTrialsMap',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'smoothOpt',0,@islogical);
addParameter(p,'round_to',1,@isnumeric);
addParameter(p,'saveMat', false, @islogical);
addParameter(p,'force', false, @islogical);

parse(p, varargin{:});
basepath = p.Results.basepath;
firingTrialsMap = p.Results.firingTrialsMap;
smoothOpt = p.Results.smoothOpt;
round_to = p.Results.round_to;
saveMat = p.Results.saveMat;
force = p.Results.force;

% Deal with inputs
prevPath = pwd;
cd(basepath);

% create data matrix (M x N x D) where M is the number of cells to analyze,
% N is the number of trials for each cell, and D is the number of time bins
for ii = 1:length(firingTrialsMap.rateMaps{1})
    data = [];
    for jj = 1:length(firingTrialsMap.rateMaps)
        data = cat(3, data, firingTrialsMap.rateMaps{jj}{ii});
    end

    data = permute(data, [3, 1, 2]);
    
    M = size(data,1);
    N = size(data,2); 
    D = size(data,3);  
    pos_info_val = zeros(M,N,D);
    a = N*D;
    % Rounding 
    if smoothOpt ~= 0
        for i = 1 : M
            for k = 1:N
                data(i,k,:) = smooth(squeeze(data(i,k,:)),smoothOpt).*smoothOpt;
            end
        end
    end
    data = round(data./round_to)*round_to;

    for i = 1 : M
      for x = 1 : D   
         for k = 1:N
                q = data(i,k,x);
                pKx = ((length(find(data(i,:,x) == q)))/N);
                pK = ((length(find(data(i,:,:) == q)))/a);
                if pK == 0 || pKx == 0 || pKx < pK
                    pos_info_val(i,k,x) = pos_info_val(i,k,x);
                else
                    pos_info_val(i,k,x) = pos_info_val(i,k,x) + (pKx*log2(pKx/pK));
                end
         end        
      end
    end

    for i = 1:M
        track_info(i,:) = sum(pos_info_val(i,:,2:end-1),2);
    end
    
    % 
    olypherInfo.(['pos_trial_info_val_map' num2str(ii)]) = pos_info_val;
    olypherInfo.(['pos_info_map' num2str(ii)]) = track_info;
    
    olypherInfo.(['max_pos_info_map' num2str(ii)]) = max(track_info,[],2);
    olypherInfo.(['median_pos_info_map' num2str(ii)]) = median(track_info,2);
    olypherInfo.(['mean_pos_info_map' num2str(ii)]) = mean(track_info,2);
end


if saveMat
   save([basenameFromBasepath(pwd) '.positionalInformation.cellinfo.mat'],'olypherInfo'); 
end

cd(prevPath);

end