function [tetrode,clust,mtint] = getTetAndClust( varargin )

if nargin == 0
    try
        [flnmroot,filepath] = uigetfile({'*.set';'*.mat'},'Select the file to load...',...
            'MultiSelect','on');
        % see if the file is a previously saved mtint structure
        if iscell(flnmroot)
            mtint = readAllDACQdata(filepath,flnmroot);
        elseif strcmpi(flnmroot(findstr(flnmroot,'.')+1:end),'mat')
           load([filepath,flnmroot],'-mat');
        elseif strcmpi(flnmroot(findstr(flnmroot,'.')+1:end),'set')
            mtint = readAllDACQdata(filepath,flnmroot);
        end
    catch
    end
elseif nargin == 1
    mtint = varargin{1};
end
ntets = numel(mtint.tetrode);
tets = zeros(1,ntets);
for itet = 1:ntets
    tets(1,itet) = mtint.tetrode(itet).id;
end
tets = cellstr(num2str(tets'));
tetrode = listdlg('PromptString', 'Select a tetrode',...
    'ListString',tets,'SelectionMode','single');
allClusters = cellstr(num2str(unique(mtint.tetrode(tetrode).cut)));
clust = listdlg('PromptString', 'Select a clust',...
    'ListString',allClusters,'SelectionMode','multiple');
clust = str2double(allClusters(clust));
end