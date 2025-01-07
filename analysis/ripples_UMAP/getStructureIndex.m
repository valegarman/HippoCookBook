function [SI] = getStructureIndex(varargin)


%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);

parse(p,varargin{:})

basepath = p.Results.basepath;

cd('C:\Users\pabad\OneDrive - imim.es\Code\structure_index');

L = [1 1 1 1 1];
myListFile = pyrunfile("structure_index.py", "L")

end