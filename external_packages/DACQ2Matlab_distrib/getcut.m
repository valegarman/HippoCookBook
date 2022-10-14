function [varargout]  = getcut(cutfile)

% Read cut file from tint
% Usage: clust = getcut(cutfile)
%
%
% Sturla Molden <sturla@molden.net>
% Centre for the Biology of Memory
% Norwegian University of Science and Technology
% http://www.cbm.ntnu.no
% 
% Copyright (C) 2003  Centre for the Biology of Memory, NTNU
% All Rights Reserved
%
% This M-file is released under Q Public License v 1.0,
% a copy of which should accompany this file.

if strcmpi(computer,'PC')
    fid = fopen(cutfile, 'rt');
else
    fid = fopen(cutfile,'r');
end
clust = [];
while ~feof(fid)
    string = fgetl(fid);
    if isempty(string)
    elseif string(1) == 'E'
        exact_text = string;
        break;
    end
end

while ~feof(fid)
  string = fgetl(fid);
  if length(string)>0
     content = sscanf(string,'%u')';
     clust = [clust content];
  end
end
fclose(fid);
clust = clust';
varargout{1} = clust;
varargout{2} = exact_text;