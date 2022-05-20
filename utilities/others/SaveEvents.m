%SaveEvents - Write events to file.
%
%  USAGE
%
%    SaveEvents(filename,events)
%
%    filename            event file name
%    events              event data
%
%  SEE
%
%    See also NewEvents, LoadEvents, SaveRippleEvents.

% Copyright (C) 2004-2006 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function SaveEvents(filename,events)

if exist(filename)
    disp('File already exists. Aborting.')

else


    file = fopen(filename,'w');
    if file == -1,
	    error(['Cannot write to ' filename]);
    end
    
    for i = 1:length(events.time),
	    fprintf(file,'%f\t%s\n',events.time(i)*1000,events.description{i}); % Convert to milliseconds
    end
    
    fclose(file);
end
