% mexSptxModal help file for mexSptxModal MEX file.
%
% Finds best 2-way split in matrix of floating point numbers along 2 of the
% dimensions/columns, hereafter referred to as measurements.
%
%
%   OUTPUT ARGUMENTS
%
%       X             - an integer for the 1st of 2 measurements where best
%                       split is found
%
%       Y             - an integer for the 2nd of 2 measurements
%
%       polygon       - the polygon of the best split found on the above 2
%                       measurements.
%
%	To get the 2nd best split, the API invoker adds in a set of the above
%	3 output arguments.  Getting the 3rd best split requires a 3rd set etc.
%
%
%   REQUIRED INPUT ARGUMENT
%       data - a table of float (AKA single) numbers normalized from 0 to 
%		1 ... any out of range numbers are censored.
%
%       
%   OPTIONAL NAME-VALUE ARGUMENTS
%
%     'balanced'   -   A boolean value indicating the split goal.  If
%                      false, then favor splits with the least
%                      weight/density along there edges. If true, then
%                      favor splits that are similarly sized.
%		               The default is true.
%
%     'W'	   -       Standard deviation of kernel.  This is the highest
%                      achievable resolution; in practice a higher value
%	                   might be used for application reasons or just
%                      performance.
%                      The default is .006.
%
%    'sigma'	    -  controls the density threshold for starting a new
%                      cluster.
%                      The default is 3.0.
%
%    'KLD_normal_1D' - Kullback-Leibler Divergence (KLD) test to determine
%                      each measurement's informativeness and whether it is
%                      worth using in split.
%                      The default is .16.
%
%    'KLD_normal_2D' - Is a particular pair of dimensions worth splitting? 
%                      The default is .16.
%
%    'KLD_exponential_1D' - Is this an exponential tail (e.g. CyTOF)?
%                      The default is .16.
%    
%    'max_clusters'  - The most clusters the graph logic should handle.
%                      The default is 12.
%    
%    'verbose_flags' - Determines console output.
%
%    'threads'       - the # of threads to use  0 means no threads
%                      and -1 means s many threads as available
%                      hardware cores on the computer coing the split  
%                      The default is -1.
%
%    'simplify_polygon' - a boolean indicqating the need to smooth the 
%			polygon using the Ramer–Douglas–Peucker algorithm
%                      The default is true.

%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause