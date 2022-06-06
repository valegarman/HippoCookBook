function uTest_ScaleTime(doSpeed)
% Automatic test: ScaleTime
% This is a routine for automatic testing. It is not needed for processing and
% can be deleted or moved to a folder, where it does not bother.
%
% uTest_ScaleTime(doSpeed)
% INPUT:
%   doSpeed: Optional logical flag to trigger time consuming speed tests.
%            Default: TRUE. If no speed test is defined, this is ignored.
% OUTPUT:
%   On failure the test stops with an error.
%
% Tested: Matlab 2009a, 2015b(32/64), 2016b, 2018b, Win7/10
% Author:  Jan Simon, Heidelberg, (C) 2009-2018 j@n-simon.de
% License: BSD

% $JRev: R-t V:019 Sum:5+4xLt3PCtZX Date:19-Oct-2020 12:32:57 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\UnitTests_\uTest_ScaleTime.m $
% History:
% 009: 20-Aug-2010 01:23, BUGFIX: isEqualTol was not included.
% 017: 12-Jul-2020 20:58, SINGLE as input.
% 019: 19-Oct-2020 11:44, Compare with griddedInterpolant.

% This is a dull test function, which creates a lot of variables without using
% them:
%#ok<*NASGU>

% Initialize: ==================================================================
ErrID = ['JSimon:', mfilename];
whichScaleTime = which('ScaleTime');

disp(['==== Test ScaleTime:  ', datestr(now, 0), char(10), ...
   'Version: ', whichScaleTime, char(10)]);

% Key for linear interpolation - the source of INTERP1 in Matlab6.5 states, that
% "*linear" is the old Matlab 5.3 method. But for large arrays the test of
% equally spaced steps is time-consuming!
% Set it to 'linear', if '*linear' does not work anymore.
linearMethod = '*linear';

if nargin == 0
   doSpeed = true;
end

% Look for lininterp1f, Umberto Picchini, Nov 2005
% http://www.mathworks.com/matlabcentral/fileexchange/8627
hasInterp1f = ~isempty(which('lininterp1f'));

% Look for qinterp1, Nathaniel Brahms,
% http://www.mathworks.com/matlabcentral/fileexchange/10286
hasQinterp1 = ~isempty(which('qinterp1'));

% Introduced in R2011b:
hasGriddedInterpolant = ~isempty(which('griddedInterpolant'));

% Do we use the MEX or M version:
[dummy1, dummy2, fileExt] = fileparts(whichScaleTime);  %#ok<ASGLU>
useMex = strcmpi(strrep(fileExt, '.', ''), mexext);

% Start tests: -----------------------------------------------------------------
for iClass = 1:2
   if iClass == 1
      T = 'double';
   else
      T = 'single';
   end
   fprintf('== Tests for class: %s\n', T);
   
   try
      desc = 'ScaleTime([], [])';
      y = ScaleTime(cast([], T), []);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isempty(y) && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime([], 1, 2, 0)';
      y = ScaleTime(cast([], T), 1, 2, 0);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isempty(y) && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime([], 1, 2, -1)';
      y = ScaleTime(cast([], T), 1, 2, -1);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isempty(y) && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime([1; 2], 1.5, 1, 1)';
      y = ScaleTime(cast([1; 2], T), 1.5, 1, 1);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isempty(y) && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime((1:10)'', [])';
      y = ScaleTime(cast((1:10)', T), []);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isempty(y) && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime(ones(10, 10), [])';
      y = ScaleTime(cast(ones(10, 10), T), []);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isempty(y) && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime((1:10)'', 1)';
      y = ScaleTime(cast((1:10)', T), 1);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isequal(y, 1) && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime(ones(10, 10), 1)';
      y = ScaleTime(cast(ones(10, 10), T), 1);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isequal(y, ones(1, 10)) && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime(ones(10, 10), 10)';
      y = ScaleTime(cast(ones(10, 10), T), 10);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isequal(y, ones(1, 10)) && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime([1; 2], 1, 1, 1)';
      y = ScaleTime(cast([1; 2], T), 1, 1, 1);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isequal(y, 1) && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime([1; 2], 1.5, 2, 1)';
      y = ScaleTime(cast([1; 2], T), 1.5, 2, 1);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isEqualTol_l(y, 1.5) && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime([1; 3], 1.6)';
      y = ScaleTime(cast([1; 3], T), 1.6);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isEqualTol_l(y, 2.2) && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime((1:10)'', 1, 10, 9)';
      y = ScaleTime(cast((1:10)', T), 1, 10, 9);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isEqualTol_l(y, (1:1.125:10)') && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime((1:10)'', 1, 10, 10)';
      y = ScaleTime(cast((1:10)', T), 1, 10, 10);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isEqualTol_l(y, (1:10)') && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime((1:10)'', 1, 10, 11)';
      y = ScaleTime(cast((1:10)', T), 1, 10, 11);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isEqualTol_l(y, (1:0.9:10)') && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime((1:10)'', 1, 9, 8)';
      y = ScaleTime(cast((1:10)', T), 1, 9, 8);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isEqualTol_l(y, (1:(8 / 7):9)') && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime((1:10)'', 1, 9, 9)';
      y = ScaleTime(cast((1:10)', T), 1, 9, 9);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isEqualTol_l(y, (1:9)') && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime((1:10)'', 1, 9, 10)';
      y = ScaleTime(cast((1:10)', T), 1, 9, 10);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isEqualTol_l(y, (1:8/9:9)') && isa(y, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime(1:10, 1, 9, 10) row vector';
      data = rand(1, 10);
      yRow = ScaleTime(cast(data, T),    1, 9, 10);
      yCol = ScaleTime(cast(data(:), T), 1, 9, 10);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isequal(yRow(:), yCol(:)) && isequal(size(yRow), [1, 10]) && ...
         isa(yRow, T) && isa(yCol, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime(3D-array)';
      data = rand(10, 20, 30);
      y1   = ScaleTime(cast(data, T), 1, 9, 10);
      y2   = ScaleTime(cast(reshape(data, 10, []), T), 1, 9, 10);
      y2R  = reshape(y2, [], 20, 30);
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isequal(y1, y2R) && isa(y2R, T)
      disp(['  ok: ', desc]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   % Compare with INTERP1: -----------------------------------------------------
   disp([char(10), '== Compare with INTERP1:']);
   try
      desc = 'ScaleTime(sin, 1, 2*pi, 100)';
      data = cast(sin(0:0.1:2*pi)', T);
      y    = ScaleTime(data, 1, 2*pi, 100);
      yM   = interp1(1:size(data, 1), data, linspace(1, 2*pi, 100)');
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isEqualTol_l(y, yM)
      disp(['  ok: ', desc, ' == INTERP1']);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime([sin, sin(down:up)], 1, 2*pi, 100)';
      data = cat(2, data, data(end:-1:1));
      y    = ScaleTime(data, 1, 2*pi, 100);
      yM   = interp1(1:size(data, 1), data, linspace(1, 2*pi, 100)');
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isEqualTol_l(y, yM)
      disp(['  ok: ', desc, ' == INTERP1']);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   
   try
      desc = 'ScaleTime(rand(100, 100), 1, 100, 57)';
      data = cast(rand(100, 100), T);
      y    = ScaleTime(data, 1, 100, 57);
      yV   = ScaleTime(data, 1:99/56:100);
      yM   = interp1(1:size(data, 1), data, linspace(1, 100, 57)');
   catch ME
      error(ErrID, ['Crash: ', desc, char(10), ME.message]);
   end
   if isEqualTol_l(y, yM, 100)
      disp(['  ok: ', desc, ' == INTERP1']);
      disp(['      max deviation: ', ...
         sprintf('%.1f * eps', max(abs(y(:) - yM(:))) / eps)]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
   if isEqualTol_l(y, yV, 100)
      disp(['  ok: ', desc, ' == ScaleTime(rand(100, 100), 1:99/56:100)']);
      disp(['      max deviation: ', ...
         sprintf('%.1f * eps', max(abs(yM(:) - yV(:))) / eps)]);
   else
      error(ErrID, ['Bad reply from: ', desc]);
   end
end

% Speed: -----------------------------------------------------------------------
if doSpeed
   TestTime = 0.5;  % sec
   fprintf('\n== Speed test (test time: %g sec):\n', TestTime);
else
   TestTime = 0.01;
   fprintf('\n== Speed test (test time: %g sec - may be inaccurate):\n', ...
      TestTime);
end

DataLenList   = [10, 100, 1000, 10000, 100000];
DataWidthList = [1, 10, 100];
InterpLenList = [10, 100, 1000, 10000];

for DataLen = DataLenList
   for DataWidth = DataWidthList
      data = rand(DataLen, DataWidth);
      
      fprintf('Data size: [%d x %d]\n', DataLen, DataWidth);
      fprintf('  Interpolation steps:    ');
      fprintf('  %7d', InterpLenList);
      
      % Matlab's INTERP1:
      fprintf('\n    INTERP1:              ');
      for InterpLen = InterpLenList
         NewT = 1:(DataLen - 1) / (InterpLen - 1):DataLen;
         
         tic;
         N = 0;
         eTime = 0;
         while eTime < TestTime || N < 5
            y = interp1(1:DataLen, data, NewT, linearMethod);
            N = N + 1;
            eTime = toc;
         end
         PrintLoop(N / eTime);
      end
      fprintf('  loops/sec\n');
      drawnow;
      
      % griddetInterpolant including its creation:
      if hasGriddedInterpolant
         fprintf('    griddedInterpolant:   ');
         for InterpLen = InterpLenList
            NewT  = 1:(DataLen - 1) / (InterpLen - 1):DataLen;
            Index = 1:DataWidth;
            
            Method = 'linear';
            ExtraP = 'none';
            if DataWidth == 1  % GI({NewT, 1}) does not work...
               tic;
               N = 0;
               eTime = toc;
               while eTime < TestTime || N < 5
                  GI = griddedInterpolant(data, Method, ExtraP);
                  y  = GI(NewT);
                  N  = N + 1;
                  eTime = toc;
               end
               
            else  % DataLen > 1: GI({NewT, Index})
               tic;
               N = 0;
               eTime = toc;
               while eTime < TestTime || N < 5
                  GI = griddedInterpolant(data, Method, ExtraP);
                  y  = GI({NewT, Index});
                  N  = N + 1;
                  eTime = toc;
               end
            end
            PrintLoop(N / eTime);
         end
         fprintf('  loops/sec  (inclusive creation of Interpolant)\n');
         drawnow;

         % griddedInterpolant without its creation:
         fprintf('    griddedInterpolant:   ');
         for InterpLen = InterpLenList
            NewT  = 1:(DataLen - 1) / (InterpLen - 1):DataLen;
            Index = 1:DataWidth;
            
            Method = 'linear';
            ExtraP = 'none';

            % Before the loop:
            GI = griddedInterpolant(data, Method, ExtraP);
               
            if DataWidth == 1  % GI({NewT, 1}) does not work...
               tic;
               N = 0;
               eTime = toc;
               while eTime < TestTime || N < 5
                  y  = GI(NewT);
                  N  = N + 1;
                  eTime = toc;
               end
               
            else  % DataLen > 1: GI({NewT, Index})
               tic;
               N = 0;
               eTime = toc;
               while eTime < TestTime || N < 5
                  y  = GI({NewT, Index});
                  N  = N + 1;
                  eTime = toc;
               end
            end
            PrintLoop(N / eTime);
         end
         fprintf('  loops/sec  (exclusive creation of Interpolant)\n');
         drawnow;
      end
      
      % Compare with local Matlab version:
      fprintf('    ScaleTime.m:          ');
      for InterpLen = InterpLenList
         NewT = 1:(DataLen - 1) / (InterpLen - 1):DataLen;
         
         tic;
         N = 0;
         eTime = toc;
         while eTime < TestTime || N < 5
            y = ScaleTime_local(data, NewT);
            N = N + 1;
            eTime = toc;
         end
         PrintLoop(N / eTime);
      end
      fprintf('  loops/sec  (no checks of inputs)\n');
      drawnow;
      
      % MEX with vector input:
      if useMex
         fprintf('    ScaleTime.mex/vector: ');
      else
         fprintf('    ScaleTime.m/vector:   ');
      end
      for InterpLen = InterpLenList
         NewT = 1:(DataLen - 1) / (InterpLen - 1):DataLen;
         
         tic;
         N = 0;
         eTime = toc;
         while eTime < TestTime || N < 5
            y = ScaleTime(data, NewT);
            N = N + 1;
            eTime = toc;
         end
         PrintLoop(N / eTime);
      end
      fprintf('  loops/sec\n');
      drawnow;
      
      % MEX with index method:
      if useMex
         fprintf('    ScaleTime.mex/index:  ');
      else
         fprintf('    ScaleTime.m/index:    ');
      end
      for InterpLen = InterpLenList
         tic;
         N = 0;
         eTime = toc;
         while eTime < TestTime || N < 5
            y = ScaleTime(data, 1, DataLen, InterpLen);
            N = N + 1;
            eTime = toc;
         end
         PrintLoop(N / eTime);
      end
      fprintf('  loops/sec\n');
      drawnow;
      
      % Compare with functions from Matlab's FEX:
      if hasInterp1f && DataWidth == 1
         fprintf('    lininterp1f:          ');
         for InterpLen = InterpLenList
            NewT = 1:(DataLen - 1) / (InterpLen - 1):DataLen;
            
            tic;
            N = 0;
            eTime = toc;
            while eTime < TestTime || N < 5
               y = lininterp1f(1:DataLen, data, NewT, []);
               N = N + 1;
               eTime = toc;
            end
            PrintLoop(N / eTime);
         end
         fprintf('  loops/sec\n');
         drawnow;
      end
      
      if hasQinterp1 && DataWidth == 1
         fprintf('    qinterp1:             ');
         for InterpLen = InterpLenList
            NewT = 1:(DataLen - 1) / (InterpLen - 1):DataLen;
            
            tic;
            N = 0;
            eTime = toc;
            while eTime < TestTime || N < 5
               y = qinterp1(1:DataLen, data, NewT, 1);
               N = N + 1;
               eTime = toc;
            end
            PrintLoop(N / eTime);
         end
         fprintf('  loops/sec\n');
         drawnow;
      end
      
      fprintf('\n');
   end
end

% Bye:
fprintf('\n== ScaleTime passed the tests.\n');

end

% ******************************************************************************
function PrintLoop(N)
if N > 10
   fprintf('  %7.0f', N);
else
   fprintf('  %7.1f', N);
end

end

% ******************************************************************************
% Local copy of ScaleTime to test speed if compiled MEX is in the path:
function Yi = ScaleTime_local(Y, Ti, Tf, Tn)
% Author: Jan Simon, Heidelberg, (C) 2009-2020 j@n-simon.de
% Was JRev: R0h V:033 Sum:98748BA6 Date:30-Sep-2009 00:50:13 $
% Was File: Tools\GLMath\ScaleTime.m $
persistent hasAutoExpand
if isempty(hasAutoExpand)
   matlabV       = [100, 1] * sscanf(version, '%d.%d', 2);
   hasAutoExpand = (matlabV > 901);
end

[nRow, nCol] = size(Y);

if nargin == 4
   if Tn < 1.0
      Yi = [];
      return;
   elseif Tn > 1
      Ti = Ti:((Tf - Ti) / (Tn - 1)):Tf;
   elseif Tf < Ti
      Ti = [];
   end
end

if isempty(Ti)
   Yi = Y([]);
   return;
end

Ti = Ti(:);
Si = Ti - floor(Ti);
Ti = floor(Ti);

d     = (Ti == nRow);
Ti(d) = Ti(d) - 1;
Si(d) = 1;

if nCol > 1 % && ~hasAutoExpand
   Si = Si(:, ones(1, nCol));
end
Yi = Y(Ti, :) .* (1 - Si) + Y(Ti + 1, :) .* Si;

end

% ******************************************************************************
function Equal = isEqualTol_l(x, y, F)
% Compare two double arrays with absolute tolerance
% Equal = isEqualTol_l(x, y, Tol)
% If the maximal difference between two double arrays is smaller than Tol, they
% are accepted as equal.
% As in Matlab's ISEQUAL, NaNs are treated as inequal, so isEqualTol_l(NaN, NaN)
% is FALSE. If you need comparable NaNs use ISEQUALWITHEQUALNANS, although this
% name is horrible and it is not existing in Matlab5.3.
%
% Author: Jan Simon, Heidelberg, (C) 2009-2020 j@n-simon.de

% Was JRev: R0c V:022 Sum:mP3jthZjqvhS Date:12-Mar-2010 01:41:33
% Was File: Tools\GLMath\isEqualTol_l.m

% Simplified version!
if nargin < 3
   F = 10;
end
Tol = F * eps(class(x));

% Try the exact and fast comparison at first:
Equal = false;
if isequal(size(x), size(y))
   % Same as "if all(abs(xMy(:)) <= Tol)", but faster:
   xMy = x - y;
   if all(or((abs(xMy) <= Tol), (x == y)))   % is FALSE for NaNs
      Equal = true;
   end
end

end
