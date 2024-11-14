function [fcsdat, fcshdr, fcsdatscaled] = fca_readfcs(filename, ...
    headerOnly, numericKeys, numericDefaults, booleanKeys,...
    booleanDefaults, stringKeys, stringDefaults, maxEvents, ...
    isSideScatterLinear,doLogicle, startAtEvent)
% [fcsdat, fcshdr, fcsdatscaled] = fca_readfcs(filename);
%
% Read FCS 2.0, 3.0 and 3.1 flow cytometry data file and put the list mode  
% parameters to the fcsdat array with size of [NumOfPar TotalEvents]. 
% Some important header data are stored in the fcshdr structure:
% TotalEvents, NumOfPar, starttime, stoptime and specific info for parameters
% as name, range, bitdepth, logscale(yes-no) and number of decades.
%
% [fcsdat, fcshdr] = fca_readfcs;
% Without filename input the user can select the desired file
% using the standard open file dialog box.
%
% [fcsdat, fcshdr, fcsdatscaled] = fca_readfcs(filename);
% Supplying the third output the fcsdatscaled array contains the scaled     
% parameters. It might be useful for logscaled parameters, but no effect 
% in the case of linear parameters. The log scaling is the following
% operation for the "ith" parameter:  
% fcsdatscaled(:,i) = ...
%   10.^(fcsdat(:,i)/fcshdr.par(i).range*fcshdr.par(i).decade;);

% Ver 2.5
% 2006-2009 / University of Debrecen, Institute of Nuclear Medicine
% Laszlo Balkay 
% balkay@pet.dote.hu
%
% 14/08/2006 I made some changes in the code by the suggestion of 
% Brian Harms <brianharms@hotmail.com> and Ivan Cao-Berg <icaoberg@cmu.edu> 
% (given at the user reviews area of Mathwork File exchage) The program should work 
% in the case of Becton EPics DLM FCS2.0, CyAn Summit FCS3.0 and FACSDiva type 
% list mode files.
%
% 29/01/2008 Updated to read the BD LSR II file format and including the comments of
% Allan Moser (Cira Discovery Sciences, Inc.)
%
% 24/01/2009 Updated to read the Partec CyFlow format file. Thanks for
% Gavin A Price
 
% if noarg was supplied
is16bitInflux=false;
if nargin == 0
     [FileName, FilePath] = uigetfile('*.*','Select fcs2.0 file');
     filename = [FilePath,FileName];
     if FileName == 0
          fcsdat = []; fcshdr = [];
          return;
     end
else
    filecheck = dir(filename);
    if size(filecheck,1) == 0
        msgWarning(['<html>Cannot find this file :' ...
            Html.FileTree(filename) '<hr></html>'], ...
            8, 'center', 'FCS reading issue..');
        fcsdat = []; fcshdr = [];
        return;
    end
end

% if filename arg. only contains PATH, set the default dir to this
% before issuing the uigetfile command. This is an option for the "fca"
% tool
[FilePath, FileNameMain, fext] = fileparts(filename);
FilePath = [FilePath filesep];
FileName = [FileNameMain, fext];
if  isempty(FileNameMain)
    currend_dir = cd;
    cd(FilePath);
    [FileName, FilePath] = uigetfile('*.*','Select FCS file');
     filename = [FilePath,FileName];
     if FileName == 0
          fcsdat = []; fcshdr = [];
          return;
     end
     cd(currend_dir);
end

%fid = fopen(filename,'r','ieee-be');
fid = fopen(filename,'r','b');
fcsheader_1stline   = fread(fid,64,'char');
fcsheader_type = char(fcsheader_1stline(1:6)');
%
%reading the header
%
if strcmp(fcsheader_type,'FCS1.0')
    msgbox('FCS 1.0 file type is not supported!','FCS reading issue','warn');
    fcsdat = []; fcshdr = [];
    fclose(fid);
    return;
elseif  strcmp(fcsheader_type,'FCS2.0') || strcmp(fcsheader_type,'FCS3.0')...
        || strcmp(fcsheader_type,'FCS3.1')% FCS2.0 or FCS3.0 types
    fcshdr.fcstype = fcsheader_type;
    FcsHeaderStartPos   = str2num(char(fcsheader_1stline(11:18)'));
    FcsHeaderStopPos    = str2num(char(fcsheader_1stline(19:26)'));
    FcsDataStartPos     = str2num(char(fcsheader_1stline(27:34)'));
    status = fseek(fid,FcsHeaderStartPos,'bof');
    fcsheader_main = fread(fid,FcsHeaderStopPos-FcsHeaderStartPos+1,'char');%read the main header
    warning off MATLAB:nonIntegerTruncatedInConversionToChar;
    fcshdr.filename = FileName;
    fcshdr.filepath = FilePath;
    % "The first character of the primary TEXT segment contains the
    % delimiter" (FCS standard)
    if fcsheader_main(1) == 12
        mnemonic_separator = 'FF';
    elseif fcsheader_main(1) == 9
        mnemonic_separator = 'TAB'; %added by RLF August 2010
    elseif fcsheader_main(1) == 10
        mnemonic_separator = 'LF';
    else
        mnemonic_separator = char(fcsheader_main(1));
    end
    if mnemonic_separator == '@'% WinMDI
      %  hm = msgbox([FileName,': The file cannot be read (Unsupported FCS type: WinMDI histogram file)'],'FcAnalysis info','warn');
        warning([FileName,': The file cannot be read (Unsupported FCS type: WinMDI histogram file)'],'FCS reading issue...','warn');
        
        fcsdat = []; fcshdr = [];
        fclose(fid);
        return;
    end
    fcshdr.TotalEvents = str2num(get_mnemonic_value('$TOT',fcsheader_main, mnemonic_separator));
    if fcshdr.TotalEvents == 0
        fcsdat = 0;
        fcsdatscaled = 0;
        return
    end
    if nargin<12
        startAtEvent=1;
        if nargin<11
            doLogicle=true;
            if nargin<10
                isSideScatterLinear=false;
                if nargin<9
                    maxEvents=0;
                end
            end
        end
    end
    fcshdr.isSideScatterLinear=isSideScatterLinear;
    fcshdr.doLogicle=doLogicle;
    totalEvents=fcshdr.TotalEvents;
    if startAtEvent>1
        totEvents=fcshdr.TotalEvents-(startAtEvent-1);
        if maxEvents>0
            if maxEvents<totEvents
                fcshdr.TotalEvents=maxEvents;
            end
        else
            fcshdr.TotalEvents=totEvents;
        end
    else
        if maxEvents>0
            if maxEvents<fcshdr.TotalEvents
                fcshdr.TotalEvents=maxEvents;
            end
        end
    end
    fcshdr.NumOfPar = str2num(get_mnemonic_value('$PAR',fcsheader_main, mnemonic_separator));
    fcshdr.Creator = get_mnemonic_value('CREATOR',fcsheader_main, mnemonic_separator);
    [fcshdr.CompLabels, fcshdr.CompMat, fcshdr.SPILL]=getKwTable(...
        fcsheader_main, mnemonic_separator, 'SPILL', fid);
    
    %%%%%%%%%%%%
    %%%%%%added by RLF to account for large files
    if FcsDataStartPos == 0
        FcsDataStartPos = str2num(get_mnemonic_value('$BEGINDATA',fcsheader_main, mnemonic_separator));
    end
    %%%%%%%%%%%%%%%%%%%%%
    plate = get_mnemonic_value('PLATE NAME',fcsheader_main,mnemonic_separator);
    if ~isempty(plate)
        fcshdr.plate=plate;
    end

        
    for i=1:fcshdr.NumOfPar
        fcshdr.par(i).name = get_mnemonic_value(['$P',num2str(i),'N'],fcsheader_main, mnemonic_separator);
        fcshdr.par(i).name2 = get_mnemonic_value(['$P',num2str(i),'S'],fcsheader_main, mnemonic_separator);
        fcshdr.par(i).range = str2num(get_mnemonic_value(['$P',num2str(i),'R'],fcsheader_main, mnemonic_separator));
        fcshdr.par(i).bit = str2num(get_mnemonic_value(['$P',num2str(i),'B'],fcsheader_main, mnemonic_separator));
%==============   Changed way that amplification type is treated ---  ARM  ==================
        par_exponent_str= (get_mnemonic_value(['$P',num2str(i),'E'],fcsheader_main, mnemonic_separator));
        if isempty(par_exponent_str)
            % There is no "$PiE" mnemonic in the Lysys format
            % in that case the PiDISPLAY mnem. shows the LOG or LIN definition
            islogpar = get_mnemonic_value(['P',num2str(i),'DISPLAY'],fcsheader_main, mnemonic_separator);
            if length(islogpar)==3 && isequal(islogpar, 'LOG')  % islogpar == 'LOG'
               par_exponent_str = '5,1'; 
            else % islogpar = LIN case
                par_exponent_str = '0,0';
            end
        end
   
       
        par_exponent= str2num(par_exponent_str);
        fcshdr.par(i).decade = par_exponent(1);
        if fcshdr.par(i).decade == 0
            fcshdr.par(i).log = 0;
            fcshdr.par(i).logzero = 0;
        else
            fcshdr.par(i).log = 1;
            if (par_exponent(2) == 0)
              fcshdr.par(i).logzero = 1;
            else
              fcshdr.par(i).logzero = par_exponent(2);
            end
        end

        gain_str = get_mnemonic_value(['$P',num2str(i),'G'],fcsheader_main, mnemonic_separator); % added by RLF 
        if ~isempty(gain_str) % added by RLF 
            fcshdr.par(i).gain=str2double(gain_str);
        else
            fcshdr.par(i).gain=1;
        end
        
%============================================================================================
    end
    fcshdr.starttime = get_mnemonic_value('$BTIM',fcsheader_main, mnemonic_separator);
    fcshdr.stoptime = get_mnemonic_value('$ETIM',fcsheader_main, mnemonic_separator);
    fcshdr.cytometry = get_mnemonic_value('$CYT',fcsheader_main, mnemonic_separator);
    fcshdr.date = get_mnemonic_value('$DATE',fcsheader_main, mnemonic_separator);
    fcshdr.byteorder = get_mnemonic_value('$BYTEORD',fcsheader_main, mnemonic_separator);
    fcshdr.datatype = get_mnemonic_value('$DATATYPE',fcsheader_main, mnemonic_separator);
    fcshdr.system = get_mnemonic_value('$SYS',fcsheader_main, mnemonic_separator);
    fcshdr.project = get_mnemonic_value('$PROJ',fcsheader_main, mnemonic_separator);
    fcshdr.experiment = get_mnemonic_value('$EXP',fcsheader_main, mnemonic_separator);
    fcshdr.cells = get_mnemonic_value('$Cells',fcsheader_main, mnemonic_separator);
    fcshdr.creator = get_mnemonic_value('CREATOR',fcsheader_main, mnemonic_separator);
    fcshdr.cytofDataShift = get_mnemonic_value('CYTOF_DATA_SHIFT',fcsheader_main, mnemonic_separator);
    fcshdr.isCytof=String.Contains(fcshdr.cytometry,'CYTOF');
    fcshdr.fcsheader_type=fcsheader_type;
    if startAtEvent>1
        errMsg=ShrinkFcs.GetErrMsg(fcshdr);
        if ~isempty(errMsg)
            error(errMsg);
        end
    end
    [fcshdr.numCols, fcshdr.fullColNames, fcshdr.compensatableColIdxs, ...
        fcshdr.logicleColIdxs, fcshdr.channelColNames, ...
        fcshdr.markerColNames, fcshdr.isLogicleCol, fcshdr.hasMarker, ...
        fcshdr.isNonReagentCol, fcshdr.numDataParameters, fcshdr.channelTypes, ...
        fcshdr.flowJoCompensatedChannels]=parseColInfo(fcshdr);
    try
        fcshdr.keywords=MatLabMap;
        if nargin>2 && ~isempty(fcshdr.keywords)
            fcshdr=copyNumbers(fcshdr, numericKeys, numericDefaults,...
                fcsheader_main, mnemonic_separator);
            fcshdr=copyBoolean(fcshdr, booleanKeys, booleanDefaults, ...
                fcsheader_main, mnemonic_separator);
            if nargin>6 && ~isempty(stringKeys)
                copyKeywords(fcshdr.keywords, stringKeys, stringDefaults, ...
                    fcsheader_main, mnemonic_separator);
            end
        end
    catch
        fcshdr.keywords=[];
    end

else
    %hm = msgbox([FileName,': The file cannot be read (Unsupported FCS type)'],'FcAnalysis info','warn');
    warning([FileName,': The file cannot be read (Unsupported FCS type)'],'FCS reading issue','warn');
    fcsdat = []; fcshdr = [];
    fclose(fid);
    return;
end
if headerOnly
    fcsdat = []; 
    fclose(fid);
    return;
end
%
%reading the events
%
fcsdat = []; 
fcshdr.isCytof32=fcshdr.isCytof && length(fcshdr.par)>1 && fcshdr.par(1).bit==32;
if strcmp(fcshdr.byteorder, '1,2,3,4')
    machineformat = 'ieee-le';
elseif strcmp(fcshdr.byteorder, '4,3,2,1')
    machineformat = 'ieee-be';
end
if startAtEvent>1
    skip=startAtEvent-1;
    jump=fcshdr.NumOfPar*skip*4;
    FcsDataStartPos=FcsDataStartPos+jump;
    fseek(fid, 0,'eof');
    eof=ftell(fid);
    if FcsDataStartPos>eof
        fcsdat=[];
        fcsdatscaled=[];
        fclose(fid);
        return;
    end
    status = fseek(fid,FcsDataStartPos,'bof');
    if skip+fcshdr.TotalEvents>totalEvents
        fcshdr.TotalEvents=totalEvents-skip;
    end
else
    status = fseek(fid,FcsDataStartPos,'bof');
end
if status ~= 0
    fcsdat=[];
    fcsdatscaled=[];
    fclose(fid);
    return;
end
if strcmp(fcsheader_type,'FCS2.0')
    if strcmp(mnemonic_separator,'\') || strcmp(mnemonic_separator,'FF')... %ordinary or FacsDIVA FCS2.0 
           || strcmp(mnemonic_separator,'/') || strcmp(mnemonic_separator,'TAB')% added by GAP 1/22/09 %added by RLF 09/02/10
        if fcshdr.par(1).bit == 16
            fcsdat = uint16(fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16')');
            if strcmp(fcshdr.byteorder,'1,2')...% this is the Cytomics data
                    || strcmp(fcshdr.byteorder, '1,2,3,4') %added by GAP 1/22/09
                fcsdat = bitor(bitshift(fcsdat,-8),bitshift(fcsdat,8));
            end
        elseif fcshdr.par(1).bit == 32
                if fcshdr.datatype ~= 'F'
                    fcsdat = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint32')');
                else % 'LYSYS' case
                    fcsdat = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32')');
                end
        else 
            bittype = ['ubit',num2str(fcshdr.par(1).bit)];
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],bittype, 'ieee-le')';
        end
    elseif strcmp(mnemonic_separator,'!')% Becton EPics DLM FCS2.0
        fcsdat_ = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16', 'ieee-le')';
        fcsdat = zeros(fcshdr.TotalEvents,fcshdr.NumOfPar);
        for i=1:fcshdr.NumOfPar
            bintmp = dec2bin(fcsdat_(:,i));
            fcsdat(:,i) = bin2dec(bintmp(:,7:16)); % only the first 10bit is valid for the parameter  
        end
    end
    fclose(fid);
elseif strcmp(fcsheader_type,'FCS3.0')|| strcmp(fcsheader_type,'FCS3.1')
    is16bitInflux=fcshdr.par(1).bit==16&& strcmp(fcshdr.datatype,'I') && ...
        String.ContainsI(fcshdr.cytometry, 'Influx');
    if is16bitInflux
        fcsdat = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16',machineformat)');
        for  i = 1 : fcshdr.NumOfPar
            if fcshdr.par(i).log
                decades = fcshdr.par(i).decade;
                range = fcshdr.par(i).range;
                antiLog = fcshdr.par(i).logzero;
                fcsdat(:,i) = antiLog*10.^(decades * fcsdat(:,i)/range);
            end
        end
    elseif ~fcshdr.isCytof32 && strcmp(mnemonic_separator,'|') ...
            && strcmp(fcshdr.datatype,'I')  % CyAn Summit FCS3.0
        fcsdat_ = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16','ieee-le')');
        fcsdat = zeros(size(fcsdat_));
        new_xrange = 1024;
        for i=1:fcshdr.NumOfPar
            fcsdat(:,i) = fcsdat_(:,i)*new_xrange/fcshdr.par(i).range;
            fcshdr.par(i).range = new_xrange;
        end
    elseif strcmp(mnemonic_separator,'/')
        if contains(lower(fcshdr.cytometry),'accuri')  % Accuri C6, this condition added by Rob Egbert, University of Washington 9/17/2010
            fcsdat = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'int32',machineformat)');
        elseif contains(lower(fcshdr.cytometry),'partec')%this block added by GAP 6/1/09 for Partec, copy/paste from above
            fcsdat = uint16(fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16',machineformat)');
            %fcsdat = bitor(bitshift(fcsdat,-8),bitshift(fcsdat,8));
        elseif contains(lower(fcshdr.cytometry),'lx') % Luminex data
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'int32',machineformat)';
            fcsdat = mod(fcsdat,1024);
        end
    %else % ordinary FCS 3.0
    end
    if isempty(fcsdat)
        %%%%%edited by RLF 06_30_10
        if fcshdr.datatype=='F'
            btype='float32';
        elseif fcshdr.datatype=='D'
            btype='double';
        elseif fcshdr.datatype=='I'
            btype='int32';
        end
        fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],...
            btype,machineformat)';      
    end
    fclose(fid);
end
if nargout==3 
    if is16bitInflux
        fcsdatscaled=fcsdat;
    else
        %
        %calculate the scaled events (for log scales)
        %
        fcsdatscaled = zeros(size(fcsdat));
        for  i = 1 : fcshdr.NumOfPar
            if ~fcshdr.par(i).log
                fcsdatscaled(:,i)  = fcsdat(:,i);
            else
                decades = fcshdr.par(i).decade;
                range = fcshdr.par(i).range;
                antiLog = fcshdr.par(i).logzero;
                fcsdatscaled(:,i) = antiLog*10.^(double(fcsdat(:,i))/range*decades);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function to = copyNumbers(to, keys, defaults, hdr, mnemonic_separator)
N=length(keys);
for i=1:N
    name=keys{i};
    value=get_mnemonic_value(name, hdr, mnemonic_separator);
    if isempty(value)
        to.(name) = defaults(i);
    else
        numericValue=str2num(value);
        if isempty(numericValue)
            to.(name) = defaults(i);
        else
            to.(name) = numericValue;
        end
    end
end

function to=copyBoolean(to, keys, defaults, hdr, mnemonic_separator)
N=length(keys);
for i=1:N
    name=keys{i};
    value=get_mnemonic_value(name, hdr, mnemonic_separator);
    if isempty(value)
        to.(name) = defaults(i);
    else
        to.(name) = strcmp(value,'1');
    end
end    



function copyKeywords(to, keys, defaults, hdr, sep)
N=length(keys);
for i=1:N
    name=keys{i};
    value=copyKeyword(to, name, true, defaults{i}, hdr, sep);
    if isempty(value)
        value=get_mnemonic_value( ['$Count.' name], hdr, sep);
        if ~isempty(value)
            N=str2num(value);
            for j=1:N
                addKeyword(to, name, j, defaults{i}, hdr,sep);
            end
        else % pre cluster genie keyword like EXPERIMENT NAME
            value=copyKeyword(to,  name, false, defaults{i}, hdr, sep);
        end
    end
end

function value=copyKeyword(to, name, dollarPrefix, default, hdr, mnemonic_separator)
if dollarPrefix
    value=get_mnemonic_value(['$' name], hdr, mnemonic_separator);
else
    value=get_mnemonic_value(name, hdr, mnemonic_separator);
end
if isempty(value)
    if ~isempty(default)
        to.set(name, default);
        value=default;
    end
else
    to.set(name, value);
end

function addKeyword(to, name, num, default, hdr, mnemonic_separator)
key=['$' num2str(num) '.' name];
value=get_mnemonic_value(key, hdr, mnemonic_separator);
if isempty(value)
    to.add(name, default);
else
    to.add(name, value);
end

function mneval = get_mnemonic_value(mnemonic_name,fcsheader,mnemonic_separator)

if strcmp(mnemonic_separator,'\')  || strcmp(mnemonic_separator,'!') ...
        || strcmp(mnemonic_separator,'|') || strcmp(mnemonic_separator,'@')...
        || strcmp(mnemonic_separator, '/') % added by GAP 1/22/08
    mnemonic_startpos = findstr(char(fcsheader'), [mnemonic_name mnemonic_separator]);
    if isempty(mnemonic_startpos)
        mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
        if isempty(mnemonic_startpos)
            mneval = [];
            return;
        end
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length;
    next_slashes = findstr(char(fcsheader(mnemonic_stoppos+1:end)'),mnemonic_separator);
    next_slash = next_slashes(1) + mnemonic_stoppos;
    
    mneval = char(fcsheader(mnemonic_stoppos+1:next_slash-1)');
elseif strcmp(mnemonic_separator,'FF')
    mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
    next_formfeeds = find( fcsheader(mnemonic_stoppos+1:end) == 12);
    if isempty(next_formfeeds)
        mneval='';
    else
        next_formfeed = next_formfeeds(1) + mnemonic_stoppos;
        mneval = char(fcsheader(mnemonic_stoppos + 1 : next_formfeed-1)');
    end
elseif strcmp(mnemonic_separator,'TAB') %added by RLF August 2010
    mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    nn=length(mnemonic_startpos);
    if nn>1
        ln=length(mnemonic_name);
        found=false;
        for i=1:nn
            if fcsheader(mnemonic_startpos(i)+ln)==9
                found=true;
                mnemonic_startpos=mnemonic_startpos(i);
                break;
            end
        end
        if ~found
            mneval = [];
            return;
        end
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
    next_formfeeds = find( fcsheader(mnemonic_stoppos+1:end) == 9);
    next_formfeed = next_formfeeds(1) + mnemonic_stoppos;
    
    mneval = char(fcsheader(mnemonic_stoppos + 1 : next_formfeed-1)');
else  %added by SWM August 2014
    mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
    next_formfeeds = find( fcsheader(mnemonic_stoppos+1:end) == mnemonic_separator);
    next_formfeed = next_formfeeds(1) + mnemonic_stoppos;
    
    mneval = char(fcsheader(mnemonic_stoppos + 1 : next_formfeed-1)');
end

function [header, table, comp]= getKwTable(...
    fcsheader_main, mnemonic_separator, Keyword, fid)
comp = get_mnemonic_value(Keyword,fcsheader_main,mnemonic_separator);
header=[];
table=[];
try
if ~isempty(comp)
    commas=strfind(comp,',');
    if isempty(commas)
        return;
    end
    nc=str2double(comp(1:commas(1)-1));  %first element is number of comped columns
    header=cell(1,nc);
    for i=1:nc
        header{i}=comp(commas(i)+1:commas(i+1)-1); %list of labels of comped columns
    end
    table=zeros(1,nc^2);
    for i=nc+1:length(commas)-1
        table(i-nc)=str2double(comp(commas(i)+1:commas(i+1)));
    end
    table(end)=str2double(comp((commas(end)+1):end));
    table=reshape(table,[nc nc]);
    for i=1:nc
        table(i,i)=1;
    end
    table=table';
end
catch
    comp = [];
    header=[];
    table=[];
    %fclose(fid);
    %throw ex;
end

