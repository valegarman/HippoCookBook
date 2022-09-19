
function channelOrderOut = channelMapGenerator(channelOrderIn, varargin)
% Manuel Valero 2022

if nargin == 1
    disp('Recording system/ probe?');
    disp(' 1/ 32ch probe from NeuroNexus to intan/open ephys');
    disp(' 1/ 64ch probe from NeuroNexus to intan/open ephys');
    option = input('Your option (numeric): ');
    
    switch option
        case 1
            headstage = '32ch_NeuroNexus2intan';
        case 2
            headstage = '64ch_NeuroNexus2intan';
    end
elseif nargin == 2
    headstage = varargin{1};
else
    error('Wrong number of inputs!');
end

if strcmpi(headstage,'32ch_NeuroNexus2intan')
    neuroscopeOrder = [29  2 27  4 26  5 25  6 24  7 28  3 30  1 23  8 19 12 11 21 10 20 13 22  9 18 14 31  0 17 15 16];
%    layoutOrder =     [ 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
end

if strcmpi(headstage,'64ch_NeuroNexus2intan')
    neuroscopeOrder = [49 48 51 50 53 52 55 54 57 56 59 58 61 60 63 62 32 33 34 36 37 38 40 41 42 44 45 46 47 43 39 35     29 25 21 17 16 19 18 20 23 22 24 27 26 28 31 30  0  1  2  3  4  5  6  7  8 60 10 11 12 13 14 15];
%     layoutOrder =     [ 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32     33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64];
end

if isnumeric(channelOrderIn) && ismatrix(channelOrderIn)
    channelOrderOut.chs_0_index = neuroscopeOrder(channelOrderIn)';
    channelOrderOut.chs_1_index = neuroscopeOrder(channelOrderIn)' + 1;
elseif iscell(channelOrderIn)
    for ii = 1:length(channelOrderIn)
        chs_0_index{ii} = neuroscopeOrder(channelOrderIn{ii})';
        chs_1_index{ii} = neuroscopeOrder(channelOrderIn{ii})' + 1;
    end
    channelOrderOut.chs_0_index = chs_0_index;
    channelOrderOut.chs_1_index = chs_1_index;
end


end




