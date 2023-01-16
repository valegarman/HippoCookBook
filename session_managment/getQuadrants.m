function [quadrants] = getQuadrants(varargin)
% Computes the quadrants that are on and the HM response (good or bad
% trial)
% INPUTS
%
%
%
%
%
%
%
%
% OUTPUT
%   'quadrants'
%
% Pablo Abad 2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. position_seq (sequence of flashes : '1: up-left', '2: bottom-right', '3: bottom-left', '4: up-right')
% 2. answers (sequence of answers : '7: Up-left', '3: Bottom-right', '1 : Bottom-left', '9: Up-right') 
% 3. time_to_respond


% Parse options
p = inputParser;

addParameter(p,'basepath',pwd,@istruct);
addParameter(p,'sr',30000,@isnumeric);

addParameter(p,'position_BL',3,@isnumeric);
addParameter(p,'position_BR',2,@isnumeric);
addParameter(p,'position_UL',1,@isnumeric);
addParameter(p,'position_UR',4,@isnumeric);

addParameter(p,'answer_BL',1,@isnumeric);
addParameter(p,'answer_BR',3,@isnumeric);
addParameter(p,'answer_UL',7,@isnumeric);
addParameter(p,'answer_UR',9,@isnumeric);

parse(p,varargin{:});

basepath = p.Results.basepath;
sr = p.Results.sr;

position_BL = p.Results.position_BL; 
position_BR = p.Results.position_BR;
position_UL = p.Results.position_UL;
position_UR = p.Results.position_UR;

answer_BL = p.Results.answer_BL;
answer_BR = p.Results.answer_BR;
answer_UL = p.Results.answer_UL;
answer_UR = p.Results.answer_UR;

%% Get quadrants and correct/incorrect trials

try
    file = dir('*pulses.events.mat');
    load(file.name);
    
    file = dir('data_quadrants.mat');
    load(file.name);
catch
    warning('Not possible to getQuadrants...');
end

% Correct and incorrect trials

quadrants = [];

count = 1;
for ii = 1:length(pulses.timestampsOn)
   if ~isempty(pulses.timestampsOn{ii})
       quadrants.ts = pulses.timestampsOn{ii};
       count = count + 1;
   end
end

quadrants.ts(1) = [];
quadrants.ts(end) = [];
quadrants.answer = answers;
quadrants.position = position_seq;

for ii = 1:length(quadrants.position)
    if quadrants.position(ii) == position_BL && quadrants.answer(ii) == answer_BL % '3': Botom-left / '1': Bottom-left
        quadrants.choice(ii) = 1; % correct;
    elseif quadrants.position(ii) == position_BL && quadrants.answer(ii) ~= answer_BL
        quadrants.choice(ii) = 0; % incorrect
    end
    
    if quadrants.position(ii) == position_BR && quadrants.answer(ii) == answer_BR % '2': Bottom-right / '3': Bottom-right
        quadrants.choice(ii) = 1; % correct
    elseif quadrants.position(ii) == position_BR && quadrants.answer(ii) ~= answer_BR
        quadrants.choice(ii) = 0; % incorrect
    end
    
    if quadrants.position(ii) == position_UL && quadrants.answer(ii) == answer_UL % '1': Up-left / '7': Up-left
        quadrants.choice(ii) = 1; % correct
    elseif quadrants.position(ii) == position_UL && quadrants.answer(ii) ~= answer_UL
        quadrants.choice(ii) = 0; % incorrect
    end
    
    if quadrants.position(ii) == position_UR && quadrants.answer(ii) == answer_UR % '4': Up-right / '9': Up-right
        quadrants.choice(ii) = 1; % correct
    elseif quadrants.position(ii) == position_UR && quadrants.answer(ii) ~= answer_UR
        quadrants.choice(ii) = 0; % incorrect
    end
end

quadrants.performance = sum(quadrants.choice) / length(quadrants.choice) * 100;

[~,fbasename,~] = fileparts(pwd);
quadrants.folder = fbasename;


end