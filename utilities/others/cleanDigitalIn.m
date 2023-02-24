TMaze = false;
OpenField = true;
LinearTrack = false;

if TMaze
    start = 6;
elseif OpenField
    start = 3;
elseif LinearTrack
    start = 5;
else
    start = 1;
end

for ii = start:length(digitalIn.timestampsOn)
    digitalIn.timestampsOn{ii} = [];
    digitalIn.timestampsOff{ii} = [];
    digitalIn.ints{ii} = [];
    digitalIn.dur{ii} = [];
    digitalIn.intsPeriods{ii} = [];
end
path = cd;
path_ = strsplit(path,filesep);
filename = path_{end};

save([filename,'.DigitalIn.events.mat'],'digitalIn');