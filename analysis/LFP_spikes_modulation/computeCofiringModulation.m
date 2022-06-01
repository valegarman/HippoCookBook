
function cofiringMod = computeCofiringModulation(varargin)
% cofiringMod = computeCofiringModulation(varargin)
% 
% Plot and returns cofiring modulation for two phase modulation structures
% (tipically theta during REM and Run, for describing REM shifing activity)
% as in Mizuseki et al, 2011.
% INPUT
% <optional>
% phasemod1     phase modulation structure in X, as generated with
%                   PhaseModulation.mat. If empty, loads thetaREMMod
% phasemod2     phase modulation structure in Y, as generated with
%                   PhaseModulation.mat If empty, loads thetaRunMod
% basepath      Default, pwd
% saveMat       Default, true
% plotOpt       Default, true
% p_valueCutoff Rayleigh test p value for discarting unmodulated cells
%                   (tipically .05).
% shiftingThreshold   
%               Distance, in radiants, for defining shifting behaviour 
%               (default pi/2 rad, 90 deg)
%
% OUTPUT
% cofiringMod
%
%%  Manuel Valero 2022
%% Dealing with inputs
p = inputParser;

addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'phasemod1',[]);
addParameter(p,'phasemod2',[]);
addParameter(p,'plotOpt',true, @islogical);
addParameter(p,'saveMat',true, @islogical);
addParameter(p,'p_valueCutoff',0.05, @isnumeric);
addParameter(p,'shiftingThreshold',pi/2, @isnumeric);

parse(p,varargin{:})

basepath = p.Results.basepath;
phasemod1 = p.Results.phasemod1;
phasemod2 = p.Results.phasemod2;
plotOpt = p.Results.plotOpt;
saveMat = p.Results.saveMat;
p_valueCutoff = p.Results.p_valueCutoff;
shiftingThreshold = p.Results.shiftingThreshold;

% dealing with inputs 
prevPath = pwd;
cd(basepath);

if isempty(phasemod1)
    targetFile = dir('*.thetaREM_6-12.PhaseLockingData.cellinfo.mat'); 
    phasemod1 = importdata(targetFile.name);
    phasemod1.epochsName = 'thetaREM-6-12';
else
    phasemod1.epochsName = 'phasemod1';
end
if isempty(phasemod2)
    targetFile = dir('*.thetaRun_6-12.PhaseLockingData.cellinfo.mat');
    phasemod2 = importdata(targetFile.name);
    phasemod2.epochsName = 'thetaRun-6-12';
else
    phasemod2.epochsName = 'phasemod2';
end

significant_cells = (phasemod1.phasestats.p < p_valueCutoff) & ...
    (phasemod2.phasestats.p < p_valueCutoff);

phase_shift = wrapToPi(circ_dist(phasemod1.phasestats.m, phasemod2.phasestats.m));
phase_shift(~significant_cells) = NaN;
phase_shift_normalized = abs(phase_shift)/pi;
shiftingUnits = abs(phase_shift) > shiftingThreshold;

cofiringMod.phasemod1 = phasemod1;
cofiringMod.phasemod2 = phasemod2;
cofiringMod.m_phase1 = phasemod1.phasestats.m;
cofiringMod.m_phase2 = phasemod2.phasestats.m;
cofiringMod.significant_units = significant_cells;
cofiringMod.phase_shift = phase_shift;
cofiringMod.phase_shift_normalized = phase_shift_normalized;
cofiringMod.shiftingUnits = shiftingUnits;
cofiringMod.shiftingThreshold = shiftingThreshold;
cofiringMod.p_valueCutoff = p_valueCutoff;

if saveMat
    disp('Saving results...');
    save([basenameFromBasepath(pwd) '.cofiringMod_', phasemod1.epochsName, '_' phasemod1.epochsName '.cellinfo.mat'],'cofiringMod');
end

if plotOpt
    figure
    hold on
    units = ~cofiringMod.significant_units;
    s1 = scatter([cofiringMod.phasemod1.phasestats.m(units) cofiringMod.phasemod1.phasestats.m(units)...
        cofiringMod.phasemod1.phasestats.m(units)+2*pi cofiringMod.phasemod1.phasestats.m(units)+2*pi],...
        [cofiringMod.phasemod2.phasestats.m(units) cofiringMod.phasemod2.phasestats.m(units)+2*pi...
        cofiringMod.phasemod2.phasestats.m(units)+2*pi cofiringMod.phasemod2.phasestats.m(units)],30,[.9 .9 .9],"filled");
    % no shifting
    units = cofiringMod.significant_units & ~cofiringMod.shiftingUnits;
    s2 = scatter([cofiringMod.phasemod1.phasestats.m(units) cofiringMod.phasemod1.phasestats.m(units)...
        cofiringMod.phasemod1.phasestats.m(units)+2*pi cofiringMod.phasemod1.phasestats.m(units)+2*pi],...
        [cofiringMod.phasemod2.phasestats.m(units) cofiringMod.phasemod2.phasestats.m(units)+2*pi...
        cofiringMod.phasemod2.phasestats.m(units)+2*pi cofiringMod.phasemod2.phasestats.m(units)],50,[.3 .3 .3],"filled");
    % shifting
    units = cofiringMod.significant_units & cofiringMod.shiftingUnits;
    s3 = scatter([cofiringMod.phasemod1.phasestats.m(units) cofiringMod.phasemod1.phasestats.m(units)...
        cofiringMod.phasemod1.phasestats.m(units)+2*pi cofiringMod.phasemod1.phasestats.m(units)+2*pi],...
        [cofiringMod.phasemod2.phasestats.m(units) cofiringMod.phasemod2.phasestats.m(units)+2*pi...
        cofiringMod.phasemod2.phasestats.m(units)+2*pi cofiringMod.phasemod2.phasestats.m(units)],50,[.8 .3 .7],"filled");
    xlim([0 4*pi]); ylim([0 4*pi]);
    x_wave = 0:0.1:4*pi; y_wave = (cos(x_wave)+1)/2;
    plot(x_wave, y_wave,'-','color',[.7 .7 .7]);
    plot(y_wave, x_wave,'-','color',[.7 .7 .7]);
    set(gca,'TickDir','out','XTick',[0:pi:4*pi],'XTickLabel',{'0','\pi','2\pi', '3\pi','4\pi'},...
        'YTick',[0:pi:4*pi],'YTickLabel',{'0','\pi','2\pi', '3\pi','4\pi'});
    xlabel([cofiringMod.phasemod1.epochsName ' (rad)']);
    ylabel([cofiringMod.phasemod2.epochsName ' (rad)']);
    
    mkdir('SummaryFigures');
    saveas(gcf,['SummaryFigures' filesep 'cofiringMod_', phasemod1.epochsName, '_' phasemod1.epochsName '.png']);

end


cd(prevPath);
end