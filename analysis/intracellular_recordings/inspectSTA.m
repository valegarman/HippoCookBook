
function inspectSTA(StaDat, u, varargin)
% StaDat = myMpSta(d, cellPos, u, fs)
%
% Perform Spike trigger average (STA) over intracell data
%
% INPUTS
% StaDat        N row structure with STA data, where N refeer to 
%   {N}.wave    Amplitude (mV) x time (ms) matrix with STA response
%   {N}.wave_df Driving force (mV, amplitude minus average) x time (ms)
%               matrix with STA response
%   {N}.curr    Current (mV) x time (ms) matrix during STA responses
%               This function uses .wave and .wave_df
% <options> Optional list of property-value pairs (see table below)
%   fs          (scalar) Sampling frequency in HZ (default 20000).
%   winSTA      2 x 1 vector with STA time window, default [-10 100]
%   sigPlot     specify what signal trace will be plot between 'raw',
%               'drivingForce' and 'filt' (5Hz hp filt, default)
%   zoomPlot    Time window on STA plots, as x-axlis limit. Default all
%               windows.
%
%   Manu Valero 2018
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default values
fs = 20000;
for ii = 1: size(StaDat{1},2)
    s(ii) = size(StaDat{1}(ii).dfilt,1);
end
winSTA = [0 s(find(s>0,1))/(fs/1000)];
subthr = 1;
zoomPlot = winSTA;

% Parse options
for ii = 1:2:length(varargin)
    if ~ischar(varargin{ii})
		error(['Parameter ' num2str(ii+3) ' is not a property (type ''help <a href="matlab:help myMpSta">myMpSta</a>'' for details).']);
    end
    switch(lower(varargin{ii}))
        case 'fs'
            fs = varargin{ii+1};
        case 'win'
            winSTA = varargin{ii+1};
            zoomPlot = winSTA;
        case 'zoomplot'
            zoomPlot = varargin{ii+1};
        case 'subthr'
            subthr = varargin{ii+1};
    end
end

disp('Inspect intracell STA...');
xt = gpuArray(linspace(winSTA(1), winSTA(2), s(find(s>0,1))));
figure
ii = 1;
while ii <= size(StaDat,2)
    jj = 1;
    while jj <= size(StaDat{ii},2)
            disp('Plotting... ');
            [xout,r] = autocorr_spikes(u.ts{jj},fs,26,1);
            clf;
            subplot(2,2,1)
            area(xout,r,'LineStyle','none');
            xlim([-25 25]); xlabel('ms'); ylabel('prob');
            text(.05,.9,strcat(num2str(numel(u.ts{jj})),' spks, ',...
                num2str(round(numel(u.ts{jj})/((u.tsC{jj}(end)-u.tsC{jj}(1))/fs),2)),'Hz'),'Units','Normalized');
            title(strcat('Intra ','{ }', num2str(ii),'{ }', 'of','{ }', num2str(size(StaDat,2)),...
                ', unit', '{ }', num2str(jj),'{ }', 'of','{ }', num2str(size(u,2))),'FontWeight','normal');

            if size(StaDat{ii}(jj).dfilt,2) > 10
                if subthr 
                    idx = find(StaDat{ii}(jj).subThre);
                else
                    idx = find(ones(size(StaDat{ii}(jj).subThre)));
                end
                if numel(idx) > 10
                    subplot(2,2,2)
                    plotFill1(xt,StaDat{ii}(jj).dfilt(:,idx));
                    xlim(zoomPlot);
                    xlabel('ms'); ylabel('mV');
                    text(.05,.9,strcat(num2str(numel(idx)),' stas'),'Units','Normalized');
                    title('STA','FontWeight','normal');
                    text(.05,.05,'STA','Units','Normalized','color',[.0 .75 .75]);
                    
                    if ~isempty(StaDat{ii}(jj).surr)
                        plot(xt, mean(StaDat{ii}(jj).surr,2),'color',[.8 .2 .2],'lineWidth',2);
                        plot(xt, mean(StaDat{ii}(jj).dfilt(:,idx),2)-mean(StaDat{ii}(jj).surr,2),...
                            'color',[.2 .2 .2],'lineWidth',2);
                        text(.3,.05,'Jit. surr','Units','Normalized','color',[.8 .2 .2]);
                        text(.65,.05,'STA-Surr','Units','Normalized','color',[.2 .2 .2]);
                    end

                    subplot(2,2,3)
                    [sortResting, idxResting] = sort(StaDat{ii}(jj).mp(idx));
                    temp = StaDat{ii}(jj).dfilt(:,idx);
                    imagesc(xt,1:length(idxResting),temp(:,idxResting)',[-2 2]);
                    xlim(zoomPlot);
                    set(gca,'YDir','normal');
                    colormap jet
                    xlabel('ms'); ylabel('STAs (#)');
                    title('All STAs','FontWeight','normal');

                    subplot(2,2,4)
                    plot(sortResting,1:length(idxResting)); ylim([1 length(idxResting)]);
                    title('All STAs membrane potential','FontWeight','normal');
                    xlabel('mV'); ylabel('STAs (#)');
                end
            end

            tog = input('Press ENTER or 6 to continue, 4 to go back: ');
            if isempty(tog) || tog == 6 
                jj = jj + 1;
            elseif tog == 4
                jj = jj - 1;
            else
                warning('Input not recognized. ');
            end
            
    end
    fprintf('\nEnd of cell #%i.\n',ii);
    toh = input('Press ENTER or 6 to continue, 4 to go back: ');
    if isempty(tog) || tog == 6 
        ii = ii + 1;
    elseif tog == 4
        ii = ii - 1;
    else
        warning('Input not recognized. ');
    end
end

end