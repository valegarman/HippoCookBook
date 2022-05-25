
function eventDetector(varargin)
% Arduino interface to threshold detect online events.
%
% USAGE
%
% eventDetector(varargin)
%
% INPUTS
%
% port            -serial port (by default try '/dev/cu.usbmodem14101' and COM1-10 )
%
%% Parse options
p = inputParser;
addParameter(p,'port',[],@isstr);

parse(p,varargin{:});

port = p.Results.port;

% to do: serial port searching!!

s1 = serial('COM7');    % define serial port % 
s1.BaudRate=9600;               % define baud rate
%set(s1,'DataBits', 8);
%set(s1,'StopBits', 1);
fopen(s1);
set(s1, 'terminator', 'LF');    % define the terminator for println
s1.ReadAsyncMode = 'continuous';

try 
    %w=fscanf(s1,'%s');              % hand shanke with arduino
    %if (w=='A')
    %    disp('Arduino is there!');
    %    fprintf(s1,'%s\n','A');     % establishContact just wants 
    %                                % something in the buffer
    %end
    t = 0;
    v = 0;
    w = 0;
    global winSize minP maxP delay thr;                 % moving windows size
    winSize = 1;
	minP = -2.;                         % set y-min
    maxP = 2;  
    delay = 0.001;                  % make sure sample faster than resolution 
    thr = 1;
    zeroPlot = 0;
    %Set up Plot
    hMain = figure;
    set(hMain,'MenuBar','none','color','w','Name','eventDetector','NumberTitle','off');
    
    plotGraph = plot(t,v,'-r');  % every AnalogRead needs to be on its own Plotgraph
    set(gca,'TickDir','out');
    hold on                            %hold on makes sure all of the channels are plotted
    plotThr = plot([zeroPlot zeroPlot+winSize],[thr thr],'--g');  % every AnalogRead needs to be on its own Plotgraph
    xlabel('s'); ylabel('V');
    box off
    axis([zeroPlot zeroPlot+winSize minP maxP]);
    
    ii = 1;
    tic
    disp('Plotting data!');
    uicontrol('Style', 'pushbutton', 'String', 'Close',...
                        'Units','normalize','Position', [.89 .93 .10 .06],...
                        'Callback', 'close(gcbf)');
    changeXMenu = uicontrol('Style', 'popupmenu', 'String', {'.5s', '1s', '5s', '10s'},...
                        'Units','normalize','Position', [.76 .925 .125 .06],'Value',3,...
                        'Callback', @changeX);
    changeYMenu = uicontrol('Style', 'popupmenu', 'String', {'+.1V','+.5V', '+1V', '+5V',...
                        '+/-.1V','+/-.5V','+/-1V','+/-5V'},'Value',4,...
                        'Units','normalize','Position', [.63 .925 .135 .06],...
                        'Callback',@changeY);
    changeDelayMenu = uicontrol('Style', 'popupmenu', 'String', {'1KHz','10KHz', '20KHz', '30KHz'},...
                        'Value',2,...
                        'Units','normalize','Position', [.49 .925 .145 .06],...
                        'Callback',@changeDelay);
    uicontrol('Style', 'text', 'String', 'Thr: ',...
                        'Units','normalize','Position', [.36 .93 .065 .05],...
                        'BackgroundColor','w',...
                        'Callback', @changeThr);
    changeThrMenu = uicontrol('Style', 'edit', 'String', '1',...
                        'Units','normalize','Position', [.415 .93 .065 .06],...
                        'Callback', @changeThr);
    % fwrite(s1, 100, 'uint8');
    while ishandle(hMain) % Loop when Plot is Active will run until plot is closed!!
        absTime = toc; % track absolute time
        t(ii)=absTime; % time plot is absolute time - time plot
        v(ii)=fscanf(s1,'%f')/204.6;        % must define the input % d, %f, %s, etc,
        w(ii)=fscanf(s1,'%f')/204.6;        % must define the input % d, %f, %s, etc, 

        %set(plotGraph,'XData',t,'YData',v); % Converting 1023 base to V (1023 is 5V)
        axis([zeroPlot zeroPlot+winSize minP maxP]);
        % disp(v(ii));
        if t(ii) > zeroPlot+winSize % max winSize
            ii = 1; % initialize ii
            plotGraph = plot(t,v,'-r');
            t = t(ii);
            v = v(ii); 
            zeroPlot = absTime;
            %delete(plotGraph);
            %plotGraph = plot(t,v,'-r');  % start again plotting!
            %delete(plotThr);
            %plotThr = plot([zeroPlot zeroPlot+winSize],[thr thr],'--g');  % every AnalogRead needs to be on its own Plotgraph
            % fwrite(s1, 5, 'int8'); % send threshold
            clear t v w
        end
        ii = ii + 1;
        pause(delay);
    end
    disp('Closing!')
    fclose(s1);
    
catch exception
    disp('Error!!!');
    fclose(s1);                 % always, always want to close s1
    throw (exception);
end 

    function [] = changeX(hObject, eventdata) % callback for x axis
        switch get(changeXMenu,'Value')   
            case 1
                winSize = .5;
            case 2
                winSize = 1;
            case 4
                winSize = 10;
            otherwise
                winSize = 5;
        end
        delete(plotThr);
        plotThr = plot([zeroPlot zeroPlot+winSize],[thr thr],'--g');  % every AnalogRead needs to be on its own Plotgraph
    end

    function [] = changeY(hObject, eventdata) % callback for x axis
        switch get(changeYMenu,'Value')   
            case 1
                minP = 0; maxP = .11;
            case 2
                minP = 0; maxP = 0.51;
            case 3
                minP = 0; maxP = 1.1;
            case 5
                minP = -.11; maxP = .11;
            case 6
                minP = -.51; maxP = .51;
            case 7
                minP = -1.1; maxP = 1.1;
            case 8
                minP = -5.1; maxP = 5.1;
            otherwise
                minP = 0; maxP = 5.1;
        end
    end

    function [] = changeDelay(hObject, eventdata) % callback for x axis
        switch get(changeDelayMenu,'Value')   
            case 1
                delay = 0.001;
            case 3
                delay = 0.00005;
            case 4
                delay = 0.000033;
            otherwise
                delay = 0.0001;
        end
    end
    function [] = changeThr(hObject, eventdata) % callback for x axis
        thr = str2double(get(changeThrMenu,'String'));
        delete(plotThr);
        plotThr = plot([zeroPlot zeroPlot+winSize],[thr thr],'--g');  % every AnalogRead needs to be on its own Plotgraph
        fwrite(s1, 100, 'int8'); % send threshold
    end
end