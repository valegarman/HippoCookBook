%% BatchScript for POTometry project

sessions_to_analyze = {'X:\unindexedSubjects\wt10\wt10_251027_sess1','X:\unindexedSubjects\wt10\wt10_251028_sess2','X:\unindexedSubjects\wt10\wt10_251029_sess3',...
    'X:\unindexedSubjects\wt10\wt10_251030_sess4','X:\unindexedSubjects\wt10\wt10_251031_sess5','X:\unindexedSubjects\wt10\wt10_251103_sess6','X:\unindexedSubjects\wt10\wt10_251104_sess7',...
    'X:\unindexedSubjects\wt10\wt10_251105_sess8','X:\unindexedSubjects\wt10\wt10_251106_sess9','Z:\wt10\wt10_251112_sess12','X:\unindexedSubjects\wt10\wt10_251114_sess14','Z:\wt10\wt10_251117_sess15',...
    'Z:\wt10\wt10_251118_sess16','X:\unindexedSubjects\wt10\wt10_251120_sess17','Z:\wt10\wt10_251125_sess20\SummaryFigures','Z:\wt10\wt10_251202_sess25','Z:\wt10\wt10_251210_sess28','Z:\wt10\wt10_251212_sess29'};

for ii = 1:length(sessions_to_analyze)

    fprintf(' > %3.i/%3.i session \n',ii, length(sessions_to_analyze));
    cd(sessions_to_analyze{ii});

     try
        ripples_fiber = fiberPhotometryModulation_pablo([],'eventType','ripples','reload_fiber',true);
    catch
        warning('No fiber recording in this session...');
     end
end

