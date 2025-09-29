session = loadSession();

fs = session.extracellular.sr;
eventTimes = round([digitalIn.timestampsOn{6}' * fs]);



fid = fopen([session.general.name,'.ttl.evt'],'w');
for i = 1:length(eventTimes)
    fprintf(fid,'%d\t%d\n',eventTimes(i),1); % 1 puede ser el código del evento
end
fclose(fid);