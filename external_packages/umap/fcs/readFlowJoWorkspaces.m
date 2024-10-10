function wsp=readFlowJoWorkspaces(file)
if nargin<1
    file=File.Home('.autoGate', 'application.properties');
end
s=File.ReadTextFile(file);
lines=strsplit(s, '\n');
N=length(lines);
wsp={};
p1=FlowJoTree.PROP_OPEN;
p2=[p1 '.count'];
for i=1:N
    if startsWith(lines{i}, p1) && ~startsWith(lines{i}, p2)
        line=lines{i};
        idx=String.IndexOf(lines{i}, '=');
        if idx>0
            wsp{end+1}=strtrim(line(idx+1:end));
        end
    end
end
end