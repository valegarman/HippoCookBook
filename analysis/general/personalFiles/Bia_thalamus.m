
clear all

HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions_thalamus.csv']); % the variable is called allSessions

neuroscope_yes = sessionsTable.Neuroscope; 

index_good = [];

for ii = 1 : length(neuroscope_yes) % we check one by one what is inside of each cells, to know if is yes or no
  if strcmp(neuroscope_yes{ii},'yes') % once we are inside, we check with this if , if the cells are cointining 'yes' 
      index_good =[index_good; ii];   % we save the index of the cells in which there is written 'yes'
  end
end

list_sessions_thalamus = sessionsTable.SessionName(index_good); % we save all the session names we need 

% in one session I look at cell metrics

cell_metrics = loadCellMetrics;

shankID = cell_metrics.shankID';
ChannelsID = cell_metrics.maxWaveformCh';
thalamus_channels = [18;20;44;42];

% start to save info of neurons we need
Candidate_neurons_shank = find(shankID == 3); %saving the ID of teh neurons that are in the shank taht i want to study 
Candidate_neurons_channels = ChannelsID(find(shankID == 3)); %saving the channels in which the neurons that I have just saved above are located 

channels_members = ismember(Candidate_neurons_channels,thalamus_channels); %we look for the channels that are part of the channel we are investigating

Putative_neurons = Candidate_neurons_shank(channels_members); % now we save ONLY the neurons that are part of thalamus




