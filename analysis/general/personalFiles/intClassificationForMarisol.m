
% interneuron classification for Marisol

% 1 extract sessions
cd([nyuShare_path '\Buzsakilabspace\Datasets\SoulaM\IED_ Manuscript']);
load('AD_AnimalsSess.mat');

% find ripples sessions that are good in AD MICE
ii=1;
for kk=1:length(fieldnames(allSessions))
 list=struct2cell(allSessions);
 if list{kk,1}.IED == 1 &&  list{kk,1}.Genotype==1   ;    %list{kk,1}.littermates==2 && list{kk,1}.littermates==3
 ADgoodsess{ii,1}= list{kk,1}.path; 
 ADmonths(ii,1)=list{kk,1}.months;
 AD_IDs(ii,1)=list{kk,1}.ID;
 ii=ii+1;
 else 
 end    
end

% keep only months 3 and 4 and 9 10 and 11 to compare. 
c = ismember(ADmonths,[2 3 4 9 10 11 12]);
indexes = find(c);
ADmonths=ADmonths(c);
e = ismember(ADmonths,[2 3 4]);
ADmonths(find(e))=1;
f = ismember(ADmonths,[9 10 11 12]);
ADmonths(find(f))=2;

ADgoodsess=ADgoodsess(c);
AD_IDs=AD_IDs(c);
% color 
ADmap = brewermap(max(round(ADmonths)),'PiYG'); 

% ADgoodsess, paths to good sessions
% ADmonths 1 is young,
basepaths = [];
for ii = 1:size(ADgoodsess,1)
    % load all cells
    basepaths{ii} = [pwd ADgoodsess{ii}];
end

cell_metrics = loadCellMetricsBatch('basepaths',basepaths);

bins = linspace(-500,500,1001);
bins_0_to_10 = bins >= 0 & bins <= 10;
bins_40_to_50 = bins >= 40 & bins <= 50;

cell_metrics.burstIndex_AFR2017 = max(cell_metrics.acg.wide(bins_0_to_10,:))-mean(cell_metrics.acg.wide(bins_40_to_50,:));

% load all sessions
is_interneurons = ismember(cell_metrics.putativeCellType, {'Narrow Interneuron','Wide Interneuron'});

figure
plot(cell_metrics.acg_refrac(is_interneurons), cell_metrics.burstIndex_Doublets(is_interneurons),'or')

figure
histogram(cell_metrics.acg_tau_rise(is_interneurons))