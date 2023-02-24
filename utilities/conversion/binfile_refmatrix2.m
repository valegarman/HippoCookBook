function [chMat,posMat,recDurSamps,RBDLen] = binfile_refmatrix(channels, RBDLen,recDurSamps,doPos,f_from)

% this creates a reference matrix for extraction of channel's contiguous samples in read_raw_bin_wav.m
% written by M Rutledge, 2/9/2011

remap_channels=[32, 33, 34, 35, 36, 37, 38, 39,0, 1, 2, 3, 4, 5, 6, 7,40, 41, 42, 43, 44, 45, 46, 47,8, 9, 10, 11, 12, 13, 14, 15,48, 49, 50, 51, 52, 53, 54, 55,16, 17, 18, 19, 20, 21, 22, 23,56, 57, 58, 59, 60, 61, 62, 63,24, 25, 26, 27, 28, 29, 30, 31 ];

chMat(1:length(channels),1:recDurSamps)=NaN; 
for channel=channels
    channeled=remap_channels(channel); 
    for i=1:3
        %   16 is packet header. 64 is gap between same channel consecutive
        %   samples within a packet. 216 is packet length
        chMat(channel,i:3:recDurSamps)=16+1+channeled+64*(i-1):216:RBDLen; 
    end          
end

wherenow = ftell(f_from);
frewind(f_from);  a = textscan(f_from, '%[ADU2]',4); 
fseek(f_from,wherenow,'bof');
if strcmp(a{1,1},'ADU2') & doPos==1 % if position data there and position analysis wanted
    for i=7:16  %pos data locations within each 432 byte (216 RBD elements) packet according to DacqUSBFileFormats doc.
        posMat(i-6:10:10*numPacks)=(i:216:RBDLen);
    end
else
    posMat=[];
end