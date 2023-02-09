function convert_raw_bin1(filenames,channels,doPos)
% converts the waveform for selected channels from Dacq raw recording BIN
% files to individual channel wav stream BIN files.
% eg. a = 'c:\BIN processing\110729_pp1_raw.bin; 'c:\BIN processing\110729_pp2_raw.bin'
% convert_raw_bin({a,[1:16],0);
%doPos either 1 or 0 for includin position data conversion
%
% see DacqUSBFileFormats document form Axona downloads page for file structure
% adapted from read_raw_bin by M Rutledge 3/9/2011
    

for file=1:length(filenames) 
    
    mkdir(filenames{file}(1:end-4)); 
    f_from=fopen(filenames{file},'r');
    fseek(f_from,216*2,'bof'); %skip first 432 byte packet (it is odd.. maybe header information?)
    [userview, ~] = memory;
    arraySize = 0.4*userview.MaxPossibleArrayBytes/16;  % assuming (size(chMat)==size(RawBinaryData)) == same number of bytes ... ***   MAY NOT BE RIGHT.... ***
    arraySize =  arraySize - (rem(arraySize,216)); % avoiding splitting 432 byte packets
   
    % In the BIN file, each 48kbits/second sample is spread accross two 8byte packets of
    % 2's complement binary format. This captures them as single integers
    RawBinaryData = fread(f_from,arraySize,'2*int16=>int16'); % disp(['f_from is at ',num2str(ftell(f_from))]);
    RBDLen=length(RawBinaryData);
    numPacks=(RBDLen/216);
    recDurSamps=3*numPacks;
    
    count=0;
    [chMat,posMat,recDurSamps,RBDLen] = binfile_refmatrix2(channels,RBDLen,recDurSamps,doPos,f_from);
    for i=channels  %  (should add one for doPos)
        ichannels_fopen1; % creates the fopen file handles; one file per individual channell
    end
    tic; 
    cdnow=cd;    cd(filenames{file}(1:end-4));  
    while feof(f_from)==0
        wavData = RawBinaryData(chMat(channels,1:recDurSamps));   j=1;
        for i=channels
            eval(['fwrite(f_to_',num2str(i),', wavData(j,:),''int16'');']);  j=j+1;
        end
        posData=RawBinaryData(posMat(1:length(posMat)));
        RawBinaryData = fread(f_from,arraySize,'2*int16=>int16');  
        RBDLen=length(RawBinaryData);
        numPacks=(RBDLen/216);
        recDurSamps=3*numPacks;
        count=count+1;%disp(count);
    end    
    fclose('all'); disp(count);  toc    
    cd(cdnow);
end