function [led_pos,post,led_pix] = rawpos(posfile,n_leds)
%
% [posx,posy,post] = rawpos(posfile,colour)
%
% Read raw positions from Axona POS-file
% (C) 2004 Sturla Molden
% CBM NTNU
%

fid = fopen(posfile,'r','ieee-be');
for t = 1:11
   string = fgetl(fid);
end
edges = ones(1,4);
for t = 1:4;
   string = fgetl(fid);
   edges(t) = sscanf(string,'%*s %u');
end 
string = fgetl(fid);
timebase = sscanf(string,'%*s %u %*s');
for t = 1:10
	string = fgetl(fid);
end
conversion = sscanf(string,'%*s %u');
string = fgetl(fid);
no_samples = sscanf(string,'%*s %u');
fseek(fid,10,0);
temp = ones(no_samples,8);
post = ones(no_samples,1);
for ind = 1:no_samples
    post(ind) = fread(fid,1,'uint32'); 
    temp(ind,1) = fread(fid,1,'uint16');  
    temp(ind,2) = fread(fid,1,'uint16');
    temp(ind,3) = fread(fid,1,'uint16');
    temp(ind,4) = fread(fid,1,'uint16');
    temp(ind,5) = fread(fid,1,'uint16');
    temp(ind,6) = fread(fid,1,'uint16');
    temp(ind,7) = fread(fid,1,'uint16');
    temp(ind,8) = fread(fid,1,'uint16');
    ind2 = find( temp(ind,:) == 1023);
    temp(ind,ind2) = NaN;    
end
fclose(fid);
if n_leds == 1
    led_pos(:,1,1) = temp(:,1);
    led_pos(:,1,2) = temp(:,2);
    led_pix(:,1) = temp(:,5);
elseif n_leds == 2;
    led_pos(:,1,1) = temp(:,1);
    led_pos(:,1,2) = temp(:,2);
    led_pix(:,1) = temp(:,5);
    led_pos(:,2,1) = temp(:,3);
    led_pos(:,2,2) = temp(:,4);
    led_pix(:,2) = temp(:,6);
end
clear temp;
post = (post-post(1))/timebase;
% index = find(post>=0);
% post = post(index);
% led_pos = led_pos(index,:,:);
led_pos(:,:,1) = led_pos(:,:,1) + edges(1);
led_pos(:,:,2) = led_pos(:,:,2) + edges(3);
% [B,I] = unique(post); 
% post = post(I);
% led_pos = led_pos(I,:,:);
% led_pix = led_pix(I);
