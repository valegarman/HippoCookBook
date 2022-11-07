function [dataN] = Notched(data,fs,F1,F2,F3,F4)
% Notch filter for 4 different frequencies ( 50, 100 , 150 , 200 Hz)
% Uses a notch filter for 50, 100 , 150 or 200 Hz, depending on the number
% of inputs
% 
%   USAGE
%   dataN = bz_Notched(data,fs,F1,F2,F3,F4);
%
%   INPUTS
%   data    lfp.data
%   fs      sample Rate
%   F1      F1 to be Notched
%   F2      F2 to be Notched
%   F3      F3 to be Notched
%   F4      F4 to be Notched
%
%   OUPUT
%   dataN   data Notched

if nargin<3
    disp('error not enough arguments, check for FS etc...')
elseif nargin==3
%     CA1 =filter(b,a,CA1);
    [Hd b a] = NotchIrFir50 (fs);
    dataN=filter(b,a,data);
elseif nargin==4
    [Hd b a] = NotchIrFir50(fs) ;
    dataN=filter(b,a,data);

    [Hd b a] = NotchIrFir100(fs);
    dataN=filter(b,a,dataN);

elseif nargin == 5
    [Hd b a] = NotchIrFir50(fs) ;
    dataN=filter(b,a,data);

    [Hd b a] = NotchIrFir100(fs);
    dataN=filter(b,a,dataN);
    
    [Hd b a] = NotchIrFir150(fs);
    dataN=filter(b,a,dataN);
elseif nargin == 6
    [Hd b a] = NotchIrFir50(fs) ;
    dataN=filter(b,a,data);

    [Hd b a] = NotchIrFir100(fs);
    dataN=filter(b,a,dataN);
    
    [Hd b a] = NotchIrFir150(fs);
    dataN=filter(b,a,dataN);
    
    [Hd b a] = NotchIrFir200(fs);
    dataN = filter(b,a,dataN);
    
    
end