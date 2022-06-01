
function [ints] = IDXtoINT(IDX, timestamps)
    % as simple as that. MV 2022 :)
    % assumes that the first second goes from 0 to 1s
    if isrow(IDX)
        IDX = IDX';
    end
    
    if isrow(timestamps)
        timestamps = timestamps';
    end
    
    dt = mode(diff(timestamps));
    timestamps = [timestamps(1)-dt; timestamps];
    IDX = [0; IDX;0];

    ints = [timestamps(find(diff(IDX)==1)) timestamps(find(diff(IDX)==-1))];
end