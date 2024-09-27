function location = nas_path(input)
% return dropbox folder path
    switch computer
       case 'MACI64'
            %computerName = char(java.net.InetAddress.getLocalHost.getHostName);
            computerName = char(java.lang.System.getProperty('user.name'));
        case 'PCWIN64'
            computerName = getenv('computername');
    end
    
    location = [];
    switch computerName
        case 'manu'
            keyboard;
        case 'IMW02691'
            location = whereIsThisNas1(input);
        case 'IMW02703' % pc de andrea
            location = whereIsThisNas1(input);
        case 'DESKTOP-BEPJ8P0'
            location = whereIsThisNas1(input);
        otherwise
            location = whereIsThisNas1(input);
    end

end

function here = whereIsThisNas1(input)
    switch lower(input)
        case 'neural'
            here = 'Z:';
        case 'neural3'
            here = 'Y:';
    end
end