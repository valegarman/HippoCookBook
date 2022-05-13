function dp = database_path
% return dropbox folder path
    switch getenv('computername')
        case 'MANULAPTOP'
            dp = 'Z:\data';
        case 'MANUPC'
            dp = 'X:\data';
        case 'MANUXPS'
            dp = 'W:';
        case 'SB13FLPC017'
            %dp = 'W';
            dp = 'Z:\data';
    end

end