function dp = nyuShare_path
% return dropbox folder path
    switch getenv('computername')
        case 'MANULAPTOP'
            dp = 'W:';
        case 'MANUPC'
            dp = 'W:';
        case 'MANUXPS'
            dp = 'W:';
    end

end