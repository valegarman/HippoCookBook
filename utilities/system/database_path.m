function dp = database_path
% return dropbox folder path
[~,computerName] = system('hostname');
    % switch computerName
    switch getenv('computername')
        case 'MANULAPTOP'
            dp = 'W:\Buzsakilabspace\Datasets\ValeroM';
        case 'MANUPC'
            dp = 'X:\data';
        case 'MANUXPS'
            dp = 'W:';
        case 'SB13FLPC017'
            %dp = 'W';
            dp = 'W:\buzsakilab\Buzsakilabspace\Datasets\ValeroM';
        otherwise
            dp = [];
    end

end