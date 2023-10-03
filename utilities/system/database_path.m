function dp = database_path
% return dropbox folder path
[~,computerName] = system('hostname');
computerName = (strtrim(computerName));
    switch computerName
    % switch getenv('computername')
        case 'MANULAPTOP'
            dp = 'W:\Buzsakilabspace\Datasets\ValeroM';
        case 'MANUPC'
            dp = 'W:\Buzsakilabspace\Datasets\ValeroM';
        case 'MANUXPS'
            dp = 'W:';
        case 'SB13FLPC017'
            %dp = 'W';
            dp = 'W:\buzsakilab\Buzsakilabspace\Datasets\ValeroM';
        case 'Manuels-MacBook-Pro.local'
            dp = '/Volumes/NEURAL';
        otherwise
            dp = [];
    end

end