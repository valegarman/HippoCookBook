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
        case 'DESKTOP-BEPJ8P0'
            dp = 'Z:\';
        case 'IMW02691'
            dp = 'Z:\';
        case 'IMW02703'
            dp = 'Z:\';
        case 'MountainJorge'
            dp = [];
        otherwise
            error('Computer name not found! Not possible to retrieve database path!');
    end

end