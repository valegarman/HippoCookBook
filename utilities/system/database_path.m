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
        case 'DESKTOP-1SBCILE'
            dp = 'Z:\';
<<<<<<< HEAD
        case 'IMW02691'
            dp = 'Z:\';
=======
        case 'MountainJorge'
            dp = [];
>>>>>>> 917727d382d59f19d76c839fcc9db74328ca608b
        otherwise
            error('Computer name not found! Not possible to retrieve database path!');
    end

end