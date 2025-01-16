function dp = onedrive_path
% return dropbox folder path
    switch computer
       case 'MACA64'
            %computerName = char(java.net.InetAddress.getLocalHost.getHostName);
            computerName = char(java.lang.System.getProperty('user.name'));
        case 'PCWIN64'
            computerName = getenv('computername');
    end
    
    dp = [];
    switch computerName
        case 'manu'
            dp = '/Users/manu/Library/CloudStorage/OneDrive-imim.es/';
        case 'IMW02691'
            dp = 'C:\Users\mvalero\OneDrive - imim.es\';
        case 'IMW02703' % pc de andrea
            dp = 'C:\Users\agallardo\OneDrive - imim.es\';
        case 'DESKTOP-BEPJ8P0'
            dp = 'C:\Users\mpicco\OneDrive - imim.es\';
        case 'IMW02838' % pc de Ane, la mejor
            dp = 'C:\Users\amartinez11\OneDrive - imim.es\';
        case 'DESKTOP-4MOGMGG' % pablo's pc
            dp = 'C:\Users\pabad\OneDrive - imim.es\';
        otherwise
            error('Not recornized computer!');
    end

end