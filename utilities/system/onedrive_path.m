function dp = onedrive_path
% return dropbox folder path
    os = computer;
    switch computer
        case 'MACI64'
            %computerName = char(java.net.InetAddress.getLocalHost.getHostName);
            computerName = char(java.lang.System.getProperty('user.name'));
        case 'PCWIN64'
            computerName = getenv('computername');
    end
    
    dp = [];
    switch computerName
        case 'manu'
            dp = '/Users/manu/Dropbox';
        case 'IMW02691'
            dp = 'C:\Users\mvalero\OneDrive - imim.es\';
        case ''% pc de andrea
            dp = '';
        otherwise
            error('Not recornized computer!');
    end

end