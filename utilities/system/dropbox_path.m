function dp = dropbox_path
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
        case 'MANULAPTOP'
            dp = 'C:\Users\valeg\Dropbox';
        case 'MANUPC'
            dp = 'D:\Dropbox';
        case 'MANUXPS'
            dp = 'C:\Users\valeg\Dropbox';
        case 'manu'
            dp = '/Users/manu/Dropbox';
        case 'DESKTOP-1SBCILE'
            dp = 'C:\Users\mvalero\Dropbox';
        case 'IMW02691'
            dp = 'D:\Dropbox';
    end

end