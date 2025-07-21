
function delete_phy_cache(baseDir)
% remove_phy_cache(baseDir)
% Recursively finds and deletes all folders named '.phy' inside baseDir,
% after listing them, issuing a warning, and asking for confirmation.
% MV 2025

    if nargin < 1
        error('Please provide a base directory.');
    end

    % Find all directories named '.phy'
    phyDirs = dir(fullfile(baseDir, '**', '.phy'));

    if isempty(phyDirs)
        fprintf('No .phy folders found in %s.\n', baseDir);
        return;
    end

    % List found folders
    fprintf('Found the following .phy folders:\n');
    for k = 1:length(phyDirs)
        fprintf('  %s\n', fullfile(phyDirs(k).folder, phyDirs(k).name));
    end

    % Warning message
    warning('This operation will permanently delete the above folders and their contents.');

    % Ask user for confirmation
    userInput = input('Do you want to proceed? (y/n): ', 's');
    if ~strcmpi(userInput, 'y')
        fprintf('Operation cancelled by user.\n');
        return;
    end

    % Proceed to delete
    for k = 1:length(phyDirs)
        folderPath = fullfile(phyDirs(k).folder, phyDirs(k).name);
        try
            rmdir(folderPath, 's');
            fprintf('Deleted: %s\n', folderPath);
        catch ME
            warning('Could not delete folder %s. Error: %s', folderPath, ME.message);
        end
    end
end