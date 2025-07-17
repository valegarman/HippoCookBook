
function delete_amplifier_files(rootDir)
    % 1. Search for all 'amplifier.dat' files in rootDir and subdirectories
    files = dir(fullfile(rootDir, '**', 'amplifier.dat'));

    if isempty(files)
        fprintf('No files named "amplifier.dat" were found.\n');
        return;
    end

    % 2. List the files in the terminal
    fprintf('The following "amplifier.dat" files were found:\n');
    for i = 1:length(files)
        fprintf('%d: %s\n', i, fullfile(files(i).folder, files(i).name));
    end

    % 3. Ask for user confirmation
    warning('This action will make these files disappear permanently!');
    resp = input('\nDo you want to delete these files? (y/n): ', 's');
    if strcmpi(resp, 'y')
        for i = 1:length(files)
            fullPath = fullfile(files(i).folder, files(i).name);
            delete(fullPath);
            fprintf('Deleted: %s\n', fullPath);
        end
        fprintf('All files have been deleted.\n');
    else
        fprintf('No files were deleted.\n');
    end
end