function build_run_epp
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    eppPath=fileparts(mfilename('fullpath'));
    autoGatePath=fileparts(eppPath);
    utilPath=fullfile(autoGatePath, 'util');
    buildPath=fullfile(fileparts(autoGatePath), 'run_epp_compiled');
    if exist(buildPath, 'dir')
        rmdir(buildPath, 's')
    end
    mkdir(buildPath);
    copyfile(fullfile(autoGatePath, 'CytoGenius.png'), buildPath);
    copyfile(fullfile(eppPath, '*.m'), buildPath);
    copyfile(fullfile(utilPath, '*.m'), buildPath);
    copyfile(fullfile(utilPath, '*.png'), buildPath);
    copyfile(fullfile(utilPath, '*.gif'), buildPath);
    copyfile(fullfile(utilPath, ['mexSptxModal.' mexext]), buildPath);
    copyfile(fullfile(utilPath, 'suh.jar'), buildPath);
    copyfile(fullfile(utilPath, 'transform.jar'), buildPath);
    if ismac
        copyfile(fullfile(eppPath, 'run_epp_mac.prj'), buildPath);
    elseif ispc
        copyfile(fullfile(eppPath, 'run_epp_pc.prj'), buildPath);
        copyfile(fullfile(utilPath, 'libfftw*.dll'), buildPath);
    else
        error('Only supported for Mac and PC');
    end
    cd(buildPath)
    if ismac
        open run_epp_mac.prj
    else
        open run_epp_pc.prj
    end

end
