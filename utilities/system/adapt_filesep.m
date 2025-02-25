function [file_path_out] = adapt_filesep(file_path_in)
% adapt file path according to SO
% Manu Valero 2022
    switch computer
        case 'MACI64'
            file_path_in(strfind(file_path_in,'\'))='/';
        case 'MACA64'
            file_path_in(strfind(file_path_in,'\'))='/';
        case 'PCWIN64'
            file_path_in(strfind(file_path_in,'/'))='\';
        case 'GLNXA64'
            file_path_in(strfind(file_path_in,'\'))='/';
    end
    file_path_out = file_path_in;
end