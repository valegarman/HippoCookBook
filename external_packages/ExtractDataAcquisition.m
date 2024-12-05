function [outputArg] = ExtractDataAcquisition(filename)
%filename is the full name of the .doric file where we want to extract all
%the data from the DataAcquisition
%
%OutputArg is a structure with all the Data contained in the .doric file

if ~contains(filename,'.doric')
    filename = [filename '.doric'];
end

%DataAcquisition = h5info(filename,'/DataAcquisition');
DataAcquisition = h5info(filename);


%Recursive function to go find and extract all the data
    function [Dataset] = getalldata(H5)
        
        if isempty(H5.Groups)
            if ~isempty(H5.Datasets)
                
                Dataset_tmp = [];
                for k = 1:length(H5.Datasets)
                    Name = [H5.Name '/' H5.Datasets(k).Name];
                    Data = h5read(filename,Name);
                    
                    Dataset_tmp(k).Name = H5.Datasets(k).Name;
                    Dataset_tmp(k).Data = Data;
                    Dataset_tmp(k).DataInfo = H5.Datasets(k).Attributes;
                end
                
                Dataset.Name = strrep(H5.Name(2:end),'/','_');
                Dataset.Data = Dataset_tmp;
            end
        else
            Dataset = [];
            for k=1:length(H5.Groups)
                Dataset = [Dataset getalldata(H5.Groups(k))];
            end
        end
        
    end


   outputArg = getalldata(DataAcquisition);

end

