
function stacked = stackSessionResult(toStack)
% stack all fields of a cell array of structures
% MV 2022

if length(toStack)<1
    error('Input data is empty');
end

names = fieldnames(toStack{1});
for ii = 1:length(names)
    for jj = 1:length(toStack)
       [dim1(jj), dim2(jj), dim3(jj), dim4(jj) dim5(jj)] = size(toStack{jj}.(names{ii})); 
       is_vector(jj) = isvector(toStack{jj}.(names{ii}));
       type_field{jj} = class(toStack{jj}.(names{ii}));
       ndims_field(jj) = ndims(toStack{jj}.(names{ii}));
    end
    dims = [max(dim1) max(dim2) max(dim3) max(dim4) max(dim5)];
    
    type_field = unique(type_field);
    if length(type_field) > 1
        error('More than one type of data per field!')
    end
    if length(unique(ndims_field)) > 1
        error('More than one dimension setting of data per field!')
    end
    dims = dims(1:unique(ndims_field)); idx_dims = idx_dims(1:unique(ndims_field));
    
    if strcmpi((names{ii}),'timestamps')
        stacked.(names{ii}) = toStack{jj}.(names{ii});
    elseif strcmpi((names{ii}),'channels')
        stacked.(names{ii}) = toStack{jj}.(names{ii});
    elseif (strcmpi(type_field,'logical') || strcmpi(type_field,'double')...
            || strcmpi(type_field,'int32')) && any(~is_vector)
        stacked.(names{ii}) = [];
        for jj = 1:length(toStack)
            % stack by the cells in first dimension
            temp = nan(dims);
            values = toStack{jj}.(names{ii});
            temp(1:size(values,1),1:size(values,2),1:size(values,3),1:size(values,4),1:size(values,5)) = values;
            temp(isnan(temp(:,1)),:,:,:,:) = [];
            stacked.(names{ii}) = cat(1,stacked.(names{ii}),temp);
        end
    end
    
end

end