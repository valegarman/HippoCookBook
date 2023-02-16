
function stacked = stackSessionResult(toStack,numCells,varargin)
% stack all fields of a cell array of structures
% MV 2022

p = inputParser;
addRequired(p,'toStack',@iscell);
addRequired(p,'numCells',@isnumeric);
addParameter(p,'exampleSession',1,@isnumeric);
addParameter(p,'activateTry',true,@islogical);

parse(p, toStack,numCells,varargin{:});

exampleSession = p.Results.exampleSession;
activateTry = p.Results.activateTry;

if length(toStack)<1
    error('Input data is empty');
end

if length(numCells)<1
    error('Number of cells field not provided');
end

if size(numCells,1) < size(numCells,2)
    numCells = numCells';
end

for mm = 1:length(toStack)
    if isempty(toStack{mm})
        toStack{mm} = NaN;
    end
end 

% In case exampleSession(1) is NaN
for ii = 1:length(toStack)
    while ~isstruct(toStack{exampleSession})
        exampleSession = exampleSession + 1;
        if isstruct(toStack{exampleSession})
            continue;
        end
    end
end

names = fieldnames(toStack{exampleSession});
clonned = zeros(size(toStack));
% if nan, clon data of same size as numCells with nan values
for mm = 1:length(toStack)
    if ~isstruct(toStack{mm}) && isnan(toStack{mm})
        % find in example session, fields with size as in numCells
        for ii = 1:length(names)
            if isnumeric(toStack{exampleSession}.(names{ii})) %any(size(toStack{exampleSession}.(names{ii}))==numCells(exampleSession))
                if ischar(toStack{exampleSession}.(names{ii}))
                    keyboard;
                    toStack_clon{mm}.(names{ii}) = toStack{exampleSession}.(names{ii});
                else
                    example_data = nan*toStack{exampleSession}.(names{ii});
                    size_ex = size(example_data);
                    size_sess = size_ex; 
                    size_sess(find(size_ex == numCells(exampleSession))) = numCells(mm);
                    toStack_clon{mm}.(names{ii}) = nan(size_sess);
                end
                clonned(mm) = 1;
                
            else % if not numeric
                toStack_clon{mm}.(names{ii}) = nan(1);
            end
        end
    end
end

% flaten structure! if there is any anidated (1 level) structure numCells si
for ii = 1:length(names)
    if isstruct(toStack{exampleSession}.(names{ii}))
        l2_names = fieldnames(toStack{exampleSession}.(names{ii}));
        for jj = 1:length(l2_names)
            if isnumeric(toStack{exampleSession}.(names{ii}).(l2_names{jj})) &&  any(size(toStack{exampleSession}.(names{ii}).(l2_names{jj})) == numCells(exampleSession))
                fprintf('Flattening %s... \n', l2_names{jj}); %\n
                for mm = 1:length(toStack)
                    if isstruct(toStack{mm}) && isfield(toStack{mm}.(names{ii}),(l2_names{jj}))
                        toStack{mm}.([names{ii} '_' l2_names{jj}]) = toStack{mm}.(names{ii}).(l2_names{jj});
                    elseif isstruct(toStack{mm}) && ~isfield(toStack{mm}.(names{ii}),(l2_names{jj})) % if doesn't exist, create as nan matrix
                        example_data = nan*toStack{exampleSession}.([names{ii} '_' l2_names{jj}]);
                        size_ex = size(example_data);
                        size_sess = size_ex; 
                        size_sess(find(size_ex == numCells(exampleSession))) = numCells(mm);
                        toStack{mm}.([names{ii} '_' l2_names{jj}]) = nan(size_sess);
                    elseif isnan(toStack{mm}) % if nan, clone
                        if  ischar(toStack{exampleSession}.([names{ii} '_' l2_names{jj}]))
                            toStack_clon{mm}.([names{ii} '_' l2_names{jj}]) = ...
                                toStack{exampleSession}.([names{ii} '_' l2_names{jj}]);
                        else
                            example_data = nan*toStack{exampleSession}.([names{ii} '_' l2_names{jj}]);
                            size_ex = size(example_data);
                            size_sess = size_ex; 
                            size_sess(find(size_ex == numCells(exampleSession))) = numCells(mm);
                            toStack_clon{mm}.([names{ii} '_' l2_names{jj}]) = nan(size_sess);
                        end
                        clonned(mm) = 1;    
                    end
                end
            elseif isstruct(toStack{exampleSession}.(names{ii}).(l2_names{jj}))
                l3_names = fieldnames(toStack{exampleSession}.(names{ii}).(l2_names{jj}));
                
                for kk = 1:length(l3_names)
                    if isnumeric(toStack{exampleSession}.(names{ii}).(l2_names{jj}).(l3_names{kk})) & any(size(toStack{exampleSession}.(names{ii}).(l2_names{jj}).(l3_names{kk})) == numCells(exampleSession))
                        fprintf('Flattening %s... \n', l3_names{kk});
                        for mm = 1:length(toStack)
                           if isstruct(toStack{mm}) && isfield(toStack{mm}.(names{ii}).(l2_names{jj}),(l3_names{kk}))
                               toStack{mm}.([names{ii} '_' l2_names{jj} '_' l3_names{kk}]) = toStack{mm}.(names{ii}).(l2_names{jj}).(l3_names{kk});
                           elseif isstruct(toStack{mm}) && ~isfield(toStack{mm}.(names{ii}).(l2_names{jj}),(l3_names{kk})) % if it doesn't exist, create as nan matrix
                               example_data = nan*toStack{exampleSession}.([names{ii} '_' l2_names{jj} '_' l3_names{kk}]);
                                size_ex = size(example_data);
                                size_sess = size_ex; 
                                size_sess(find(size_ex == numCells(exampleSession))) = numCells(mm);
                                toStack{mm}.([names{ii} '_' l2_names{jj} '_' l3_names{kk}]) = nan(size_sess);
                           end
                        end
                    end
                end
                
            elseif strcmpi('timestamps',l2_names{jj})
                for mm = 1:length(toStack)
                    if isstruct(toStack{mm}) && isstruct(toStack{mm}.(names{ii}))
                        toStack{mm}.([names{ii} '_' l2_names{jj}]) = toStack{mm}.(names{ii}).(l2_names{jj});
                    % elseif isnan(toStack{mm})
                    else
                        toStack_clon{mm}.([names{ii} '_' l2_names{jj}]) = toStack{exampleSession}.(names{ii}).(l2_names{jj});
                    end
                end
            end
        end
    end
end

clonned = find(clonned);
for ii = 1:length(find(clonned))
    toStack(clonned(ii)) = toStack_clon(clonned(ii));
end

names = fieldnames(toStack{exampleSession});
for ii = 1:length(names)
    try
    
    for jj = 1:length(toStack)
       size_field(jj,:) = size(toStack{jj}.(names{ii})); 
       is_vector(jj) = isvector(toStack{jj}.(names{ii}));
       type_field{jj} = class(toStack{jj}.(names{ii}));
       ndims_field(jj) = ndims(toStack{jj}.(names{ii}));
    end
    [dims] = max(size_field);
    type_field = unique(type_field);
    
    for jj = 1:size(size_field,2)
        dim_cells_data(jj) = isequal(size_field(:,jj),numCells);
    end

    if length(type_field) > 1
        % keyboard;
        if any(strcmpi(type_field,'struct')) && ~isempty(clonned) % unless a structure was clonned
            type_field = 'struct';
        elseif any(strcmpi(type_field,'cell')) && ~isempty(clonned) % unless a structure was clonned
            type_field = 'struct';
        elseif any(strcmpi(type_field,'logical')) && ~isempty(clonned) % unless a structure was clonned
            type_field = 'double';
        else
            type_field = 'multiple';
            warning('More than one type of data per field!')
        end
    end
    if length(unique(ndims_field)) > 1
        error('More than one dimension setting of data per field!')
    end
    dim_sorted = 1:unique(ndims_field);
    dim_sorted = [dim_sorted(dim_cells_data) dim_sorted(~dim_cells_data)];
    dims = dims(dim_sorted);
    
    if ~isempty(find(dim_cells_data)) && strcmpi(type_field,'double')
        stacked.(names{ii}) = [];
        for jj = 1:length(toStack)
            % stack by the cells in first dimension
            values = permute(toStack{jj}.(names{ii}),dim_sorted);
            temp = nan(dims);
            temp(1:size(values,1),1:size(values,2),1:size(values,3),1:size(values,4),1:size(values,5),1:size(values,6)) = values;
            temp = temp(1:numCells(jj),:,:,:,:,:); %temp(isnan(temp(:,1)),:,:,:,:,:) = [];
            stacked.(names{ii}) = cat(1,stacked.(names{ii}),temp);
        end
        clear temp values
    elseif ~isempty(find(dim_cells_data)) && strcmpi(type_field,'cell') && any(dims==1) && length(dims)<3
        stacked.(names{ii}) = [];
        for jj = 1:length(toStack)
            % stack by the cells in first dimension
            values = permute(toStack{jj}.(names{ii}),dim_sorted);
            stacked.(names{ii}) = cat(1,stacked.(names{ii}),values);
        end
    elseif strcmpi((names{ii}),'timestamps') || contains((names{ii}),'timestamps')
        stacked.(names{ii}) = toStack{exampleSession}.(names{ii});
    elseif strcmpi(type_field,'double') && any(dims==1) && length(dims)<3
        % multiply for cell numbers to stack
        stacked.(names{ii}) = [];
        for jj = 1:length(toStack)
            values = toStack{jj}.(names{ii});
            if size(values,1) < size(values,2)
                values = values';
            end
            values = (values * ones(1, numCells(jj)))';
            temp = nan(numCells(jj),max(dims));
            temp(1:size(values,1),1:size(values,2)) = values;
            stacked.(names{ii}) = cat(1,stacked.(names{ii}),temp);
        end
    elseif strcmpi(type_field,'struct') || strcmpi(type_field,'cell') 
        fprintf('Field %s can not be stacked \n', names{ii}); %\n
    end
    clear size_field is_vector is_vector ndims_field dims dim_sorted dim_cells_data type_field
    catch
        warning('One field failed to be stacked!!');
    end
end
stacked.sessionNumber = [];
for jj = 1:length(numCells)
    stacked.sessionNumber = cat(1,stacked.sessionNumber,jj*ones(numCells(jj),1));
end

end