function colors = getColors(varargin)

if nargin < 1
    % Cell type colors
    cmap = brewermap(256,'Spectral');
    colors.camk2 = cmap(5,:);
    colors.camk2_dark =  cmap(5,:)/1.5;

    colors.sup = cmap(20,:);
    colors.deep = cmap(1,:);

    colors.id2 = cmap(100,:);
    colors.id2_dark = cmap(100,:)/1.5;

    colors.sncg = cmap(75,:);
    
    colors.sst = cmap(200,:);
    colors.sst_dark = cmap(200,:)/1.5;
    
    colors.pv = cmap(240,:);
    colors.pv_dark = cmap(240,:)/2;

    colors.axo = cmap(215,:);
    
    cmap = brewermap([],'PiYG');
    colors.vip = cmap(50,:);
    colors.vip_dark = cmap(50,:)/1.5;
    
    colors.vip_ww = cmap(50,:);
    colors.vip_ww_dar = cmap(50,:)/1.5;

    cmap = brewermap([],'PRGn');
    colors.vip_nw = cmap(50,:);
    colors.vip_nw_dark = cmap(50,:)/1.5;
    
    cmap = brewermap(100,'RdBu');
    colors.pyr = cmap(25,:);
    colors.pyr_dark = cmap(5,:);
    colors.pyr_light = cmap(40,:);
    
    colors.int = cmap(75,:);
    colors.int_dark = cmap(95,:);
    colors.int_light = cmap(60,:);

    colors.sncg = cmap(75,:);
else
    colors = [];
    input_colors = varargin{1};
    if ~iscell(input_colors)
        temp = input_colors;
        input_colors = cell(1);
        input_colors{1} = temp;
    end

    for ii = 1:length(input_colors)
        colors = [colors; selectColor(input_colors{ii})];
    end
end
end

function output = selectColor(input)
switch lower(input)
    case 'camk2_deep'
        cmap = brewermap(256,'Spectral');
        output = cmap(1,:);
    case 'camk2_sup'
        cmap = brewermap(256,'Spectral');
        output = cmap(20,:);
    case 'pv+'
        cmap = brewermap(256,'Spectral');
        output = cmap(240,:);
    case 'pv'
        cmap = brewermap(256,'Spectral');
        output = cmap(240,:);
    case 'pvalb'
        cmap = brewermap(256,'Spectral');
        output = cmap(240,:);
    case 'sst+'
        cmap = brewermap(256,'Spectral');
        output = cmap(200,:);
    case 'sst'
        cmap = brewermap(256,'Spectral');
        output = cmap(200,:);
    case 'id2+'
        cmap = brewermap(256,'Spectral');
        output = cmap(100,:);
    case 'id2'
        cmap = brewermap(256,'Spectral');
        output = cmap(100,:);
    case 'vip+'
        cmap = brewermap([],'PiYG');
        output = cmap(50,:);
    case 'vip'
        cmap = brewermap([],'PiYG');
        output = cmap(50,:);
    case 'camk2+'
        cmap = brewermap(256,'Spectral');
        output = cmap(5,:);
    case 'camk2'
        cmap = brewermap(256,'Spectral');
        output = cmap(5,:);
    case 'pyr'
        cmap = brewermap(256,'Spectral');
        output = cmap(5,:);
    case 'int'
        cmap = brewermap(256,'Spectral');
        output = cmap(220,:);
    case 'pyramidal cell'
        cmap = brewermap(256,'Spectral');
        output = cmap(5,:);
    case 'narrow interneuron'
        cmap = brewermap(256,'Spectral');
        output = cmap(220,:);
    case 'wide interneuron'
        cmap = brewermap(256,'Spectral');
        output = cmap(220,:);
    otherwise 
        output = [.5 .5 .5];
end

end