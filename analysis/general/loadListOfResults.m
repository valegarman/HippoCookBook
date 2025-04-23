function list = loadListOfResults(project)

list = [];
switch lower(project)
    case 'bibliocampus'
        list = {'optogeneticResponse','averageCCG','ripples_psth','slowOscillations_psth','theta_*.PhaseLockingData','thetaREM*.PhaseLockingData',...
            'thetaRun*.PhaseLockingData','lgamma*.PhaseLockingData','hgamma*.PhaseLockingData','ripple*.PhaseLockingData','spatialModulation','placeFields','behavior.cellinfo','ACGPeak',...
            'speedCorr','uLEDResponse.cellinfo'};


    case 'lightininh'
        list = {'optogeneticResponse','averageCCG','ripples_psth','slowOscillations_psth','theta_*.PhaseLockingData','thetaREM*.PhaseLockingData',...
        'thetaRun*.PhaseLockingData','lgamma*.PhaseLockingData','hgamma*.PhaseLockingData','ripple*.PhaseLockingData','spatialModulation','placeFields','behavior.cellinfo','ACGPeak',...
        'speedCorr.cellinfo','uLEDResponse.cellinfo','lightSpikeCollisions','spikeCCGchange','lightSpikeCollisions_control'};

    otherwise
        list = {'ripples_psth','slowOscillations_psth','theta_*.PhaseLockingData','lgamma*.PhaseLockingData','hgamma*.PhaseLockingData','ripple*.PhaseLockingData','spatialModulation','placeFields','behavior.cellinfo','ACGPeak'};
end
end