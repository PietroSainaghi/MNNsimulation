function [xInit, randomizedIndices] = setInitialConditions(LinkPropertiesStruct,LatticeGeometryStruct,OptimizerDataStruct,icIter)

%% function inputs

% link properies structure - LinkPropertiesStruct
kLinMax = LinkPropertiesStruct.kLinMax;
kLinMin = LinkPropertiesStruct.kLinMin;

% optimization problem structure - OptimizerDataStruct
discORcont = OptimizerDataStruct.discORcont;
icType = OptimizerDataStruct.icType;
PossibleStiffnessArray = OptimizerDataStruct.PossibleStiffnessArray;
deterministicRNG = OptimizerDataStruct.deterministicRNG;

% lattice structure - LatticeGeometryStruct
Nbeams = LatticeGeometryStruct.Nbeams;

%% computations

switch discORcont
    
    % create ICs for discrete dataset
    case 'discrete'
        if icType == 1
            % all initial stiffness states at max value
            xInit = ones(Nbeams,1)*max(PossibleStiffnessArray);
            randomizedIndices = ones(Nbeams,1)*length(PossibleStiffnessArray);
        elseif icType == 2
            % all initial stiffness states at min value
            xInit = ones(Nbeams,1)*min(PossibleStiffnessArray);
            randomizedIndices = ones(Nbeams,1);
            
        elseif icType == 3
            % randomized initial state
            
            % enable deterministic RNG
            if deterministicRNG == 1
                rng(icIter*3);
            end
            randomizedIndices = randi(length(PossibleStiffnessArray),Nbeams,1);
            xInit = zeros(Nbeams,1);
            for eachx = 1:Nbeams
                xInit(eachx) = PossibleStiffnessArray(randomizedIndices(eachx));
            end
        end
    case 'continuous'
        
        % set array to zero since it is not used
        randomizedIndices = zeros(Nbeams,1);
        
        if icType == 1
            % all initial stiffness states at max value
            xInit = ones(Nbeams,1)*kLinMax;
            
        elseif icType == 2
            % all initial stiffness states at min value
            xInit = ones(Nbeams,1)*kLinMin;
            
        elseif icType == 3
            % randomized initial state
            
            % enable deterministic RNG
            if deterministicRNG == 1
                rng(icIter*3);
            end
            
            xInit = kLinMin + (kLinMax - kLinMin) * rand(Nbeams,1);
            
        end
end