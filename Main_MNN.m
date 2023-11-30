clear
clc
close all

% Author: Pietro Sainaghi
% code modified from: Binary_MinStiffness by Ryan Lee
% works with the other functions in this folder, written by Reinier Kuppens
% and Ryan Lee

% this code studies a given MNN lattice ability to learn behaviors with
% respect to various lattice parameters

%% Add Function Folders

addpath([pwd,'\Main_Functions'])

%% Individual Link Parameters

% these should be set based on CADs and FEA for the links tobe used in the
% MNN

% length of each beam
L1                  = 0.156; % m
%X-component of beam under 60 degrees / 1/3 pi rad
L2                  = L1*cosd(60); 
%Y-component of beam under 60 degrees / 1/3 pi rad
L3                  = L1*cosd(30);     

% maximum link elongation
MaxLinkElongation   = 0.0025; % m

% maximum link compression
MaxLinkCompression   = 0.0025; % m

% stiffness values
    % Force acting perpendicular while rotation is fixed: utcross=0;
k55                  = 237.1; %N/m
    % Force acting perpendicular with rotation of tip allowed
k56                  = 2.11; %N/rad
    % Moment acting on end while perpendicular translation of tip is allowed 
k66                  = 25.3; %N*m/rad

% passive axial stiffness
kLinPassive             = 2120; %N/m

% maximum controllable stiffness
kLinMax                 = 2300; %N/m

% minimum controllable stiffness
kLinMin                 = -2000; %N/m

% assemble LinkPropertiesStruct
LinkPropertiesStruct.L1 = L1;
LinkPropertiesStruct.L2 = L2;
LinkPropertiesStruct.L3 = L3;
LinkPropertiesStruct.MaxLinkElongation = MaxLinkElongation;
LinkPropertiesStruct.MaxLinkCompression = MaxLinkCompression;
LinkPropertiesStruct.k55 = k55;
LinkPropertiesStruct.k56 = k56;
LinkPropertiesStruct.k66 = k66;
LinkPropertiesStruct.kLinPassive = kLinPassive;
LinkPropertiesStruct.kLinMax = kLinMax;
LinkPropertiesStruct.kLinMin = kLinMin;


%% Lattice Gemetry Parameters

% select type of lattice to simulate 
latticeType          = 1; %1) triangular 2) square

% number of nodes that receive forces and measure displacements
    % set as array to test multiple ones at once
NinputANDoutputArray      = [8];

% array of numbers of layers
    % set as array to test multiple ones at once
NlayersArray              = [8];

% dimensions in space
DOI = 3;

% assemble LatticeGeometryStruct
LatticeGeometryStruct.latticeType = latticeType;
LatticeGeometryStruct.NinputANDoutputArray = NinputANDoutputArray;
LatticeGeometryStruct.NlayersArray = NlayersArray;
LatticeGeometryStruct.DOI = DOI;

%% Behavior Configuration

% type of study 
caseType        = 3;  % 1) sinusoid behavior     2) random forces       3) saved behavior
% GeneratePresetBehavior will create .mat files compatible with caseType 3

% number of force behaviors in randomized behaviors 
    % set as array to test multiple ones at once
    % can only be 2 for caseType = 1
    % will be overwritten for caseType = 3
NcasesArray             = [2]; 

% minimum force allowed as fraction of maximum force
threshold       = 0.3;

% maximum allowed elongation of random behavior
    % set as array to test multiple ones at once
    % will be overwritten for caseType = 3
Elongation_maxArray  = [0.0025]; % units of m

% aplitude of sine wave for sinusoid behavior
    % set same value as Elongation_maxArray for looping convenience, since
    % random and sinusoid are not tested in the same run
    % will be overwritten for caseType = 3
dxArray              = Elongation_maxArray; % units of m

% maximum allowed input force
    % set as array to test multiple ones at once
    % will be overwritten for caseType = 3
MaxForceArray = [8]; % units of N

% type of RNG for randomized behavior
RNGtype = 'deterministic'; % options: 'dateandtime' 'deterministic' 'nocontrol'


% assemble BehaviorStruct
BehaviorStruct.caseType = caseType;
BehaviorStruct.NcasesArray = NcasesArray;
BehaviorStruct.threshold = threshold;
BehaviorStruct.Elongation_maxArray = Elongation_maxArray;
BehaviorStruct.dxArray = dxArray;
BehaviorStruct.MaxForceArray = MaxForceArray;
BehaviorStruct.RNGtype = RNGtype;


%% Optimizer Configuration

% force scaling
ForceScaling        = false;

% number of runs for each set of behaviors
startPts            = 1;

% select treshlod for desired precision
errorChangeThreshold  = 1e-15;

% type of initial condition
icType          = 3; %1) all max, 2) all min, 3) rand

% wether to use MSE or scaled normalized MSE
doScaleMSE = 0; % 0) use MSE, 1) use unitless MSE

% select units for MSE
    % 1 = m, 100 = cm, 1000 = mm, 1000000 = um
MSEunits = 1000;

% whether to run a discrete or continuous optimization
    % discrete only works with GA
discORcont = 'discrete'; % 'discrete' or 'continuous'

% set diiscrete set of stiffness values
    % unused if discORcont = 'continuous'
PossibleStiffnessArray = linspace(-2000, 2300, 4300); %N/m

% type of optimizer to use
    % set as cell to test multiple ones at once
optimizerArray = {'GA'}; % 'GA', 'SQP'

% whether to use deterministic random initial conditions
deterministicRNG = 0;

% whether to set computational cache in separate folder
% will be same as results if false
% set to true if using GA or other algorithm with lots of function
% evaluations
SpecifyCache = 1;


% assemble OptimizerDataStruct
OptimizerDataStruct.ForceScaling = ForceScaling;
OptimizerDataStruct.startPts = startPts;
OptimizerDataStruct.errorChangeThreshold = errorChangeThreshold;
OptimizerDataStruct.icType = icType;
OptimizerDataStruct.doScaleMSE = doScaleMSE;
OptimizerDataStruct.discORcont = discORcont;
OptimizerDataStruct.optimizerArray = optimizerArray;
OptimizerDataStruct.deterministicRNG = deterministicRNG;
OptimizerDataStruct.PossibleStiffnessArray = PossibleStiffnessArray;
OptimizerDataStruct.MSEunits = MSEunits;
OptimizerDataStruct.SpecifyCache = SpecifyCache;


%% Plotting Options Configuration

% scale width of deformed parts
    % Not yet implemented
scalingterm = 1;

% plot undeformed lattice
plotUndeformed = 1;

% plot stiffness combinations in lattice
plotColoredLattice = 0;

% plot deformed lattice
    % set to 1 only if the loop inludes one lattice and one set of behaviors
plotDeformed = 0;

% plot endpoints
    % set to 1 only if the loop inludes lattice and one set of behaviors
plotEndpoints = 1;

% assemble plotOptionsStruct
plotOptionsStruct.scalingterm = scalingterm;
plotOptionsStruct.plotUndeformed = plotUndeformed;
plotOptionsStruct.plotDeformed = plotDeformed;
plotOptionsStruct.plotEndpoints = plotEndpoints;
plotOptionsStruct.plotColoredLattice = plotColoredLattice;


%% Save File

% flag to save learning output
IWantToSaveOutput = true;

if IWantToSaveOutput
    fileName = 'PlanarMeso_learningAbilityStudy';
    fileType = '.mat';
    t= datetime;
    t.Format = 'yyyy-MM-dd_@_HHmm-ss';
    t = char(t);
    fullName = [fileName,'_','run@',t, fileType] ;
    
    savePath = uigetdir(pwd,'Select output save location');
    folderName = 'MNN_study';
    fullPath = [savePath,'\',folderName,'\', fullName];
    
    mkdir([savePath,'\',folderName]);
end

%% select pregen behaviors

% only happens if caseType = 3
if caseType == 3
    
    % TODO
    % To improve
    [fileList, dataPath] = uigetfile('*.mat','Select preset behavior','MultiSelect', 'on');
    NumBehFiles = 1;
    load([dataPath,fileList]);
    BehaviorStruct.fileList = fileList;
    BehaviorStruct.dataPath = dataPath;
    NcasesArray = [PregenBehStruct.Ncases];
    BehaviorStruct.NcasesArray = [PregenBehStruct.Ncases];
    
    % check if preset behavior and lattice geometry setting set in code
    % here are compatible
    if (length(LatticeGeometryStruct.NinputANDoutputArray) > 1) || (length(LatticeGeometryStruct.NinputANDoutputArray) > 1)
        error('You must only specify one lattice size if you are using preset behaviors')
    end
    if LatticeGeometryStruct.NinputANDoutputArray ~= PregenBehStruct.NinputANDoutput
        error('Preset behavior structure is incompatible with your current lattice geometry, wrong NinputANDoutput')
    end
    if LatticeGeometryStruct.NlayersArray ~= PregenBehStruct.Nlayers
        error('Preset behavior structure is incompatible with your current lattice geometry, wrong Nlayers')
    end
else
    NumBehFiles = 1;
end

%% Select Computation Cache

if OptimizerDataStruct.SpecifyCache == 0
    
    fileNameCache = 'PlanarMeso_learningAbilityStudy';
    fileType = '.mat';
    t= datetime;
    t.Format = 'yyyy-MM-dd_@_HHmm-ss';
    t = char(t);
    cacheName = [fileNameCache,'_','run@',t, '_chachedata',fileType] ;
    
    savePathCache = savePath;
    folderNameCache = 'MNN_study_cache';
    
    CachePath = [savePathCache,'\',folderNameCache,'\', cacheName];
    
    mkdir([savePathCache,'\',folderNameCache]);
    
elseif OptimizerDataStruct.SpecifyCache == 1
    
    fileNameCache = 'PlanarMeso_learningAbilityStudy';
    fileType = '.mat';
    t= datetime;
    t.Format = 'yyyy-MM-dd_@_HHmm-ss';
    t = char(t);
    cacheName = [fileNameCache,'_','run@',t, '_chachedata',fileType] ;
    
    savePathCache = uigetdir(pwd,'Select cache location');
    folderNameCache = 'MNN_study_cache';
    
    CachePath = [savePathCache,'\',folderNameCache,'\', cacheName];
    
    mkdir([savePathCache,'\',folderNameCache]);
end



save('CacheDirectory.mat','CachePath');


%% Main Loops

% overall run counter
NoptRuns = 0;

% loop through number of inputs and output nodes
for NinputANDoutputIter = 1:length(NinputANDoutputArray)
    LatticeGeometryStruct.NinputANDoutput = LatticeGeometryStruct.NinputANDoutputArray(NinputANDoutputIter);
    
    % loop through number of layers
    for NlayersIter = 1:length(NlayersArray)
        LatticeGeometryStruct.Nlayers = LatticeGeometryStruct.NlayersArray(NlayersIter);
        
        % loop through output amplitudes
        for dxIter = 1:length(dxArray)
            BehaviorStruct.dx = BehaviorStruct.dxArray(dxIter);
            BehaviorStruct.Elongation_max = BehaviorStruct.Elongation_maxArray(dxIter);
            
            % loop through input force magnitudes
            for MaxForceIter = 1:length(MaxForceArray)
                BehaviorStruct.MaxForce = BehaviorStruct.MaxForceArray(MaxForceIter);
                
                % loop through number of cases
                for NcasesIter = 1:length(NcasesArray)
                    BehaviorStruct.Ncases = BehaviorStruct.NcasesArray(NcasesIter);
                    
                    % loop though optimizers
                    for optimizerIter = 1:length(OptimizerDataStruct.optimizerArray)
                        OptimizerDataStruct.optimizer = OptimizerDataStruct.optimizerArray{optimizerIter};
                        
                        % Loop through number of behaviors
                        for behaviorSetIter = 1:NumBehFiles
                            
                            % load behavior if from file
                            if caseType == 3
                                load([dataPath,fileList])
                                BehaviorStruct.Target = PregenBehStruct.Target;
                                BehaviorStruct.forces = PregenBehStruct.forces;
                                BehaviorStruct.Ncases = PregenBehStruct.Ncases;
                                BehaviorStruct.Elongation_max = PregenBehStruct.Elongation_max;
                                BehaviorStruct.MaxForce = PregenBehStruct.MaxForce;
                                BehaviorStruct.dx = PregenBehStruct.dx;
                            end
                            
                            % set overall run number
                            NoptRuns = NoptRuns + 1;
                            % debug study display
                            disp(['Design study run number: ', num2str(NoptRuns),' / ',num2str(length(optimizerArray)*length(NcasesArray)*length(MaxForceArray)*length(dxArray)*length(NlayersArray)*length(NinputANDoutputArray))])
                            disp(['Number of input/output nodes: ', num2str(LatticeGeometryStruct.NinputANDoutput)])
                            disp(['Number of layers: ', num2str(LatticeGeometryStruct.Nlayers)])
                            if BehaviorStruct.caseType == 1
                                disp(['Behavior type: sinusoid'])
                            elseif BehaviorStruct.caseType == 2
                                disp(['Behavior type: random',', ',BehaviorStruct.RNGtype, ' RNG'])
                            elseif BehaviorStruct.caseType == 3
                                disp(['Behavior type: saved'])
                            end
                            disp(['Number of sets of behaviors: ', num2str(NumBehFiles)])
                            disp(['Behavior Output Amplitude: ', num2str(BehaviorStruct.Elongation_max)])
                            disp(['Behavior max force: ', num2str(BehaviorStruct.MaxForce)])
                            disp(['Number of behaviors to optimize: ', num2str(BehaviorStruct.Ncases)])
                            disp(['Type of optimizer: ', (OptimizerDataStruct.optimizer)])
                            disp(['Optimization space: ', (OptimizerDataStruct.discORcont)])
                            
                            % set deterministic RNG for random behavior
                            switch BehaviorStruct.RNGtype
                                case 'deterministic'
                                    BehaviorStruct.RNG_seed_det = NoptRuns*3;
                            end
                            
                            
                            % generate MNN lattice structure
                            LatticeGeometryStruct = packageLattice(LinkPropertiesStruct, LatticeGeometryStruct);
                            
                            
                            if (caseType == 1) || (caseType == 2)
                                % generate behavior
                                BehaviorStruct = packageBehavior(LatticeGeometryStruct,BehaviorStruct);
                            end



                            % compute MSE scale factor
                            % unused if OptimizerDataStruct.doScaleMSE = 0
                            if OptimizerDataStruct.doScaleMSE == 1
                                MSEscale = computeMSEScaleFactor(LatticeGeometryStruct,BehaviorStruct);
                                % store in structure OptimizerDataStruct
                                OptimizerDataStruct.MSEscale = MSEscale;
                            else
                                % store in structure OptimizerDataStruct
                                % zero means unused
                                OptimizerDataStruct.MSEscale = 0;
                            end

                            
                            
                            % set up finite truss structure (finite element used equivalently)
                            [DOFFinal,Final, F, DOFnodes,k_base,R6,RT6,DOF,Degrees_per_element] = FEM_setup(LinkPropertiesStruct,LatticeGeometryStruct,BehaviorStruct);
                            
                            %store in structure FEMStruct
                            FEMStruct.DOFFinal = DOFFinal;
                            FEMStruct.Final = Final;
                            FEMStruct.F = F;
                            FEMStruct.DOFnodes = DOFnodes;
                            FEMStruct.k_base = k_base;
                            FEMStruct.R6 = R6;
                            FEMStruct.RT6 = RT6;
                            FEMStruct.DOF = DOF;
                            FEMStruct.Degrees_per_element = Degrees_per_element;
                            
                            % run MNN with given specifications
                            ResultsStruct = MNN_run(...
                                LinkPropertiesStruct,...
                                LatticeGeometryStruct,...
                                BehaviorStruct,...
                                OptimizerDataStruct,...
                                FEMStruct,...
                                plotOptionsStruct);
                            
                            % save to individual cell
                            SweepOutput{...
                                NinputANDoutputIter,...
                                NlayersIter,...
                                dxIter,...
                                MaxForceIter,...
                                NcasesIter,...
                                optimizerIter} = ResultsStruct;
                            
                            
                            disp(newline)
                            disp(newline)
                            
                        end % end behavior set loop
                    end % end optimizer loop
                end % end number of behaviors loop
            end % end input force loop
        end % end output amplitude loop
    end % end number of layers loop
end % end number of input nodes loop

%% Save Ouput to File

if IWantToSaveOutput
    
    save(fullPath,'LinkPropertiesStruct','LatticeGeometryStruct','BehaviorStruct','OptimizerDataStruct','FEMStruct','plotOptionsStruct','SweepOutput')
end














