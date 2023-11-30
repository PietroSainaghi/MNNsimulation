function ResultsStruct = MNN_run(LinkPropertiesStruct,LatticeGeometryStruct,BehaviorStruct,OptimizerDataStruct,FEMStruct,plotOptionsStruct)


% Author: Pietro Sainaghi
% code modified from: Binary_MinStiffness by Ryan Lee
% works with the other functions in this folder, written by Reinier Kuppens
% and Ryan Lee

% this function run the optimization study for a single MNN, given a single
% set of behaviors, for multiple initial conditions

%% function structures

LinkPropertiesStruct; % stores all info about the individual links used in the MNN
LatticeGeometryStruct; % stores all info about the lattice used
BehaviorStruct; % stores all info about the force-displacemnt behavior
OptimizerDataStruct; % stores all info about the optimization study
FEMStruct; % storees all info about the finite truss method (FEM used interchangeably)
plotOptionsStruct; % plot options


%% function inputs

% optimization problem structure - OptimizerDataStruct
startPts = OptimizerDataStruct.startPts; % [1] number of initial conditions for a set of behaviors
optimizer = OptimizerDataStruct.optimizer;
discORcont = OptimizerDataStruct.discORcont;
PossibleStiffnessArray = OptimizerDataStruct.PossibleStiffnessArray;
errorThreshold = OptimizerDataStruct.errorChangeThreshold;


%% MNN optimization

% loop through initial conditions
for icIter = 1:startPts
    disp(['Initial Conditions Iteration Number: ',num2str(icIter),' / ',num2str(startPts)])
    % set randomized initial conditions
    [xInit, randomizedIndices] = setInitialConditions(LinkPropertiesStruct,LatticeGeometryStruct,OptimizerDataStruct,icIter);
    
    % run optimization study
    switch optimizer
        
        
        % genetic algorithm
        case 'GA'
            
            switch discORcont
                
                % discrete optimization
                case 'discrete'
                    % optimizer options and parameters
                    minIDX = 1;
                    maxIDX = length(PossibleStiffnessArray);
                    nvars = length(xInit);
                    population = 1000;
                    generations = 2000000;
                    stallgenerations = 50000;
                    % 'FunctionTolerance',errorThreshold,...
                    options = optimoptions('ga',...
                        'PlotFcn', @gaplotbestf,...
                        'PopulationSize',population, ...
                        'InitialPopulationMatrix', randomizedIndices',...
                        'MaxGenerations', generations,...
                        'ConstraintTolerance',errorThreshold,...
                        'FunctionTolerance',errorThreshold,...
                        'MaxStallGenerations', stallgenerations,...
                        'UseParallel', true,...
                        'OutputFcn',@gaoutfun); % saves output and population progression to .mat file(JANKY SOLUTION CAN BE IMPROVED)
                    
                    % function to be optimized
                    fun = @(x) ERROR_discrete(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x,PossibleStiffnessArray);
                    
                    % constraint equation
                    constraint = @(x) ELONGATION_discrete(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct, x,PossibleStiffnessArray);
                    
                    % optimization using GA
                    disp('Starting Optimization')
                    % [x_idx,fval,exitflag,output,~,scores] = ga(fun,nvars,[],[],[],[],ones(1,nvars)*minIDX,ones(1,nvars)*maxIDX,constraint,[1:1:nvars], options);
                    [x_idx,fval,exitflag,output,~,scores] = ga(fun,nvars,[],[],[],[],ones(1,nvars)*minIDX,ones(1,nvars)*maxIDX,[],[1:1:nvars],options);
                    
                    % optimal stiffness values
                    for eachx = 1:length(x_idx)
                        stiffnesscombofinal(eachx) = PossibleStiffnessArray(x_idx(eachx));
                    end
                    
                    % validate
                    finalerror = ERROR_discrete(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x_idx,PossibleStiffnessArray);
                    disp(['Final Error: ',num2str(finalerror)])
                    [nonhomogeneous, homogeneous] = ELONGATION_discrete(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct, x_idx,PossibleStiffnessArray);
                    coorddeformed=FEM(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,stiffnesscombofinal);
                    
                    % extract iteration progression data from mat file created by gaoutfun
                    load('CacheDirectory.mat')
                    load(CachePath)
                    
                    % output data
                    ResultsStruct(icIter).x_idx = x_idx;
                    ResultsStruct(icIter).x = stiffnesscombofinal;
                    ResultsStruct(icIter).elongationDiff = nonhomogeneous;
                    ResultsStruct(icIter).finalerror = finalerror;
                    ResultsStruct(icIter).itercount = output.generations;
                    ResultsStruct(icIter).funevaluations = output.funccount;
                    % ResultsStruct(icIter).GA_POPhistory = GA_POPhistory;
                    ResultsStruct(icIter).GA_scorehistory = GA_scorehistory;
                    ResultsStruct(icIter).GA_ERRhistory = GA_ERRhistory;
                    ResultsStruct(icIter).coord_deformed = coorddeformed;
                    
                    % discrete optimization
                case 'continuous'
                    
                    % upper and lower stiffness bounds
                    kLinMax = LinkPropertiesStruct.kLinMax;
                    kLinMin = LinkPropertiesStruct.kLinMin;
                    
                    % optimizer options and parameters
                    nvars = length(xInit);
                    population = 100;
                    generations = 1000;
                    stallgenerations = 50;
                    % 'FunctionTolerance',errorThreshold,...
                    options = optimoptions('ga',...
                        'PlotFcn', @gaplotbestf,...
                        'PopulationSize',population, ...
                        'InitialPopulationMatrix', randomizedIndices',...
                        'MaxGenerations', generations,...
                        'ConstraintTolerance',errorThreshold,...
                        'FunctionTolerance',errorThreshold,...
                        'MaxStallGenerations', stallgenerations,...
                        'UseParallel', true,...
                        'OutputFcn',@gaoutfun); % saves output and population progression to .mat file(JANKY SOLUTION CAN BE IMPROVED)
                    
                    % constraint equation
                    constraint = @(x) ELONGATION_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct, x);
                    
                    % function to be optimized
                    fun = @(x) ERROR_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
                    
                    % optimization using GA
                    disp('Starting Optimization')
                    [x,fval,exitflag,output,population,scores]  = ga(fun,nvars,[],[],[],[],kLinMin,kLinMax,constraint,options);
                    
                    
                    % validate
                    finalerror = ERROR_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
                    disp(['Final Error: ',num2str(finalerror)])
                    [nonhomogeneous, homogeneous] = ELONGATION_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct, x);
                    coorddeformed=FEM(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
                    
                    % extract iteration progression data from mat file created by gaoutfun
                    load('CacheDirectory.mat')
                    load(CachePath)
                    
                    % output data
                    ResultsStruct(icIter).x = x;
                    ResultsStruct(icIter).elongationDiff = nonhomogeneous;
                    ResultsStruct(icIter).finalerror = finalerror;
                    ResultsStruct(icIter).itercount = output.generations;
                    ResultsStruct(icIter).funevaluations = output.funccount;
                    % ResultsStruct(icIter).GA_POPhistory = GA_POPhistory;
                    ResultsStruct(icIter).GA_scorehistory = GA_scorehistory;
                    ResultsStruct(icIter).GA_ERRhistory = GA_ERRhistory;
                    ResultsStruct(icIter).coord_deformed = coorddeformed;
                    
            end % end discrete continuous switch
            
            
        case 'SQP'
            
            switch discORcont
                case 'continuous'
                    
                    % upper and lower stiffness bounds
                    kLinMax = LinkPropertiesStruct.kLinMax;
                    kLinMin = LinkPropertiesStruct.kLinMin;
                    
                    % optimization parameters
                    options = optimoptions('fmincon',...
                        'Algorithm','sqp',...
                        'Display','iter',...
                        'PlotFcn','optimplotfval',...
                        'MaxIterations',400,...,
                        'MaxFunctionEvaluations',1e6,...,
                        'TolConSQP',1e-5);
                    
                    % constraint equation
                    constraint = @(x) ELONGATION_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct, x);
                    
                    % function to be optimized
                    fun = @(x) ERROR_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
                    
                    % optimization using SQP
                    disp('Starting Optimization')
                    [x,fval,exitflag,output] = fmincon(fun,xInit,[],[],[],[],kLinMin,kLinMax,constraint,options);
                    
                    % validate
                    finalerror = ERROR_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
                    disp(['Final Error: ',num2str(finalerror)])
                    [nonhomogeneous, homogeneous] = ELONGATION_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct, x);
                    coorddeformed=FEM(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
                    
                    
                    % output data
                    ResultsStruct(icIter).x = x;
                    ResultsStruct(icIter).elongationDiff = nonhomogeneous;
                    ResultsStruct(icIter).finalerror = finalerror;
                    ResultsStruct(icIter).coord_deformed = coorddeformed;
            end
            
    end % end algorithm switch
    
    % plot endpoints
    if plotOptionsStruct.plotEndpoints == 1
        plotEndpointsFunction(LatticeGeometryStruct,BehaviorStruct,ResultsStruct(icIter));
    end
    % plot colored lattice
    if plotOptionsStruct.plotColoredLattice == 1
        plotColoredStiffnessValues(LinkPropertiesStruct, LatticeGeometryStruct, ResultsStruct(icIter).x, 'Main_Functions\colormap.mat');
    end
    
end % end initial conditions loop
