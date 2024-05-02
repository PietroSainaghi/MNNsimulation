function [x] = analyticalGradientMethod(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct, xInit)

% Algorithm by Jiaji Chen to work with mass-spring model, rewritten for finite truss model by Pietro Sainaghi

% this function uses an analytical gradient descent method to optimize a MNN

%% Explanation of variables


%% input unpacking

% optimizer structure - OptimizerDataStruct
% iterMAX = OptimizerDataStruct.AGDmaxIterations;
iterMAX = 1;
learningRate = 0.01; % Learning Rate
Decay_interval = 50; % Learning rate decay interval
Decay_factor = 5; % Learning rate decay factor
% Batch_num = 1; % Minibatch size UNUSED but shows up in original algorithm
MSEscale = 20; % scaling term for MSE used in algorthm

% behavior structure - BehaviorStruct
Ncases = BehaviorStruct.Ncases; % [1] number of force behaviors to optimize towards
Target = BehaviorStruct.Target; % [NinputANDoutput, 2, Ncases] desired x and y coordinates for each output node
Elongation_max = BehaviorStruct.Elongation_max; % [1] maximum node displacement in output nodes

% lattice structure - LatticeGeometryStruct
connectivity = LatticeGeometryStruct.connectivity; % [Nbeams,2] connectivity of nodes for each link
Outputline = LatticeGeometryStruct.outputNodes; % [NinputANDoutput] indices of output nodes
NinputANDoutput = LatticeGeometryStruct.NinputANDoutput; % [1] number of input and output nodes in a given lattice
coord_initial = LatticeGeometryStruct.coord_initial; % [Nnodes,3] x,y,z coordinates of each node 
Nnodes = size(coord_initial,1); % [1] number of nodes in the lattice
DOI = LatticeGeometryStruct.DOI; % [1] dimensions in space = 3

% finite truss structure - FEMStruct
F_all = FEMStruct.F; % [NDOF, Ncases] location and value of forces within degrees of freedom format for each behavior
DOFnodes = FEMStruct.DOFnodes; % [] indices of nodes that are not grounded

%% generate K to x inversion function

[K_contribution] = invertStiffnessMatrix2AxialStiffness(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct);


%% initialize time history storage 

% list of loss function value at each iteration
Loss_hist = []; 

% history of stiffness combinations
Stiffness_hist = [];


%% Optmizer Loop

% assign initial condition on stiffness
x = xInit; % column vector

for iter = 1:iterMAX

    % update learning rate based on decay factor
    if mod(iter, Decay_interval == 0)
        learningRate = learningRate/Decay_factor;
    end

    % store stiffness combination for each iteration
    Stiffness_hist = [Stiffness_hist x];

    % compute loss function (scaled MSE)
    % works with both unitless MSE and unit MSE, but may perform differently based on MSEscale
    MSE = ERROR_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
    L = MSEscale^2/2 * MSE;

    % store loss in vector to plot
    Loss_hist = [Loss_hist L];

    % extract FEM properties used in algorithm
    [coorddeformed,Kdof,Fdof,udof]=FEM(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);


    % store all displacements, including grounds
    % there are easier faster ways of doing this, but this is good for debugging as it shows all nodes
    u_all = zeros(Nnodes,DOI,Ncases);
    for behIDX=1:Ncases
        nodeIDX = 1;
        for freeIDX=DOFnodes %loop through the nodes that are not fixed
            u_all(freeIDX,:,behIDX)=transpose(udof(DOI*(nodeIDX-1)+1:DOI*(nodeIDX-1)+DOI,behIDX));
            nodeIDX = nodeIDX + 1;
        end
    end


   
    


end











