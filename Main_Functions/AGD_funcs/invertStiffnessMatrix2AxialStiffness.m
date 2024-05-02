function [K_contribution] = invertStiffnessMatrix2AxialStiffness(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct)


% author: Pietro Sainaghi
% using architecture written by Erwin Muelder and Ryan Lee

% this function evaluates the matrix terms needed to invert a given resulting MNN stiffness matrix
% to compute axial stiffness values for each beam

%% function inputs

% lattice structure - LatticeGeometryStruct
Nbeams = LatticeGeometryStruct.Nbeams; % [1] number of beams in lattice
coord_initial = LatticeGeometryStruct.coord_initial; % [Ncoord,3] x,y,z coordinates of each node 
DOI = LatticeGeometryStruct.DOI; % [1] dimensions in space = 3
Ncoord = LatticeGeometryStruct.Ncoord;

% behavior structure - BehaviorStruct
Ncases = BehaviorStruct.Ncases;

% finite truss structure - FEMStruct
DOFnodes = FEMStruct.DOFnodes; % [] indices of nodes that are not grounded
DOFFinal = FEMStruct.DOFFinal;
F = FEMStruct.F; % [NDOF, Ncases] location and value of forces within degrees of freedom format for each behavior
k_base = FEMStruct.k_base; % [6,6] stiffness matrix of individual element
R6 = FEMStruct.R6;
RT6 = FEMStruct.RT6;
DOF = FEMStruct.DOF; % [1] number of degrees of freedom in system
Degrees_per_element = FEMStruct.Degrees_per_element;


% optimization problem structure - OptimizerDataStruct
% ForceScaling = OptimizerDataStruct.ForceScaling; % [boolean] yes or no to force scaling

%% function outputs

% K_contribution [DOF,DOF,Nbeams] for each tunable beam, the imaginary terms show what entries are changed by the tuning, used to contruct an invertible matrix


%% create stiffness contribution matrices

% generates a 

% store contribution to stiffness matrix from each individual axial stiffness
K_contribution=zeros(DOF,DOF,Nbeams);
% stiffness matrix for individual beam with control stiffness value
k_beam=(zeros(2*DOI,2*DOI,Nbeams));

for e=1:Nbeams
    %scale the stiffnes accoding to x
    k_beam(:,:,e)=k_base(:,:,e)*diag([sqrt(-1) 1 1 sqrt(-1) 1 1]); % USE mtimesx?
    %Combine stiffness matrix for solving
    K_contribution(Degrees_per_element(e,:),Degrees_per_element(e,:),e)=RT6(:,:,e)*k_beam(:,:,e)*R6(:,:,e);
end

%% identify indices contribution for each 






end


