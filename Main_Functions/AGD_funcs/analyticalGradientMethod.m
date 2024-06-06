function [x] = analyticalGradientMethod(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct, xInit)

% Algorithm by Jiaji Chen to work with mass-spring model, rewritten for finite truss model by Pietro Sainaghi

% this function uses an analytical gradient descent method to optimize a MNN



%% input unpacking

% optimizer structure - OptimizerDataStruct
% iterMAX = OptimizerDataStruct.AGDmaxIterations;

learningRate = 1; % Learning Rate
Decay_interval = 400; % Learning rate decay interval
Decay_factor = 1; % Learning rate decay factor
% Batch_num = 1; % Minibatch size UNUSED but shows up in original algorithm
MSEscale = 1; % scaling term for MSE used in algorthm

% termination conditions
iterMAX = 200000; % maximum iterations
MSEChangethreshold = 1e-7; % minimum change in MSE
MSEthreshold = 0.001; % target MSE
averageWindow = 400;
averageWindowFrequency = 5000;
averageWindowThreshold = 1e-4;


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
DOI = LatticeGeometryStruct.DOI; % [1] dimensions in plane, corresponding to x, y, thetaz
Nbeams = LatticeGeometryStruct.Nbeams; % [1] number of tunable beams

% finite truss structure - FEMStruct
F_all = FEMStruct.F; % [NDOF, Ncases] location and value of forces within degrees of freedom format for each behavior
DOFnodes = FEMStruct.DOFnodes; % [] indices of nodes that are not grounded

% beam properties structure - LinkPropertiesStruct
kLinMax = LinkPropertiesStruct.kLinMax;
kLinMin = LinkPropertiesStruct.kLinMin;


%% Initialize length parameters

% NDOFnodes [] number of nodes that are not grounded
NDOFnodes = length(DOFnodes);
% NDOF [] numebr of degrees of freedom
NDOF = NDOFnodes*DOI;



%% generate K to x inversion function

% initialize stiffness matrix derivative
dKdx_matrix = zeros(NDOF,NDOF,length(xInit));

[stiffnessMatrixDerivative] = computeStiffnessMatrixDerivative(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct);
if size(stiffnessMatrixDerivative) == size(dKdx_matrix)
    dKdx_matrix = stiffnessMatrixDerivative;
else
    error('Issue with computing derivative of stiffness wrt tunable beams')
end

% repackage in matrix form
dKdx = zeros(length(xInit),NDOF^2);
for xIDX = 1:size(dKdx_matrix,3)
    for KIDX = 1:size(dKdx_matrix,1)
        dKdx(xIDX, NDOF*(KIDX-1)+1 : NDOF*(KIDX) ) = dKdx_matrix(KIDX,:,xIDX);
    end
end


%% contruct relationship matrix between DOF and output matrices

% identify indices in NDOFnodes format corresponding to output nodes
outINDOF_indices = [];
for outNodeIDX = 1:length(Outputline)
    for dofNodeIDX = 1:length(DOFnodes)
        if Outputline(outNodeIDX) == DOFnodes(dofNodeIDX)
            outINDOF_indices = [outINDOF_indices dofNodeIDX];
        end
    end
end


%% initialize time history storage 

% list of loss function value at each iteration
MSE_hist = [1]; 

% history of stiffness combinations
Stiffness_hist = [];


%% Optmizer Loop

% assign initial condition on stiffness
x = xInit; % column vector

% iteration counter
iter = 1;
% start change in error
MSEchange = 10;

figure('Name','AGD')

while MSEchange > MSEChangethreshold

    % update iteration
    iter = iter + 1;

    % store stiffness combination for each iteration
    Stiffness_hist = [Stiffness_hist x];

    % compute loss function (scaled MSE)
    % works with both unitless MSE and unit MSE, but may perform differently based on MSEscale
    MSE = ERROR_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
    L = MSEscale^2/2 * MSE;
    % disp(['MSE: ', num2str(MSE)])

    % store loss in vector to plot
    MSE_hist = [MSE_hist MSE];

    % compute change in MSE
    MSEchange = abs(MSE_hist(end-1) - MSE_hist(end));
    if MSEchange < MSEChangethreshold
        disp('Termination: Reached error change threshold')
    end



    % extract FEM properties used in algorithm
    [coorddeformed,Kdof,Fdof,udof]=FEM(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
    F = Fdof;
    u = udof;

    % repackage stiffness matrix into vector
    for KIDX = 1:size(dKdx_matrix,1)
        K(NDOF*(KIDX-1)+1 : NDOF*(KIDX) ,1) = Kdof(KIDX,:);
    end

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

    % restructure output node displacements and target displacements
    u_out = zeros(NinputANDoutput*size(Target,2),Ncases);
    t_out = zeros(NinputANDoutput*size(Target,2),Ncases);
    for behIDX = 1:Ncases
        for UxIDX = 1:length(Outputline)

            % actual displacements
            u_out((UxIDX-1)*2+1,behIDX) = udof( outINDOF_indices(UxIDX)*3 - 2, behIDX );
            u_out((UxIDX-1)*2+2,behIDX) = udof( outINDOF_indices(UxIDX)*3 - 1, behIDX );

            % target displacements
            t_out((UxIDX-1)*2+1,behIDX) = Target(UxIDX,1,behIDX) - coord_initial(Outputline(UxIDX),1);
            t_out((UxIDX-1)*2+2,behIDX) = Target(UxIDX,2,behIDX) - coord_initial(Outputline(UxIDX),2);

            % debugging
            if u_out((UxIDX-1)*2+1,behIDX) ~= u_all(Outputline(UxIDX),1,behIDX)
                error('Mismatch in displacement vector x')
            end
            if u_out((UxIDX-1)*2+2,behIDX) ~= u_all(Outputline(UxIDX),2,behIDX)
                error('Mismatch in displacement vector x')
            end
        end
    end

    % compute derivative of cost function with respect to output displacements
    dLdu = zeros(NinputANDoutput*size(Target,2),Ncases);
    for behIDX = 1:Ncases
        for outDOFIDX = 1:(NinputANDoutput*size(Target,2))
            dLdu(outDOFIDX,behIDX) = MSEscale^2 * ( u_out(outDOFIDX,behIDX) - t_out(outDOFIDX,behIDX) );
        end
    end

    % compute derivative of force wrt displacement
    dFdu = zeros(NDOF,NDOF);
    % dFdu = Kdof'; % TEST
    dFdu = Kdof;

    % compute derivative of force wrt stiffness
    dFdK = zeros(NDOF,NDOF^2,Ncases);
    for behIDX = 1:Ncases
        for dofIDX = 1:(NDOF)
            for kIDX = 1:(NDOF^2)
                % identity of entry is nonzero
                rowdivisor = ceil(kIDX/NDOF);
                if rowdivisor == dofIDX
                    % select specific entry in u
                    nonzeroIDX = mod(kIDX,NDOF);
                    if nonzeroIDX == 0
                        nonzeroIDX = 24;
                    end
                    selectionVec = zeros(1,NDOF);
                    selectionVec(nonzeroIDX) = 1;
                    % compute nonzero entry
                    dFdK(dofIDX,kIDX,behIDX) = u(nonzeroIDX,behIDX);
                end
            end
        end
    end

    % compute derivative of displacement with respect to stiffness
    dudK = zeros(NDOF,NDOF^2,Ncases);
    for behIDX = 1:Ncases
        dudK(:,:,behIDX) = dFdu\-dFdK(:,:,behIDX);
    end

    % compute derivative of output displacements with respect to stiffnesssi
    duoutdK = zeros(NinputANDoutput*size(Target,2),NDOF^2,Ncases);
    for outIDX = 1:length(outINDOF_indices)
        duoutdK((outIDX-1)*2+1,:) = dudK(outINDOF_indices(outIDX)*3 - 2,:);
        duoutdK((outIDX-1)*2+2,:) = dudK(outINDOF_indices(outIDX)*3 - 1,:);
    end

    % compute derivative of cost function with respect to stiffness
    dLdK_store = zeros(NinputANDoutput*size(Target,2),NDOF^2,Ncases);
    dLdK = zeros(NDOF^2,1,Ncases);
    for behIDX = 1:Ncases
        dLdK_store(:,:,behIDX) = duoutdK(:,:,behIDX).*dLdu(:,behIDX);
        dLdK(:,:,behIDX) = sum(dLdK_store(:,:,behIDX))';
    end

    % compute derivative of cost function with respect to beam value
    dLdx = zeros(length(xInit),1,Ncases);
    for behIDX = 1:Ncases
        dLdx(:,:,behIDX) = dKdx*dLdK(:,:,behIDX);
    end

    % compute unit scaling term
    largestterm = max(nonzeros(dLdx));
    smallestterm = min(nonzeros(dLdx));
    unitscaler = 1/largestterm;

    % compute gradient descent
    for behIDX = 1:Ncases
        x = x - unitscaler * learningRate * dLdx(:,:,behIDX);
    end

    % enforce maximum and minimum stiffness
    for xIDX = 1:length(x)
        if x(xIDX) > kLinMax
            x(xIDX) = kLinMax*0.9;
        elseif x(xIDX) < kLinMin
            x(xIDX) = kLinMin*0.9;
        end
    end

    
    plot(MSE_hist)
    xlabel('Iteration')
    ylabel('MSE')
    title([num2str(MSE)])
    drawnow

    if iter > iterMAX
        MSEchange = 0;
        disp('Termination: Maximum Number of Iterations')
    end

    if MSE < MSEthreshold
        MSEchange = 0;
        disp('Termination: Target MSE reached')
    end

    if mod(iter, averageWindowFrequency) == 0
        valueWindows = MSE_hist((end-averageWindow):end);
        MSEAVwindchange = max([abs(max(valueWindows)-mean(valueWindows)) abs(min(valueWindows)-mean(valueWindows))]);
        if MSEAVwindchange < averageWindowThreshold
            MSEchange = 0;
            disp('Termination: Variation in Iteration Window is less than threshold')
        end
    end


end

disp('Finished Optimization')









