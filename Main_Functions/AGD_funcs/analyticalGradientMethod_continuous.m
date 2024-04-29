function [x] = analyticalGradientMethod_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct)

% Algorithm by Jiaji Chen

%% input unpacking


% extract optimization parameters
Iter_max = OptimizerDataStruct.AGDmaxIterations;
eta = 0.01; % Learning Rate
Decay_interval = 50; % Learning rate decay interval
Decay_factor = 5; % Learning rate decay factor
Batch_num = 1; % Minibatch size
scale = 20;
dL  = 0.3;






%% Training
Loss_hist = [];

for iter = 1:Iter_max
    if mod(iter, Decay_interval == 0)
        eta = eta/Decay_factor;
    end

    % compute loss function (MSE multiplied by scale^2 in our case)
    [error]=ERROR_continuous(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
    Loss = 0.5 * scale^2 * error;

    % compute loss function derivative (absolute error multiplied by scale^2 in our case)
    [coorddeformed,Kf2u,Ku2f]=FEM(LinkPropertiesStruct, LatticeGeometryStruct, BehaviorStruct,FEMStruct,OptimizerDataStruct,x);
    dLdu = scale^2 * (Output - Target);

    % Plot current configuration
    FreePosition = u + Nodes;
    figure(2)
    xlim([-1, sqrt(L0^2 - (L0/2)^2)*Nx*2+1])
    ylim([-2, Ny*L0+1])
    pbaspect([Nx Ny 1])
    hold on
    gplot(Edges, FreePosition);
    plot(Target_x, Target_y)
    hold off


    % % Initialize dFdU and dFdk
    % dFdu = zeros(N_Free * 2, N_Free * 2);
    % dFdk = zeros(N_Free * 2, N_k);
    % 
    % for i = 1:N_Free
    %     Ind1 = NCN_Free(i);
    %     for Ind2 = 1:N_nodes
    % 
    %         if Edges(Ind1,Ind2) == 1
    %             ii = 2 * i - 1;
    %             kk = K(Ind1, Ind2);
    % 
    %             [dFx_dx1, dFx_dx2, dFx_dy1, dFx_dy2, dFy_dx1, dFy_dx2, dFy_dy1, dFy_dy2, dFx_dk, dFy_dk] ...
    %                 = force_and_derivative(K_array(kk), L0, u(Ind1,:), u(Ind2,:), Nodes(Ind1,:), Nodes(Ind2,:));
    % 
    %             dFdu(ii, ii) = dFdu(ii, ii) + dFx_dx1;
    %             dFdu(ii, ii + 1) = dFdu(ii, ii + 1) + dFx_dy1;
    %             dFdu(ii + 1, ii) =  dFdu(ii + 1, ii) + dFy_dx1;
    %             dFdu(ii + 1, ii + 1) = dFdu(ii + 1, ii + 1) + dFy_dy1;
    % 
    %             if ismember(Ind2, NCN_Free)
    %                 j = find(NCN_Free == Ind2);
    %                 jj = 2 * j -1;
    %                 dFdu(ii, jj) = dFdu(ii, jj) + dFx_dx1;
    %                 dFdu(ii, jj + 1) = dFdu(ii, jj + 1) + dFx_dy1;
    %                 dFdu(ii + 1, jj) = dFdu(ii + 1, jj) + dFy_dx1;
    %                 dFdu(ii + 1, jj + 1) =  dFdu(ii + 1, jj + 1) + dFy_dy1;
    %             end
    % 
    %             dFdk(ii, kk) =  dFdk(ii, kk) + dFx_dk;
    %             dFdk(ii + 1, kk) = dFdk(ii + 1, kk) + dFy_dk;
    %         end
    % 
    %     end
    % end
    % 
    % OutputNodesAlt_x = [];
    % for i = 1:length(OutputNodes)
    %     outputnodes = find(NCN_Free == OutputNodes(i));
    %     OutputNodesAlt_x = [OutputNodesAlt_x; 2 * outputnodes - 1];
    % end
    % OutputNodesAlt_y = OutputNodesAlt_x + 1;
    % 
    % 
    % dudk = dFdu\-dFdk;
    % dudk_output = dudk([   OutputNodesAlt_y     ], :);
    % dLdk = dudk_output .* dLdu;
    % dLdk = sum(dLdk)';
    % 
    % % Update K_array
    % K_array = K_array - eta * dLdk;

    dFdU = zeros(N_Free*2, N_Free*2);
    dFdK = zeros(N_Free*2, N_k);

    for i = 1:N_Free
        Ind1 = NCN_Free(i);
        for Ind2 = 1:N_nodes
            if Edges(Ind1,Ind2) == 1
                [dfx_du1x, dfx_du1y, dfx_du2x, dfx_du2y, dfy_du1x, dfy_du1y, dfy_du2x, dfy_du2y, dfx_dk, dfy_dk]...
                    = force_and_derivative_backup(u(Ind1,:), u(Ind2,:), K(Ind1,Ind2), K_array, L0, Nodes(Ind1,:), Nodes(Ind2,:));
                dFdU(i,i) = dFdU(i,i) + dfx_du1x;
                dFdU(i,i + N_Free) = dFdU(i,i + N_Free) + dfx_du1y;
                dFdU(i + N_Free,i) = dFdU(i + N_Free,i) + dfy_du1x;
                dFdU(i + N_Free,i + N_Free) = dFdU(i + N_Free,i + N_Free) + dfy_du1y;

                dFdK(i,K(Ind1,Ind2)) = dFdK(i,K(Ind1,Ind2)) + dfx_dk;
                dFdK(i + N_Free,K(Ind1,Ind2)) = dFdK(i,K(Ind1,Ind2)) + dfy_dk;

                if ismember(Ind2, NCN_Free)
                    j = find(NCN_Free == Ind2);
                    dFdU(i,j) = dFdU(i,j) + dfx_du2x;
                    dFdU(i,j + N_Free) = dFdU(i,j + N_Free) + dfx_du2y;
                    dFdU(i + N_Free,j) = dFdU(i + N_Free,j) + dfy_du2x;
                    dFdU(i + N_Free,j + N_Free) = dFdU(i + N_Free,j + N_Free) + dfy_du2y;
                end
            end
        end
    end

    dUdK = dFdU\-dFdK;

    OutputNodesAlt = [];
    for i = 1:length(OutputNodes)
        outputnodes = find(NCN_Free == OutputNodes(i));
        OutputNodesAlt = [OutputNodesAlt; outputnodes];

    end

    dUdK_output = dUdK([OutputNodesAlt+N_Free],:);
    dLdK = dUdK_output .* dLdu;
    dLdK = sum(dLdK)';

    K_array = K_array - eta * dLdK;


    clc
    fprintf('Loss Function: %4.3f\n', Loss)
    fprintf('Iteration: %4.3f\n', iter)

end










