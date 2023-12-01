function plotEndpointsFunction(LatticeGeometryStruct,BehaviorStruct,ResultsStruct, axesAmplitude)

%% inputs

%% plots

for numEndPoints = 1 : LatticeGeometryStruct.NinputANDoutput
    for numBeh = 1 : size(BehaviorStruct.forces,3)
    figure(numBeh*100+numEndPoints)
    hold on
    plot(LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes(numEndPoints),1),LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes(numEndPoints),2),'k*')
    plot(ResultsStruct.coord_deformed(LatticeGeometryStruct.outputNodes(numEndPoints),1,numBeh),ResultsStruct.coord_deformed(LatticeGeometryStruct.outputNodes(numEndPoints),2,numBeh),'b*')
    plot(BehaviorStruct.Target(numEndPoints,1,numBeh),BehaviorStruct.Target(numEndPoints,2,numBeh),'r*')
    axis([LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes(numEndPoints),1)-axesAmplitude LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes(numEndPoints),1)+axesAmplitude LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes(numEndPoints),2)-axesAmplitude LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes(numEndPoints),2)+axesAmplitude])
    title(['Point ',num2str(numEndPoints),' Behavior ',num2str(numBeh)])
    grid on
    
%     figure(200+numEndPoints)
%     hold on
%     plot(LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes(numEndPoints),1),LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes(numEndPoints),2),'k*')
%     plot(ResultsStruct.coord_deformed(LatticeGeometryStruct.outputNodes(numEndPoints),1,2),ResultsStruct.coord_deformed(LatticeGeometryStruct.outputNodes(numEndPoints),2,2),'b*')
%     plot(BehaviorStruct.Target(numEndPoints,1,2),BehaviorStruct.Target(numEndPoints,2,2),'r*')
%     axis([LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes(numEndPoints),1)-3e-6 LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes(numEndPoints),1)+3e-6 LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes(numEndPoints),2)-3e-6 LatticeGeometryStruct.coord_initial(LatticeGeometryStruct.outputNodes(numEndPoints),2)+3e-6])
%     title(['Point ',num2str(numEndPoints),' Behavior 2'])
%     grid on
    end
    
end