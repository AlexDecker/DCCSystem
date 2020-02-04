%Generates a random instance of the NPortSourcingProblem with a guaranteed solution "sol"
%nt: number of transmitters
%nr: number of receivers
%nSlots: number of time slots
%dt: time slot duration (s)
%maxV: maximum transmitting voltage (V)
%sampleSize: size of the sample used for the greedy voltage progression (see randomConstraints)
%sparsity and dynamicity: for controlling the timeLine (see randomTimeLine)
%ffModel: empty instance of the chosen feasible future model
function [inst,sol] = randomSourcingInstance(nt,nr,nSlots,dt,maxV,sampleSize,sparsity,dynamicity,ffModel)
    [rlTable,convTable,chargeTable,maxId,maxIn] = randomLookupTables();

    deviceData = DeviceData(rlTable,convTable,chargeTable);

    success=false;
    while ~success
        timeLine = randomTimeLine(nt,nr,nSlots,maxId*ones(nr,1),sparsity,dynamicity);
        [constraints,chargeData,sol,success] = randomConstraints(deviceData,timeLine,...
            dt,sampleSize,maxV,maxIn);
        %the failure can occour due to a single reason: the discharge current is too
        %large for the considered maxIn
        maxId = maxId/1.5;
    end 

    devList = []; 
    for r=1:nr
        devList = [devList; deviceData];
    end
    inst = NPortSourcingProblem(timeLine,dt,chargeData,devList,constraints,ffModel); 
end
