%This script tests the 'verify' function from NPortChargingProblem.

clear all;

rng('shuffle');

err = 1e-5;%tolerated error
dt=0.25;%integration interval

%generate continuously new test instances
while true
    nr = 5;
    nt = 5;

    maxPact = 10*rand+5;%maximum active power
    maxPapp = maxPact+20*rand;%maximum apparent power

    minCurr = 0.5*rand(nr,1);%min current to charge
    maxCurr = 50*rand(nr+nt,1)+1;%maximum current supported by devices

    chargeData.minimum = rand(nr,1);%lower bound
    chargeData.initial = chargeData.minimum+rand(nr,1);%in the beginnig
    chargeData.threshold = chargeData.initial+rand(nr,1);%required in the end
    chargeData.maximum = chargeData.threshold+rand(nr,1);%complete charge

    R = 3*diag(rand(nt+nr,1));
    
    MT = rand(nt);MR = rand(nr);
    
    ZT = R(1:nt,1:nt)-1i*(MT+MT'-2*diag(diag(MT)))/2;
    ZR = R(nt+1:end,nt+1:end)-1i*(MR+MR'-2*diag(diag(MR)))/2;
    M  = -1i*rand(nr,nt);

    Z = [ZT, M.';
         M, ZR];

    timeSlot.Z = Z;
    timeSlot.Id = rand(nr,1);
    
    rlCellArray = cell(0);
    convCellArray = cell(0);
    for i=1:nr
        %data regarding the load resistance
        rlCellArray{end+1} = [0, 0; 1, chargeData.maximum(i)];
        %current conversion function
        convCellArray{end+1} = [0, 0;
                                minCurr(i)-1e-9, 0;
                                minCurr(i), minCurr(i);
                                maxCurr(nt+i),maxCurr(nt+i)];
    end
    deviceData.rlCellArray = rlCellArray;
    deviceData.convCellArray = convCellArray;
    
    constraints.maxCurr = maxCurr;
    constraints.maxPapp = maxPapp;
    constraints.maxPact = maxPact;
    
    hashSize = ceil(10000*rand);
    nLevels = 1+ceil(10*rand);
    
    feasiblePastModel = FPDummie(hashSize, nLevels, chargeData.minimum, chargeData.maximum);
    
    tic
    target = generateTarget(feasiblePastModel, chargeData);
    disp('Generating Target: SUCCESS');
    toc
    
    times = [];%execution times for hits
    times1= [];%execution times for misses
    for i=1:length(target.bag)
        for j=1:length(target.bag{i})
            tic;
            %searching a vector that is in the target set
            q = search(target,target.bag{i}(j).charge);
            time = toc;
            if isempty(q)
                error('The search method has an error (1)');
            else
                times = [times; time];
            end
            %the vectors in target are always spaced by multiples of
            %(maximum-threshold)/(nlevels-1)
            errorVec = (chargeData.maximum-chargeData.threshold)./(target.nLevels-1);
            %searching a vector that is not in the target set
            charge = target.bag{i}(j).charge+(0.85*rand(nr,1)+0.05).*errorVec;
            q = search(target,charge);
            time = toc;
            if ~isempty(q)
                error('The search method has an error (2)');
            else
                times1 = [times1; time];
            end
        end
    end
    if isempty(times) || isempty(times1)
        error('It seems the bag is empty.');
    end
    disp('Searching for Vectors: SUCCESS');
    disp(['times for found vectors (min, average, max): ',...
        num2str(min(times)),', ',num2str(mean(times)),', ',num2str(max(times))]);
    disp(['times for not found vectors (min, average, max): ',...
        num2str(min(times1)),', ',num2str(mean(times1)),', ',num2str(max(times1))]);

end
