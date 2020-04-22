clear all;

%random lookup tables for load resistance, current conversion and charge conversion
[rlTable,convTable,chargeTable] = randomLookupTables();

nt = randi(4)+1;
nr = randi(4)+1;
timeLine_size = randi(10);
nSegments = 10;
sample_size = 1000;

deviceData = [];
for r = 1:nr
	%manager for the lookup tables
    deviceData = [deviceData; DeviceData(rlTable,convTable,chargeTable)];
end

[success, solution, chargeData, constraints, timeLine, dt] = generateFeasibleNPPPInstance(deviceData, nt, nSegments, timeLine_size, sample_size);

if success
	ffModel = FeasibleFuture();  
	inst = NPortSourcingProblem(timeLine,dt,chargeData,deviceData,constraints,ffModel);
end

%Q = zeros(nr,timeLine_size);
%Q(:,1) = chargeData.initial;
%for time = 1:timeLine_size
	%Q(:,time+1) = solution(time).Q;
%end

%figure;
%plot(Q(1,:),Q(2,:));
