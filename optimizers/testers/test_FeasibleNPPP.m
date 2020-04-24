clear all;
clc;

nt = randi(4)+1;
nr = randi(4)+1;
timeLine_size = randi(10);
nSegments = 10;
sample_size = 1000;

deviceData = [];
for r = 1:nr
	%random lookup tables for load resistance, current conversion and charge conversion
	[rlTable,convTable,chargeTable] = randomLookupTables();
	%manager for the lookup tables
    deviceData = [deviceData; DeviceData(rlTable,convTable,chargeTable)];
end

success=false;

while ~success
	[success, solution, chargeData, constraints, timeLine, dt] = generateFeasibleNPPPInstance(deviceData, nt, nSegments, timeLine_size, sample_size);
end

if success
	ffModel = FeasibleFuture();  
	inst = NPortSourcingProblem(timeLine,dt,chargeData,deviceData,constraints,ffModel);
	[code,Q_list] = inst.verify([solution.V]);
	if(code~=0)
		disp(Q_list);
		disp([chargeData.initial,[solution.Q]]);
		error(['FIM: ',num2str(code)]);
	%elseif max(max(abs(Q_list-[chargeData.initial,[solution.Q]])>1e-6))>0
	%	disp(Q_list-[chargeData.initial,[solution.Q]]);
	%	error('FIM');
	end
end
