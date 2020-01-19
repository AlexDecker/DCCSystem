nr = 5;
nt = 6;
nSlots = 2;
dt = 1;
maxV = 15;
sampleSize = 100;%the higher this value, the higher the difficulty

[rlTable,convTable,chargeTable,maxId,maxIn] = randomLookupTables();

deviceData = DeviceData(rlTable,convTable,chargeTable);

success=false;
while ~success
	timeLine = randomTimeLine(nt,nr,nSlots,maxId*ones(nr,1));
	[constraints,chargeData,sol,success] = randomConstraints(deviceData,timeLine,...
		dt,sampleSize,maxV,maxIn);
	%the failure can occour due to a single reason: the discharge current is too
	%large for the considered maxIn
	maxId = maxId/1.5;
end

devList = [];
for i=1:nr
	devList = [devList; deviceData];
end
P = NPortPowerProblems(timeLine,dt,chargeData,devList,constraints,FeasibleFuture());

P.verify(sol);