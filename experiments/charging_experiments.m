clc;

rng('shuffle')

found_solutions = zeros(5,1);
errors = zeros(5,1);
failures = zeros(5,1);
i=1;

while true

	result.nt = randi(3)+1;
	result.nr = randi(3)+1;
	result.timeLine_size = randi(10);
	result.nSegments = 5+randi(15);
	result.sample_size = 2+randi(10);

	result.deviceData = [];
	for r = 1:result.nr
		%random lookup tables for load resistance, current conversion and charge conversion
		[rlTable,convTable,chargeTable] = randomLookupTables();
		%manager for the lookup tables
		result.deviceData = [result.deviceData; DeviceData(rlTable,convTable,chargeTable)];
	end

	success=false;

	while ~success
		[success, solution, result.chargeData, result.constraints, result.timeLine, result.dt] = generateFeasibleNPPPInstance(...
			result.deviceData, result.nt, result.nSegments, result.timeLine_size, result.sample_size);
	end

	%1:FFDummie
	%2:FFDoubleDummie
	%3:FFPareto
	%4:FFGreedy
	%5:FFRobin
	algs = [];
	algs = [algs, struct('name','FFDummie','ffModel',...
		FFDummie(1000, result.nSegments, 10000, 100, inf, inf, 250, 40, 5, 10, result.nt, result.nr))];
	algs = [algs, struct('name','FFDoubleDummie','ffModel',...
		FFDoubleDummie(1000, result.nSegments, 10000, 100, inf, inf, 250, 80, 5, result.nt, result.nr))];
	algs = [algs, struct('name','FFPareto','ffModel',...
		FFPareto(1000, result.nSegments, 10000, 250, inf, 750, 60, result.nt, result.nr))];
	algs = [algs, struct('name','FFGreedy','ffModel',...
		FFGreedy(1000, result.nSegments, 10000, 10000, result.nt, result.nr))];
	algs = [algs, struct('name','FFRobin','ffModel',...
		FFRobin(1000, result.nSegments, 10000, 10000, result.nt, result.nr))];
	
	result.ret = [];
	
	for f = 1:5
		
		ret.name = algs(f).name;
	
		P = NPortChargingProblem(result.timeLine,result.dt,result.chargeData,result.deviceData,result.constraints,algs(f).ffModel);

		tic;
		[ret.success, ret.solution] = P.solve();
		ret.time = toc;
		disp(['Time: ',num2str(ret.time),'s']);
	
		w = whos('P');
		ret.size_MB = w.bytes/1000000;
	
		ret.code = -1;
	
		if ret.success
			[ret.code, ~] = P.verify(ret.solution.V);
			if ret.code~=0
				errors(f) = errors(f) + 1;
				disp([algs(f).name,': error number ',num2str(ret.code)]);
			else
				found_solutions(f) = found_solutions(f) + 1;
				disp([algs(f).name,': success']);
			end
		else
			disp([algs(f).name,': failure']);
			failures(f) = failures(f) + 1;
		end
		
		result.ret = [result.ret, ret];
	end
	
	for f=1:5
		disp([algs(f).name,'> Successes: ',num2str(found_solutions(f)),'; failures: ', num2str(failures(f)), '; errors: ', num2str(errors(f))]);
	end
	
	save(['exp_charging_',num2str(i),'.mat'], 'result');
	
	i=i+1;
end
