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
	result.nSegments = 5;
	result.sample_size = 100;
	
	%Li_Ion_Battery_LIR18650
	rlTable = ...
	[0.0, 0.0;
	0.004134969, 2.81087123;
	0.006134969, 3.10877193;
	0.007668712, 3.196992481;
	0.010202454, 3.245112782;
	0.024539877, 3.393483709;
	0.032208589, 3.445614035;
	0.044478528, 3.521804511;
	0.059815951, 3.581954887;
	0.081288344, 3.610025063;
	0.098159509, 3.618045113;
	0.116564417, 3.626065163;
	0.147239264, 3.634085213;
	0.170245399, 3.654135338;
	0.260736196, 3.690225564;
	0.346625767, 3.714285714;
	0.401840491, 3.726315789;
	0.469325153, 3.746365915;
	0.535276074, 3.77443609;
	0.630368098, 3.814536341;
	0.679447853, 3.842606516;
	0.743865031, 3.894736842;
	0.771472393, 3.910776942;
	0.851226994, 4.0182;
	0.8828, 4.0608;
	0.9379, 4.1347;
	0.9621, 4.1659;
	0.9759, 4.1815;
	0.9828, 4.1858;
	1.0000, 4.2000];
	
	convTable = ...
	[0         0;
    0.0300    0.0180;
    0.0900    0.0675;
    0.1500    0.1226;
    0.2000    0.1680;
    0.2400    0.2040;
    0.2800    0.2366;
    0.3500    0.3421;
    1.2000    1.1250];
	
	
	chargeTable = ...
	[-0.5, -0.64;
	0, 0;
	0.094, 0;
	0.164, 0.04;
	0.995, 0.95;
	1.024, 0.977;
	1.057, 0.993;
	1.1, 1;
	1.125, 1];
	
	result.deviceData = [];
	for r = 1:result.nr
		%manager for the lookup tables
		result.deviceData = [result.deviceData; DeviceData(rlTable,convTable,chargeTable)];
	end

	success=false;

	while ~success
		[success, solution, result.chargeData, result.constraints, result.timeLine, result.dt] = generateFeasibleNPPPInstance(...
			result.deviceData, result.nt, result.nSegments, result.timeLine_size, result.sample_size);
	end

	algs = [];
	algs = [algs, struct('name','FFGreedyCurrent','ffModel',...
		FFGreedyCurrent(1000, result.nSegments, 1000, 1000, result.nt, result.nr))];
	algs = [algs, struct('name','FFGreedy','ffModel',...
		FFGreedy(1000, result.nSegments, 1000, 1000, result.nt, result.nr))];
	algs = [algs, struct('name','FFRobin','ffModel',...
		FFRobin(1000, result.nSegments, 1000, 1000, result.nt, result.nr))];
	
	result.ret = [];
	
	for f = 1:3
		
		ret.name = algs(f).name;
	
		P = NPortSourcingProblem(result.timeLine,result.dt,result.chargeData,result.deviceData,result.constraints,algs(f).ffModel);

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
	
	for f=1:3
		disp([algs(f).name,'> Successes: ',num2str(found_solutions(f)),'; failures: ', num2str(failures(f)), '; errors: ', num2str(errors(f))]);
	end
	
	save(['exp_',num2str(i),'.mat'], 'result');
	
	i=i+1;
end
