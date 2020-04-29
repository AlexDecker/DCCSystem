clc;

rng('shuffle')

found_solutions = zeros(5,1);
errors = zeros(5,1);
failures = zeros(5,1);
i=0;

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
	ffModel = FFDummie(1000, result.nSegments, 10000, 5000, inf, inf, 250, 20, 20, 50, result.nt, result.nr);
	ffModel = [ffModel, FFDoubleDummie(1000, result.nSegments, 10000, 5000, inf, inf, 250, 20, 20, result.nt, result.nr)]; 
	ffModel = [ffModel, FFPareto(1000, result.nSegments, 10000, 250, inf, 750, 60, result.nt, result.nr)];
	ffModel = [ffModel, FFGreedy(1000, result.nSegments, 10000, 10000, result.nt, result.nr)];
	ffModel = [ffModel, FFRobin(1000, result.nSegments, 10000, 10000, result.nt, result.nr)];
	
	names = [];
	names = [names, struct('txt','FFDummie')];names = [names, struct('txt','FFDoubleDummie')];
	names = [names, struct('txt','FFPareto')];names = [names, struct('txt','FFGreedy')];
	names = [names, struct('txt','FFRobin')];
	
	result.ret = [];
	
	for f = 1:5
		P = NPortSourcingProblem(result.timeLine,result.dt,result.chargeData,result.deviceData,result.constraints,ffModel(f));

		tic;
		[ret.success, ret.solution] = P.solve();
		ret.time = toc;
	
		w = whos('P');
		ret.size_MB = w.bytes/1000000;
	
		ret.code = -1;
	
		if ret.success
			[ret.code, ~] = P.verify(ret.solution.V);
			if ret.code~=0
				errors(f) = errors(f) + 1;
				disp([names(f).txt,': error number ',num2str(ret.code)]);
			else
				found_solutions(f) = found_solutions(f) + 1;
				disp([names(f).txt,': success']);
			end
		else
			disp([names(f).txt,': failure']);
			failures(f) = failures(f) + 1;
		end
		
		result.ret = [result.ret, ret];
		disp(['Successes: ',num2str(found_solutions(f)),'; failures: ', num2str(failures(f)), '; errors: ', num2str(errors(f))]);
	end
	
	save(['exp_',num2str(i),'.mat'], 'result');
	
	i=i+1;
end
