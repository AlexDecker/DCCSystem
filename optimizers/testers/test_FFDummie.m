clc;

rng('shuffle')

found_solutions = 0;
errors = 0;
failures = 0;
i=0;

while true

	result.nt = randi(3)+1;
	result.nr = randi(3)+1;
	result.timeLine_size = randi(10);
	result.nSegments = 5+randi(15);
	result.sample_size = 2+randi(10);
	result.max_iterations = 1000;

	%arguments for the FFDummie
	result.hashSize = 1000;
	result.maxSize = ceil((result.nSegments^result.nr)/10) + randi(result.nSegments^result.nr-ceil((result.nSegments^result.nr)/10));
	result.thr_top = 10+randi(5);
	result.thr = 10+randi(5);
	result.thr_down = 10+randi(5);
	result.ttl_top = result.maxSize;
	result.ttl = result.maxSize;
	result.ttl_down = 10+randi(5);
	result.ttl_fineAdjustment = 10+randi(5);

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

	ffModel = FFDummie(result.hashSize, result.nSegments, result.maxSize, result.thr_top, result.thr, result.thr_down, result.ttl_top,...
		result.ttl, result.ttl_down, result.ttl_fineAdjustment, result.nt, result.nr); 
		
	P = NPortSourcingProblem(result.timeLine,result.dt,result.chargeData,result.deviceData,result.constraints,ffModel);

	tic;
	[success, solution] = P.solve();
	result.time_main = toc;
	
	code = -1;
	
	result.first_success = success;
	result.first_solution = solution;
	
	n_iterations_list = [];
	
	if success
		
		tic;
		[success, solution, n_iterations_list] = P.recover_voltage_progression(solution, result.max_iterations);
		result.time_recover = toc;
		
		if success
			%code = P.plot(solution,false);
			[code, ~] = P.verify([solution.V]);
			if code~=0
				errors = errors + 1;
				disp(['SolveCharging: error number ',num2str(code)]);
			else
				found_solutions = found_solutions + 1;
				disp(['Instance #',num2str(i),' completed: success']);
			end
		else
			disp(['Instance #',num2str(i),' completed: failure']);
			failures = failures + 1;
		end
	else
		disp(['Instance #',num2str(i),' completed: failure']);
		failures = failures + 1;
	end
	
	result.success = success;
	result.solution = solution;
	result.n_iterations_list = n_iterations_list;
	result.code = code;
	
	save(['result',num2str(i),'.mat'], 'result');
	
	clearvars -except found_solutions failures errors i
	
	disp(['Successes: ',num2str(found_solutions),'; failures: ', num2str(failures), '; errors: ', num2str(errors)]);
	i=i+1;
end
