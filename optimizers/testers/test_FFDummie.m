clc;

nt = 2;%randi(4)+1;
nr = randi(4)+1;
timeLine_size = 1;%randi(10);
nSegments = 10;
sample_size = 10;
max_iterations = 1000;

%arguments for the FFDummie
hashSize = 1000;
maxSize = nSegments^nr;
thr_top = 1000;
thr = 100;
thr_down = 100;
ttl_top = maxSize;
ttl = maxSize;
ttl_down = 20;
ttl_fineAdjustment = 30;

found_solutions = 0;
errors = 0;
failures = 0;
i=1;

while true
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

	ffModel = FFDummie(hashSize, nSegments, maxSize, thr_top, thr, thr_down, ttl_top, ttl, ttl_down, ttl_fineAdjustment, nt, nr); 
	P = NPortSourcingProblem(timeLine,dt,chargeData,deviceData,constraints,ffModel);

	[solveable, solution] = P.solve();

	if solveable
		
		[success, solution] = P.recover_voltage_progression(solution, max_iterations);
		
		if success
			%result = P.plot(solution,false);
			[result, ~] = P.verify([solution.V]);
			if result~=0
				errors = errors + 1;
				disp(['SolveCharging: error number ',num2str(result)]);
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
	i=i+1;
end
	
disp(['Effectiveness (0-1): ', num2str(found_solutions/(found_solutions+errors+failures))])
disp(['Errors (0-1): ', num2str(errors/(found_solutions+errors+failures))])
disp(['Failures (0-1): ', num2str(failures/(found_solutions+errors+failures))])
