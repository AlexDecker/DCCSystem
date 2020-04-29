clc;
clear all;

n_files = 133;
folder = 'tunning_double2';
sufix = '_double';

greedy_found_solutions = 0;
greedy_failures = 0;
greedy_errors = 0;
robin_found_solutions = 0;
robin_failures = 0;
robin_errors = 0;

for i=0:n_files-1
	load([folder,'\result', sufix, num2str(i),'.mat']);
	
	%GREEDY-------------------------------------------------------------------------------------------------------------------------------
	ffModel_greedy = FFGreedy(result.hashSize, result.nSegments, result.maxSize, 10000, result.nt, result.nr);
		
	P = NPortSourcingProblem(result.timeLine,result.dt,result.chargeData,result.deviceData,result.constraints,ffModel_greedy);
		
	tic;
	[output.greedy.first_success, output.greedy.first_solution] = P.solve();
	output.greedy.time = toc;
	
	w = whos('P');
    output.greedy.size_MB = w.bytes/1000000;
	
	output.greedy.code = -1;
	
	output.greedy.n_iterations_list = [];
	
	if output.greedy.first_success
		
		tic;
		[output.greedy.second_success, output.greedy.second_solution, output.greedy.n_iterations_list] = P.recover_voltage_progression(...
			output.greedy.first_solution, result.max_iterations);
		output.greedy.time_recover = toc;
		
		if output.greedy.second_success
			%code = P.plot(solution,false);
			[output.greedy.code, ~] = P.verify(output.greedy.second_solution.V);
			if output.greedy.code~=0
				greedy_errors = greedy_errors + 1;
				disp(['Greedy: SolveCharging: error number ',num2str(output.greedy.code)]);
			else
				greedy_found_solutions = greedy_found_solutions + 1;
				disp(['Greedy Instance #',num2str(i),' completed: success']);
			end
		else
			disp(['Greedy Instance #',num2str(i),' completed: failure']);
			greedy_failures = greedy_failures + 1;
		end
	else
		disp(['Greedy Instance #',num2str(i),' completed: failure']);
		greedy_failures = greedy_failures + 1;
	end
	
	clear P ffModel;
	
	%ROBIN-------------------------------------------------------------------------------------------------------------------------------
	ffModel_robin = FFRobin(result.hashSize, result.nSegments, result.maxSize, 10000, result.nt, result.nr);
		
	P = NPortSourcingProblem(result.timeLine,result.dt,result.chargeData,result.deviceData,result.constraints,ffModel_robin);
		
	tic;
	[output.robin.first_success, output.robin.first_solution] = P.solve();
	output.robin.time = toc;
	
	w = whos('P');
    output.robin.size_MB = w.bytes/1000000;
	
	output.robin.code = -1;
	
	output.robin.n_iterations_list = [];
	
	if output.robin.first_success
		
		tic;
		[output.robin.second_success, output.robin.second_solution, output.robin.n_iterations_list] = P.recover_voltage_progression(...
			output.robin.first_solution, result.max_iterations);
		output.robin.time_recover = toc;
		
		if output.robin.second_success
			%code = P.plot(solution,false);
			[output.robin.code, ~] = P.verify(output.robin.second_solution.V);
			if output.robin.code~=0
				robin_errors = robin_errors + 1;
				disp(['Robin: SolveCharging: error number ',num2str(output.robin.code)]);
			else
				robin_found_solutions = robin_found_solutions + 1;
				disp(['Robin Instance #',num2str(i),' completed: success']);
			end
		else
			disp(['Robin Instance #',num2str(i),' completed: failure']);
			robin_failures = robin_failures + 1;
		end
	else
		disp(['Robin Instance #',num2str(i),' completed: failure']);
		robin_failures = robin_failures + 1;
	end
	
	save([folder,'\output', sufix, num2str(i),'.mat'], 'output');
	
	disp(['Greedy: Successes: ',num2str(greedy_found_solutions),'; failures: ', num2str(greedy_failures), '; errors: ', num2str(greedy_errors)]);
	disp(['Robin: Successes: ',num2str(robin_found_solutions),'; failures: ', num2str(robin_failures), '; errors: ', num2str(robin_errors)]);
	
	clear result P ffModel output;
end