clc;
clear all;

n_files = 133;
folder = 'tunning_double2';
sufix = '_double';

found_solutions = 0;
failures = 0;
errors = 0;

for i=0:n_files-1
	load([folder,'\result', sufix, num2str(i),'.mat']);
	
	ffModel = FFPareto(result.hashSize, result.nSegments, 10000, 250, inf,...
		750, 60, result.nt, result.nr); 
		
	P = NPortSourcingProblem(result.timeLine,result.dt,result.chargeData,result.deviceData,result.constraints,ffModel);
		
	tic;
	[output.first_success, output.first_solution] = P.solve();
	output.time = toc;
	
	w = whos('P');
    output.size_MB = w.bytes/1000000;
	
	output.code = -1;
	
	output.n_iterations_list = [];
	
	if output.first_success
		
		%tic;
		%[output.second_success, output.second_solution, output.n_iterations_list] = P.recover_voltage_progression(...
		%	output.first_solution, result.max_iterations);
		%output.time_recover = toc;
		
		%if output.second_success
			%code = P.plot(solution,false);
			%[output.code, ~] = P.verify(output.second_solution.V);
			[output.code, ~] = P.verify(output.first_solution.V);
			if output.code~=0
				errors = errors + 1;
				disp(['SolveCharging: error number ',num2str(output.code)]);
			else
				found_solutions = found_solutions + 1;
				disp(['Instance #',num2str(i),' completed: success']);
			end
		%else
		%	disp(['Instance #',num2str(i),' completed: failure']);
		%	failures = failures + 1;
		%end
	else
		disp(['Instance #',num2str(i),' completed: failure']);
		failures = failures + 1;
	end
	
	save([folder,'\output_pareto', sufix, num2str(i),'.mat'], 'output');
	
	disp(['Successes: ',num2str(found_solutions),'; failures: ', num2str(failures), '; errors: ', num2str(errors)]);
	
	clear result P ffModel output;
end