%some experiments failed due an error from recover_voltage_progression. However, the result files
%store recover_voltage_progression input. This script aims fixing the data with the new version of
%recover_voltage_progression

clc;
clear all;

n_files = 420;

folder = 'tunning4';

for i=0:n_files-1
	load([folder,'\result', num2str(i),'.mat']);
	
	if result.code~=0 && result.code~=-1
		ffModel = FFDummie(result.hashSize, result.nSegments, result.maxSize, result.thr_top, result.thr, result.thr_down, result.ttl_top,...
			result.ttl, result.ttl_down, result.ttl_fineAdjustment, result.nt, result.nr); 
		P = NPortSourcingProblem(result.timeLine,result.dt,result.chargeData,result.deviceData,result.constraints,ffModel);
		
		Q = result.first_solution.Q;

		if min(Q(:,end)>=result.chargeData.threshold)==0
			disp(['Arquivo ',num2str(i),' falha']);
		else
			[success, solution, ~] = P.recover_voltage_progression(result.first_solution, 1000);%result.max_iterations);
			if success
				[code, QList] = P.verify(solution.V);
				if code~=0
					QList
					disp('solution.Q = ');
					solution.Q
					error(['Arquivo ',num2str(i),' da erro ', num2str(code)]);
				else
					disp(['Arquivo ',num2str(i),' da certo']);
				end
			else
				disp(['Arquivo ',num2str(i),' falha do recover']);
			end
		end
	end
	clear result P ffModel;
end