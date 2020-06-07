%similar to parser, but generate the data online from a feasible future model and the instaces

function [x_nr_unique,y_nr_mean,y_nr_err, x_nt_unique,y_nt_mean,y_nt_err, x_t_unique,y_t_mean,y_t_err,success_vector] = solve_and_plot_charging(nFiles, folder, prefix, alg, use_vector)

	x_nr=[];
	y_nr=[];
	
	x_nt=[];
	y_nt=[];
	
	x_t=[];
	y_t=[];
	
	total = 0;
	total_success = 0;
	success_vector = [];
	
	for i=1:nFiles-1
		total = total+1;
		
		if mod(i,10)==0
			disp(['starting instance ',num2str(i)]);
		end
		load([folder, '\', prefix, num2str(i),'.mat']);
		
		if alg == 3
			ffModel = FFGreedyCurrent(1, result.nSegments, 1, 1000, result.nt, result.nr);
		elseif alg == 4
			ffModel = FFMaxPower(1, result.nSegments, 1, 1000, result.nt, result.nr);
		elseif alg==5
			ffModel = FFMultiSpot(1, result.nSegments, 1, result.nt, result.nr);
		else
			error('invalid alg');
		end
		
		P = NPortChargingProblem(result.timeLine,result.dt,result.chargeData,result.deviceData,result.constraints,ffModel);

		[success, solution] = P.solve();
		
		[code, ~] = P.verify(solution.V);
		
		success_vector = [success_vector,success && code==0];
		
		if success && code==0 && use_vector(i)
			len = length(result.timeLine);
			[~,sol] = size(solution.Q);
			
			x_nr=[x_nr;result.nr];
			y_nr=[y_nr;sol/len];
			
			x_nt=[x_nt;result.nt];
			y_nt=[y_nt;sol/len];
			
			x_t=[x_t;result.timeLine_size];
			y_t=[y_t;sol/len];
			
			total_success = total_success+1;
		end
		
		clear result;
	end

	p = kruskalwallis(y_nr,x_nr,'off');
	disp(['Kruskal-Wallis test (nr): p-value (H0: independency): ',num2str(p)]);
	
	p = kruskalwallis(y_nt,x_nt,'off');
	disp(['Kruskal-Wallis test (nt): p-value (H0: independency): ',num2str(p)]);
	
	p = kruskalwallis(y_t,x_t,'off');
	disp(['Kruskal-Wallis test (t): p-value (H0: independency): ',num2str(p)]);
	
	disp(['Successes: ',num2str(total_success)]);
	disp(['Total: ',num2str(total)]);
	
	x_nr_unique = unique(x_nr);
	x_nr_unique = sort(x_nr_unique);
	
	x_nt_unique = unique(x_nt);
	x_nt_unique = sort(x_nt_unique);
	
	x_t_unique = unique(x_t);
	x_t_unique = sort(x_t_unique);
	
	y_nr_mean = x_nr_unique;
	y_nr_err = x_nr_unique;
	for k=1:length(x_nr_unique)
		y_nr_mean(k) = mean(y_nr(x_nr==x_nr_unique(k)));
		t = tinv(0.9,sum(x_nr==x_nr_unique(k))-1);
		s = std(y_nr(x_nr==x_nr_unique(k)));
		sqrt_n = sqrt(sum(x_nr==x_nr_unique(k)));
		y_nr_err(k) = t*s/sqrt_n;
	end
	
	y_nt_mean = x_nt_unique;
	y_nt_err = x_nt_unique;
	for k=1:length(x_nt_unique)
		y_nt_mean(k) = mean(y_nt(x_nt==x_nt_unique(k)));
		t = tinv(0.9,sum(x_nt==x_nt_unique(k))-1);
		s = std(y_nt(x_nt==x_nt_unique(k)));
		sqrt_n = sqrt(sum(x_nt==x_nt_unique(k)));
		y_nt_err(k) = t*s/sqrt_n;
	end
	
	y_t_mean = x_t_unique;
	y_t_err = x_t_unique;
	for k=1:length(x_t_unique)
		y_t_mean(k) = mean(y_t(x_t==x_t_unique(k)));
		t = tinv(0.9,sum(x_t==x_t_unique(k))-1);
		s = std(y_t(x_t==x_t_unique(k)));
		sqrt_n = sqrt(sum(x_t==x_t_unique(k)));
		y_t_err(k) = t*s/sqrt_n;
	end
end