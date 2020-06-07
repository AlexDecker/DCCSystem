%similar to parser, but generate the data online from a feasible future model and the instaces

function [X_nr,Y_nr,ERR_nr, X_nt,Y_nt,ERR_nt, X_t,Y_t,ERR_t] = solve_and_plot(nFiles, folder, prefix, alg, marker)
	variable_nr = [];
	number_of_success_nr = [];
	number_of_occurrences_nr = [];
	
	variable_nt = [];
	number_of_success_nt = [];
	number_of_occurrences_nt = [];
	
	variable_t = [];
	number_of_success_t = [];
	number_of_occurrences_t = [];

	x_nr=[];
	y_nr=[];
	
	x_nt=[];
	y_nt=[];
	
	x_t=[];
	y_t=[];

	for i=1:nFiles-1
		disp(['starting instance ',num2str(i)]);
		load([folder, '\', prefix, num2str(i),'.mat']);
		
		if alg == 4
			ffModel = FFMaxPower(1, result.nSegments, 1, 1000, result.nt, result.nr);
		elseif alg==5
			ffModel = FFMultiSpot(1, result.nSegments, 1, result.nt, result.nr);
		else
			error('invalid alg');
		end
		
		P = NPortSourcingProblem(result.timeLine,result.dt,result.chargeData,result.deviceData,result.constraints,ffModel);

		[success, solution] = P.solve();
		
		[code, ~] = P.verify(solution.V);
		
		x_nr=[x_nr;result.nr];
		y_nr=[y_nr;success];
		
		x_nt=[x_nt;result.nt];
		y_nt=[y_nt;success];
		
		x_t=[x_t;result.timeLine_size];
		y_t=[y_t;success];
		
		ind_nr = find(variable_nr==result.nr);
		if isempty(ind_nr)
			variable_nr = [variable_nr,result.nr];
			number_of_success_nr = [number_of_success_nr,0];
			number_of_occurrences_nr = [number_of_occurrences_nr,0];
			ind_nr = length(variable_nr);
		end
		
		ind_nt = find(variable_nt==result.nt);
		if isempty(ind_nt)
			variable_nt = [variable_nt,result.nt];
			number_of_success_nt = [number_of_success_nt,0];
			number_of_occurrences_nt = [number_of_occurrences_nt,0];
			ind_nt = length(variable_nt);
		end
		
		ind_t = find(variable_t==result.timeLine_size);
		if isempty(ind_t)
			variable_t = [variable_t,result.timeLine_size];
			number_of_success_t = [number_of_success_t,0];
			number_of_occurrences_t = [number_of_occurrences_t,0];
			ind_t = length(variable_t);
		end
		
		if success && code==0
			number_of_success_nr(ind_nr) = number_of_success_nr(ind_nr)+1;
			number_of_success_nt(ind_nt) = number_of_success_nt(ind_nt)+1;
			number_of_success_t(ind_t) = number_of_success_t(ind_t)+1;
		end
		
		number_of_occurrences_nr(ind_nr) = number_of_occurrences_nr(ind_nr)+1;
		number_of_occurrences_nt(ind_nt) = number_of_occurrences_nt(ind_nt)+1;
		number_of_occurrences_t(ind_t) = number_of_occurrences_t(ind_t)+1;
		
		clear result;
	end

	[variable_nr, ind_nr] = sort(variable_nr);
	number_of_success_nr = number_of_success_nr(ind_nr);
	number_of_occurrences_nr = number_of_occurrences_nr(ind_nr);
	
	[variable_nt, ind_nt] = sort(variable_nt);
	number_of_success_nt = number_of_success_nr(ind_nt);
	number_of_occurrences_nt = number_of_occurrences_nt(ind_nt);
	
	[variable_t, ind_t] = sort(variable_t);
	number_of_success_t = number_of_success_t(ind_t);
	number_of_occurrences_t = number_of_occurrences_t(ind_t);

	%chi-squared test https://www.mathworks.com/help/stats/crosstab.html
	
	[~,chi,p]=crosstab(x_nr,y_nr);
	disp('NR:');
	disp(['chi-square statistic: ', num2str(chi)]);
	disp(['p-value (h0: independence): ', num2str(p)]);
	
	disp('NT:');
	[~,chi,p]=crosstab(x_nt,y_nt);
	disp(['chi-square statistic: ', num2str(chi)]);
	disp(['p-value (h0: independence): ', num2str(p)]);
	
	disp('T:');
	[~,chi,p]=crosstab(x_t,y_t);
	disp(['chi-square statistic: ', num2str(chi)]);
	disp(['p-value (h0: independence): ', num2str(p)]);

	ERR_nr = norminv(1-0.1/2) ./...
		(number_of_occurrences_nr.*sqrt(number_of_occurrences_nr)) .*...
		sqrt(number_of_success_nr.*(number_of_occurrences_nr-number_of_success_nr));
	
	ERR_nt = norminv(1-0.1/2) ./...
		(number_of_occurrences_nt.*sqrt(number_of_occurrences_nt)) .*...
		sqrt(number_of_success_nt.*(number_of_occurrences_nt-number_of_success_nt));
		
	ERR_t = norminv(1-0.1/2) ./...
		(number_of_occurrences_t.*sqrt(number_of_occurrences_t)) .*...
		sqrt(number_of_success_t.*(number_of_occurrences_t-number_of_success_t));
	
	X_nr = variable_nr;
	X_nt = variable_nt;
	X_t = variable_t;
	
	Y_nr = number_of_success_nr./number_of_occurrences_nr;
	Y_nt = number_of_success_nt./number_of_occurrences_nt;
	Y_t = number_of_success_t./number_of_occurrences_t;
end