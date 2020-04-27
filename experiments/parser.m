% parser(133,'tunning_double2','_double','nr',true)
% parser(420,'tunning4','','nt',true)

function parser(n_files, folder, sufix, variable_of_interest, categorical)
	variable = [];
	number_of_success = [];
	number_of_occurrences = [];

	x=[];
	y=[];

	first_success = 0;
	second_success = 0;
	third_success = 0;
	failures = 0;

	for i=0:n_files-1
		load([folder,'\result',sufix, num2str(i),'.mat']);
		
		x=[x;result.(variable_of_interest)];
		y=[y;result.success];
		
		ind = find(variable==result.(variable_of_interest));
		if isempty(ind)
			variable = [variable,result.(variable_of_interest)];
			number_of_success = [number_of_success,0];
			number_of_occurrences = [number_of_occurrences,0];
			ind = length(variable);
		end
		
		if result.success
			number_of_success(ind) = number_of_success(ind)+1;
		end
		
		number_of_occurrences(ind) = number_of_occurrences(ind)+1;
		
		failures = failures + (1-result.success);
		first_success = first_success + result.first_success;
		second_success = second_success + result.success;
		third_success = third_success + result.success*(result.code==0);
		
		clear result;
	end

	[variable, ind] = sort(variable);
	number_of_success = number_of_success(ind);
	number_of_occurrences = number_of_occurrences(ind);

	disp(['Initial success: ', num2str(first_success)]);
	disp(['Success after V rect: ', num2str(second_success)]);
	disp(['Errors after V rect: ', num2str(third_success)]);
	disp(['Failures after V rect: ', num2str(failures)]);

	if ~categorical
		s = corr(x,y,'type','Spearman');
		p = corr(x,y,'type','Pearson');
		disp(['Pearson correlation: ', num2str(p)]);
		disp(['Spearman correlation: ', num2str(s)]);
	else
		%chi-squared test https://www.mathworks.com/help/stats/crosstab.html
		[~,chi,p]=crosstab(x,y);
		disp(['chi-square statistic: ', num2str(chi)]);
		disp(['p-value (h0: independence): ', num2str(p)]);
	end

	plot(variable,number_of_success./number_of_occurrences);
end