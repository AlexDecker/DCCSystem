%verify if a given measurement is dependent on a given variable_of_interest for
%a given algorithm (from 1 to 5)

function parser(n_files, folder, prefix, algorithm, measurement,...
	variable_of_interest, boxplots, marker)
	
	variable = [];
	number_of_success = [];
	number_of_occurrences = [];

	x=[];
	y=[];

	%for i=0:n_files-1
	for i=1:n_files-1
		load([folder,'\',prefix, num2str(i),'.mat']);
		
		x=[x;result.(variable_of_interest)];
		y=[y;result.ret(algorithm).(measurement)];
		
		ind = find(variable==result.(variable_of_interest));
		if isempty(ind)
			variable = [variable,result.(variable_of_interest)];
			number_of_success = [number_of_success,0];
			number_of_occurrences = [number_of_occurrences,0];
			ind = length(variable);
		end
		
		if result.ret(algorithm).success && result.ret(algorithm).code==0
			number_of_success(ind) = number_of_success(ind)+1;
		end
		
		number_of_occurrences(ind) = number_of_occurrences(ind)+1;
		
		clear result;
	end

	[variable, ind] = sort(variable);
	number_of_success = number_of_success(ind);
	number_of_occurrences = number_of_occurrences(ind);

	if ~strcmp('success',measurement)
		s = corr(x,y,'type','Spearman');
		p = corr(x,y,'type','Pearson');
		disp(['Pearson correlation: ', num2str(p)]);
		disp(['Spearman correlation: ', num2str(s)]);
		h = kruskalwallis(y,x,'off');
		disp(['Kruskal-Wallis test: p-value (H0: independency): ',num2str(h)]);
		if boxplots
			boxplot(y,x);
		else
			y_variable = variable;
			e_variable = variable;
			for k=1:length(variable)
				y_variable(k) = mean(y(x==variable(k)));
				t = tinv(0.9,sum(x==variable(k))-1);
				s = std(y(x==variable(k)));
				sqrt_n = sqrt(sum(x==variable(k)));
				e_variable(k) = t*s/sqrt_n;
			end
			errorbar(variable,y_variable,e_variable,marker,'linewidth',3);
		end
		ylim([0,inf])		
	else
		%chi-squared test https://www.mathworks.com/help/stats/crosstab.html
		[~,chi,p]=crosstab(x,y);
		disp(['chi-square statistic: ', num2str(chi)]);
		disp(['p-value (h0: independence): ', num2str(p)]);
		err = norminv(1-0.1/2) ./...
			(number_of_occurrences.*sqrt(number_of_occurrences)) .*...
			sqrt(number_of_success.*(number_of_occurrences-number_of_success));
		errorbar(variable,number_of_success./number_of_occurrences,...
			err,marker, 'linewidth', 3);
		ylim([0,1])
	end
	xlim([min(variable),max(variable)])
	set(gca,'fontsize', 12);
	set(gca,'XTick',min(variable):max(variable));
end