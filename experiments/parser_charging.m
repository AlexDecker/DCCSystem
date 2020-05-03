%verify if the normalized output of the charging algorithms is dependent
%on a given variable_of_interest for a given algorithm (from 1 to 5)

function parser_charging(n_files, folder, prefix, algorithm,...
	variable_of_interest, marker)

	x=[];
	y=[];

	for i=1:n_files-1
		load([folder,'\',prefix, num2str(i),'.mat']);
		
		if result.ret(algorithm).success && result.ret(algorithm).code==0
		%if min([result.ret.success])==1 && max(abs([result.ret.code]))==0
			
			len = length(result.timeLine);
			
			x=[x;result.(variable_of_interest)];
			
			[~,sol] = size(result.ret(algorithm).solution.Q);
			
			y=[y;sol/len];
		end
		
		clear result;
	end
	x_unique = unique(x);
	x_unique = sort(x_unique);
	
	p = kruskalwallis(y,x,'off');
	disp(['Kruskal-Wallis test: p-value (H0: independency): ',num2str(p)]);
	
	y_mean = x_unique;
	y_err = x_unique;
	for k=1:length(x_unique)
		y_mean(k) = mean(y(x==x_unique(k)));
		t = tinv(0.9,sum(x==x_unique(k))-1);
		s = std(y(x==x_unique(k)));
		sqrt_n = sqrt(sum(x==x_unique(k)));
		y_err(k) = t*s/sqrt_n;
	end
	errorbar(x_unique,y_mean,y_err,marker,'linewidth',3);
end