function print_vector(label, vector, spacing)
	ret = [];
	for i = 1:spacing
		ret = [ret, '...'];
	end
	
	ret = [ret, label, ': | '];
	
	for i = 1:length(vector)
		ret = [ret, num2str(vector(i)), ' | '];
	end
	
	disp(ret);
end