%https://staff.aub.edu.lb/~bm05/ENMG500/Set_3_revised_simplex.pdf
%This function generates two linear function (y = a1*x+b1 and y = a2*x+b2)
%The first function is a ceil for a given function defined by a lookup table
%and the secon is a floor for the fist given function. We admit x is sorted,
%do not have repetitions and have sufficient number of samples.
%-x: the vector with the abciss values
%-y: the vector with the 
function [a1,b1,a2,b2] = sandwich(x,y)
	%spliting the points: negative-y points must receive special care
	x_neg = x(y<0); y_neg = y(y<0);
	x_pos = x(y>=0); y_pos = y(y>=0);
	
    %objective function: minimize the difference between the sandwich functions
	%so, {a1,b1,a2,b2} = arg min {sum((a1-a2)*x+b1-b2)} = arg min{(a1-a2)sum(x)+|x|*(b1-b2)}
	%or c1.'*[a1;b1;a2;b2] = c1.'*z
	c1 = [sum(x);length(x);-sum(x);-length(x)];
	
    %constraints: guarantee the functions are respectively a ceil and a floor 
	%so, a1*x+b1>=y and a2*x+b2<=y. As the simplex method requires non-negative literals,
	%the constraints are rewritten as  -a1*x_neg-b1<=-y_neg, a2*x_pos+b2<=y_pos,
	%a1*x_pos+b1>=y_pos and -a2*x_neg-b2>=-y_neg. Simplex admits only equality constraints,
	%so slack, surplus and artificial variables are required. Furthermore, all variables
	%must be non-negative, so we use z = z1-z2
	%Thus, the constraints can be rewritten as A*[z1;z2;s1;s2;a] = b where z1>=0, z2>=0,
	%s1 (slack vector) >=0, s2 (surplus vector) >=0, a (artificial vector) >=0
	A1 = [-x_neg, -ones(length(x_neg),1), zeros(length(x_neg),2);
		  zeros(length(x_pos),2), x_pos, ones(length(x_pos),1);
		  x_pos, ones(length(x_pos),1), zeros(length(x_pos),2);
		  zeros(length(x_neg),2), -x_neg, -ones(length(x_neg),1)];
	
	%the coefficients for v = [z1; z2; s1; s2; a]; 
	A = [A1, -A, [eye(length(x)), zeros(length(x)), zeros(length(x));
				  zeros(length(x)), -eye(length(x)), eye(length(x))]];
	
	b = [-y_neg; y_pos; y_pos; -y_neg];
	
	c = [c1, -c1, zeros(2*length(x)), 100*max(c1)*ones(length(x))];
	
	%v vecotor has all variables. there are
	z1.n = 4; %four variables in z1 (one for each original variable)
	z2.n = 4; %four variables in z2 (one for each original variable)
	s1.n = length(y); %one slack variable for each <= constraint
	s2.n = length(y); %one surplus variable for each >= constraint
	a.n = length(y); %one artificial variable for each >= constraint
	%so, v is a n_z1+n_z2+n_s1+n_s2+n_a positions vector
	v = zeros(z1.n+z2.n+s1.n+s2.n+a.n,1);
	%the following indices are used to map v
	%for example, a = v(a.begin:a.end)
	z1.begin = 1;		   z1.end = z1.begin + z1.n - 1;
	z2.begin = z1.end + 1; z2.end = z2.begin + z2.n - 1;
	s1.begin = z2.end + 1; s1.end = s1.begin + s1.n - 1;
	s2.begin = s1.end + 1; s2.end = s2.begin + s2.n - 1;
	a.begin  = s2.end + 1; a.end = a.begin + a.n - 1;
	
	%the basic and non-basic variables sets: initially, the slack and
	%artificial variables are the basic ones.
	vb = [s1.begin:s1.end, a.begin:a.end];
	vn = [z1.begin:z1.end, z2.begin:z2.end; s2.begin:s2.end];
		
	while true %for the penalty method (regarding artificial variables)
		%simplex main loop
		while true
			%the coefficient matrix for the basic variables
			B = A(:,vb);
			%the coefficient matrix for the non-basic variables
			N = A(:,vb);
			%the cost vector for the basic variables
			cb = c(vb);
			%the cost vector for the non-basic variables
			cn = c(vn);
			
			%inverse of B
			iB = eye(length(B)) / B;
			
			%useful variables
			cb_iB = cb * iB;
			iB_b = iB * b; %the values of the basic variables
			
			%choosing the entering variable
			enters = -1; %invalid, initially
			max_gap = 0;
			for vn_ = vn
				z = cv_iB * N(:, vn_);
				gap = z - cn(vn_);
				if gap > max_gap
					enters = vn_;
					max_gap = gap;
				end	
			end
			
			%stop condition
			if enters==-1
				%if there is no entering variable, stop
				%the minimized value is cb_iB * b.
				v(vb) = iB_b; %the basic variables' values (the rest is already 0)
				break;
			end
			
			%one more useful variable
			iB_N_enters = iB * N(:, enters);
			
			%now, who leaves the basic set?
			leaves = -1; %invalid, initially
			min_ratio = inf;
			for vb_ = vb
				denominator = iB_N_enters(vb_);
				if denominator > 0
					numerator = iB_b(vb_);
					ratio = numerator / denominator;
					if ratio < min_ratio
						leaves = vb_;
						min_ratio = ratio;
					end
				end
			end
			
			%stop condition. unbounded. something got really wrong
			if enters==-1
				error('This linear program should not be unbounded!!');
			end
			
			%switch the variables
			vb = [vb(vb~=leaving), entering];
			vn = [vn(vn~=entering), leaving];
		end
		
		%the artificial variables must always be zero
		a = v(a.begin:a.end);
		if sum(a>0) == 0
			%ok
			break;
		else
			%penalty
			c(a.begin:a.end) = 100*c(a.begin:a.end);
		end
	end
	
	%get the variables from the vectors
	z = v(z1.begin:z1.end) - v(z2.begin:z2.end);
	a1 = z(1); b1 = z(2); a2 = z(3); b2 = z(4);
	
end
