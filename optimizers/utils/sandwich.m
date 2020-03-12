%https://staff.aub.edu.lb/~bm05/ENMG500/Set_3_revised_simplex.pdf
%This function generates two linear function (y = a1*x+b1 and y = a2*x+b2)
%The first function is a ceil for a given function defined by a lookup table
%and the secon is a floor for the fist given function.
%-x: the vector with the abciss values
%-y: the vector with the 
function [a1,b1,a2,b2] = sandwich(x,y)
	
    %objective function: minimize the difference between the sandwich functions
	%so, {a1,b1,a2,b2} = arg min {sum((a1-a2)*x+b1-b2)} = arg min{(a1-a2)sum(x)+|x|*(b1-b2)}
	%or c.'*[a1;b1;a2;b2] = c.'*z
	c = [sum(x);length(x);-sum(x);-length(x)];
	
    %constraints: guarantee the functions are respectively a ceil and a floor 
	%so, -a1*x-b1<=y and a2*x+b2<=y. Let us consider only the positive coeficients.
	%The problem is equivalent to arg min {c.'*z | A*z<=b and x>=0} where
	A = [-x, -ones(length(x),1), zeros(length(x),2);
		 zeros(length(x),2), x, ones(length(x),1)];
	b = [-y;y];
    L = eye(length(b));%the coefficients of the slack variables
    
    %initial basic solution
    z = zeros(length(b),1);%original variables
    s = b;%slack variables
    d = 0;%

	while true
        %
    end
end
