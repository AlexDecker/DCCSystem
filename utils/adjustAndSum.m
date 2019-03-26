%adjust the temporal series and then sum the data in A with the data in B.
%Each serie must be a two line matrix (the first for data and the other for time)
function R = adjustAndSum(A, B)
	%eliminate duplicates
	[Ax,index] = unique(A(2,:));
	Ay = A(1,index);
	[Bx,index] = unique(B(2,:));
	By = B(1,index);
	%get the limits for the adjusted serie
	minX = max(A(2,1),B(2,1));
	maxX = min(A(2,end),B(2,end));
	len = min(length(A(2,:)),length(B(2,:)));
	%the new time axis
	X = linspace(minX,maxX,len);
	%adjusting both series
	Ay = interp1(Ax,Ay,X);
	By = interp1(Bx,By,X);
	%the sum itself
	R = [Ay+By;X];
end
