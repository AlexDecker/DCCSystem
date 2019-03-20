n = 3;
A = zeros(n);
B = zeros(n);
for i=1:n
	for j=1:n
		A(i,j) = rand;
	end
end

B(ceil(n*rand),ceil(n*rand)) = rand;

ba = B*inv(A);
g = trace(ba);

inv(A+B)
inv(A+B)*(A+B)
inv(A)-1/(1+g)*inv(A)*ba
(inv(A)-1/(1+g)*inv(A)*ba)*(A+B)
