n = 3;
Z = zeros(n);
RL = zeros(n);
for i=1:n
	for j=1:n
		Z(i,j) = rand;
	end
end

RL(end,end) = rand;

iZ = inv(Z);
M = [];
for i=1:n
	M = [M;iZ(i,end)*iZ(end,:)];
end

inv(Z+RL)
inv(Z+RL)*(Z+RL)
iZ-RL(end,end)/(1+RL(end,end)*iZ(end,end))*M
