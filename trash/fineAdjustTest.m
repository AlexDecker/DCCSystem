clear all;

nt = randi(10);
nr = randi(10);
timeLine = randomTimeLine(nt,nr,1,zeros(nr,1),rand,rand);
Z = timeLine.Z;
V = 15*rand(nt,1);
I = Z\[V;zeros(nr,1)];
IT = I(1:nt);
IR = I(nt+1:end);

iZ = eye(nt+nr)/Z;
ZT = iZ(1:nt, 1:nt);
ZR = iZ(nt+1:end, 1:nt);

dV = 0.01*max(V)*(rand(nt,1)-0.5);
x = [dV; 0; rand(nt,1)];

dIT = ZT*dV;
dIR = ZR*dV;

for t=1:nt
	z = [ZT(t,:), 0, zeros(1,nt)];
	(IT(t)'*z + (IT(t).')*((z').'))*x + x.'*z'*z*x + abs(IT(t))^2 - abs(IT(t)+dIT(t))^2
end

disp('P:');
Z3 = [ZT, zeros(nt,nt+1);
	zeros(nt+1,nt), zeros(nt+1)];

V.'*IT + [V.'*ZT+IT.',0,zeros(1,nt)]*x + x.'*Z3*x - (V+dV).'*(IT+dIT)
