%Let Z be the impedance matrix such that [v.',zeros].' = Z [it.', ir.'], where v is the TX voltage vector,
%it is the transmittng current vector and ir is the receiving current vector. This function calculates a
%possible dv such that [(v+dv).',zeros].' = Z [(it+dit).', (ir+dir).'].' and:
% * abs(it+dit)<=It, 
% * abs(ir+dir)=Ir,
% * (v+dv).'*(it+dit)<=P
% * dv real
%This function is based on Newton-Raphson method, which is good since all invelved equations are quadratic
%and the solution dv=0 is near the optimum by hypothesis
function dv_values = fine_adjustment(Z, v, it, ir, It, Ir, P, tolerance, max_iterations)
    
	%the return values of voltage adjustment: [] means no found result
    dv_values = [];
    
	%number of TX and RX elements
    nt = length(it);
    nr = length(ir);
	
	iZ = eye(nt+nr)/Z;
	
	%dit = Z1*dv
	Z1 = iZ(1:nt, 1:nt);
	%dir = Z2*dv
	Z2 = iZ(nt+1:end, 1:nt);
	
	%decision variables
	x = [zeros(nt,1);P - v.'*it; It.^2 - abs(it).^2];
	
	%error vector
	F = zeros(nt+nr+1,1);
	
	%Jacobian
	J = zeros(nt+nr+1,nt+1+nt);
		
	for i = 1:max_iterations
	
		%calculating next error vector
		%tx current
		for t=1:nt
			z = [Z1(t,:), 0, zeros(1,nt)];
			F(t) = (it(t)'*z + (it(t).')*((z').'))*x + x.'*(z'*z + diag([zeros(nt,1);0;ones(nt,1)]))*x + abs(it(t))^2 - It(t)^2;
		end
		
		%rx current
		for r=1:nr
			z = [Z2(r,:), 0, zeros(1,nt)];
			F(nt+r) = (ir(r)'*z + (ir(r).')*((z').'))*x + x.'*(z'*z)*x + abs(ir(r))^2 - Ir(r)^2;
		end
		
		Z3 = [Z1, zeros(nt,nt+1);
			  zeros(nt+1,nt), zeros(nt+1)];
		
		%active power
		F(end) = v.'*it + [v.'*Z1+it.',0,zeros(1,nt)]*x + x.'*(Z3+diag([zeros(nt,1);1;zeros(nt,1)]))*x - P;
		
		if max(abs(F))<tolerance
			dv_values = x(1:nt);
			break;
		end
		
		%calculating next Jacobian
		%tx current
		for t=1:nt
			z = [Z1(t,:), 0, zeros(1,nt)];
			J(t,:) = it(t)'*z + (it(t).')*((z').') + x.'*(z'*z + (z'*z).' + 2*diag([zeros(nt,1);0;ones(nt,1)]));
		end
		
		%rx current
		for r=1:nr
			z = [Z2(t,:), 0, zeros(1,nt)];
			J(nt+r,:) = ir(r)'*z + (ir(r).')*((z').') + x.'*(z'*z + (z'*z).');
		end
		
		J(end,:) = [v.'*Z1+it.',0,zeros(1,nt)] + x.'*(Z3+Z3.'+2*diag([zeros(nt,1);1;zeros(nt,1)]));
		
		%next solution
		x = x - pinv(J)*F;
	end

end
