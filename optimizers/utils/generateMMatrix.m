%generate M purely imaginary such that M*it = -ZR*ir
function [M,it] = generateMMatrix(ZR, ir, nt)
	nr = length(ir);
	if nt<=1 || nt>nr
		disp('generateMMatrix: nt must be larger than 1 and no larger than nr');
		M = [];
		return;
	end
	
	%Create the first rows in order to generate a feasible it
	m1 = 2i*(betarnd(2,2,nt,nt)-0.5);
	it = -(eye(nt)/m1)*ZR(1:nt,:)*ir;
	
	%m2*it =-ZR(nt+1:end,:)*ir
	m2 = zeros(nr-nt,nt);
	for i = 1:nr-nt
		A = [real(it).';imag(it).'];
		b = [-(1i)*imag(ZR(nt+i,:)*ir);(1i)*real(ZR(nt+i,:)*ir)];
		m = pinv(A)*b;
		%m.'*it+ZR(nt+i,:)*ir=0
		m2(i,:) = m.';
	end
	
	M = [m1;m2];
	
	if max(abs(M*it + ZR*ir))>1e-10
		disp('generateMMatrix: failure');
		M = [];
	end
end