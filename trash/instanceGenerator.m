%This function generates a feasible instance for both N-Port Power Problems. So, it generates
%a sequence of system states, the corresponding constraints, the timeLine, the duration of
%each slot and a pre-calculated voltage sequence which solves the problem.

%This function requires statistics toolbox. If you don't have it, please use the less elaborated
%version, randomSourcingInstance

function [solution, chargeData, constraints, timeLine, dt] = chargegenerateFeasibleNPPPInstance(deviceData, nt, nSegments, timeLine_size)
	
	%melhorar essa parametrização!!!!
	noise_factor = 2/timeLine_size; %deve ser menor que 1
	alpha = 5;
	beta = 0.5;
	tolerance = 1e-6;
	
	%validate input data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	nr = length(deviceData);
	if nt<=0 || nr==0 || timeLine_size<=0
		error('chargegenerateFeasibleNPPPInstance: invalid input');
	end
	
	for r=1:nr
		if ~deviceData(r).check()
			error(['DeviceData error at device ', num2str(r)]);
		end
	end	
	
	%Resistance vectors ([0.01..5])
	RT = 0.01 + 4.99*betarnd(2,2,nt,1);
	RR = 0.01 + 4.99*betarnd(2,2,nr,1);
	
	%generate the succession of states%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Q = zeros(nr, timeLine_size);
	
	%batteries capacity
	chargeData.maximum = 1000*rand(nr,1);
	
	%final charge vector (random value between 0 and the maximum)
	Q(:, end) = chargeData.maximum .* betarnd(alpha,beta,nr,1);
	
	%initial charge vector (uniformelly distributed value between 0 and the maximum
	chargeData.initial = chargeData.maximum.*rand(nr,1);
	
	chargeData.minimum = chargeData.initial;
	
	%creating the intermediate charge vectors
	for time = 1:timeLine_size-1
		%interpolate
		Q(:,time) = (chargeData.initial*(timeLine_size-time) + Q(:,end)*time)/timeLine_size;
		%random noise
		noise = zeros(nr,1);
		for r = 1:nr
			%using 68-95-99.7 rule
			sigma = noise_factor*(chargeData.maximum(r)-chargeData.initial(r))/4;
			noise(r) = normrnd(0,sigma);
		end
		Q(:,time) = Q(:,time) + normrnd(0,sigma,nr,1);
		%guarantee the vectors are valid
		Q(:,time) = max(zeros(nr,1), min(chargeData.maximum, Q(:,time)));
		
		%update the minimum
		chargeData.minimum = min(chargeData.minimum, Q(:,time));
	end
	
	%decrease the minimum charge (relaxion)
	chargeData.minimum = chargeData.minimum .* (0.98*betarnd(alpha,beta,nr,1)+0.01);
	
	%cloud as a discretization helper
	cloud = CloudHash(1, nSegments, chargeData.minimum, chargeData.maximum, 1, nt);
	
	%fixing the state sequency in order to include only centers of the hyper-cubes in the cloud
	D = cloud.discretize(chargeData.initial);
	[min_Q, max_Q] = cloud.dediscretize(D);
	chargeData.initial = (min_Q + max_Q) / 2;
	for time = 1:timeLine_size
		D = cloud.discretize(Q(:,time));
		[min_Q, max_Q] = cloud.dediscretize(D);
		Q(:,time) = (min_Q + max_Q) / 2;	
	end
	
	maxId = zeros(nr,1); %maximum discharge current
	%calculating dt in order to every state to be feasible
	dt = 0;
	%for each device
	for r = 1:nr
		%this is the charging/discharging range the device may provide
		[min_Ic, max_Ic] = deviceData(r).domain_iEffectiveChargeCurrent();
		maxId(r) = max(maxId, -min_Ic - tolerance);
		
		%for the first timeslot
		
		%the min dt for max_Ic to be respected
		dt1 = (Q(r,1)-chargeData.initial(r))/max_Ic;
		%the min dt for min_Ic to be respected
		dt2 = (Q(r,1)-chargeData.initial(r))/min_Ic;
		%take the maximum between them, so every limit is satisfied
		dt = max([dt, dt1, dt2]);
		
		%for the other timeslots:
		for time = 2:timeLine_size
			dt1 = (Q(r,time)-Q(r,time-1))/max_Ic;
			dt2 = (Q(r,time)-Q(r,time-1))/min_Ic;
			dt = max([dt, dt1, dt2]);
		end
	end
	
	%increase the timeSlot duration
	dt = dt / (0.98*betarnd(alpha,beta)+0.01);
	
	%creating the solution
	solution = [];
	for time = 1:timeLine_size
		%required charging current
		if time > 1
			ic = (Q(:,time)-Q(:,time-1))/dt;
		else
			ic = (Q(:,1)-chargeData.initial)/dt;
		end
		%the current range whose resulting charging current is ic
		%so, y1 <= abs(ir) <= y2
		y1 = zeros(nr,1);
		y2 = zeros(nr,1);
		for r = 1:nr
			[x1,x2] = deviceData(r).iEffectiveChargeCurrent(ic(r));
			[y1(r),~] = deviceData(r).iConvACDC(x1);
			[~,y2(r)] = deviceData(r).iConvACDC(x2+maxId(r));
		end
		
		for s1 = 1:sample_size
		
			%getting a valid abs(ir)
			coeff = rand(nr,1);
			abs_ir = y1.*coeff + y2.*(1-coeff);
			
			%getting a valid ir
			theta = rand(nr,1);
			ir = cos(theta).*abs_ir + sin(theta).*abs_ir;
			
			%generating a valid Z matrix for which ir is feasible and, therefore,
			%Zr*ir = M*pinv(M)*Zr*ir
			
			%wM is assumed to be in [0..1], so 
			wMR = betarnd(2,2,nr);
			wMT = betarnd(2,2,nt);
			%M must be symmetrical and with main diagona equal to zero
			wMR = wMR-diag(diag(wMR)); wMR = wMR+wMR.';
			wMT = wMT-diag(diag(wMT)); wMT = wMT+wMT.';
			
			RL = zeros(nr,1);
			for r = 1:nr
                SOC = Q(r,time)/chargeData.maximum(r);
				RL(r) = deviceData(r).getRLfromSOC(SOC);
			end
			ZR = diag(RR + RL) - (1i)*wMR;
			ZT = diag(RT) - (1i)*wMT;
			
			%obtaining a valid V vector
			
			%using moore-penrose psudo-inverse formula for solving linear systems
			%w is the free complex-vector parameter. It is chosen in order to
			%guarantee V is real.
			
			%V = ZT*iT + M.'*ir
			%0 = M*it + ZR*ir <=> it = -pinv(M)*ZR*ir + (eye(nt)-pinv(M)*M)*w
			%iff ZR*ir = M*pinv(M)*ZR*ir
			
			%M has nr rows and nt cols
			if nr < nt
				%a random M will probabily lead to a solveable undeterminated system
				M = -(1i)*betarnd(2,2,nr,nt);
				while (max(max(abs(ZR*ir-M*pinv(M)*ZR*ir)))>tolerance)
					disp('Searching for a good M...');
					M = -(1i)*betarnd(2,2,nr,nt);
				end
				
				%V = (-ZT*pinv(M)*ZR*+M.')*ir + (eye(nt)-pinv(M)*M)*w
				%V = A + B*w
				A = (-ZT*pinv(M)*ZR*+M.')*ir;
				B = eye(nt)-pinv(M)*M;
				
				%imag(V) = imag(A) + imag(B)*real(w) + real(B)*imag(w) = 0
				%imag(w) = real(B)\(-imag(A) - imag(B)*real(w))
				
				%the real part of the parameter is chosen with large std since we do not
				%have information about it
				w_r = normrnd(0,100);
				w_i = imag(w) = real(B)\(-imag(A) - imag(B)*w_r);
				w = w_r + (1i)*w_i;
				
				%real() is just to delete some eventual imaginary noise
				V = real(A + B*w);
				if max(abs(imag(A + B*w)))>tolerance
					error('real V generation has large error');
				end
			elseif nr == nt
				%a random M will probabily be invertible. Therefore, one ir leads to
				%one V, which MUST be real. So, M must be carefully built in order to
				%guarantee that.
			else
				%M must be carefully buil
			end
		end
		
		V = zeros(nt,1);
		solution = [solution, struct('Q',Q(:,time),'V',V)];
	end
	constraints = [];
	timeLine = [];
end