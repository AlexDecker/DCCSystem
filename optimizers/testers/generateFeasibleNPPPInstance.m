%This function generates a feasible instance for both N-Port Power Problems. So, it generates
%a sequence of system states, the corresponding constraints, the timeLine, the duration of
%each slot and a pre-calculated voltage sequence which solves the problem.

%This function requires statistics toolbox. If you don't have it, please use the less elaborated
%version, randomSourcingInstance

%I CURRENTLY DO NOT HAVE TIME ENOUGH TO FIX ALL BUGS, SO I ATTACHED A VERIFIER AT THE END OF THE
%FUNCTION. THUS, IT WILL RETURN SUCCESS ONLY IF A CORRECT SOLUTION WAS FOUND

function [success, solution, chargeData, constraints, timeLine, dt] = chargegenerateFeasibleNPPPInstance(deviceData, nt, nSegments, timeLine_size, sample_size)
	
	%melhorar essa parametrização!!!!
	noise_factor = 2/timeLine_size; %deve ser menor que 1
	
	betadist.alpha = 5;
	betadist.beta = 0.5;
	
	betadist_sim.alpha = 2;
	betadist_sim.beta = 2;
	
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
	
	%Fixed Receiving Resistance vector ([0.01..5])
	RR = tolerance + (5-tolerance)*betarnd(betadist_sim.alpha,betadist_sim.beta,nr,1);
	
	success = true; %default value
	
	%Z matriz and discharge current for each timeslot
	timeLine = [];
	
	solution = [];
	constraints.maxCurr = zeros(nt+nr,1);
	constraints.maxPact = 0;
	
	dt = 0;%default
	
	%generate the succession of states%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Q = zeros(nr, timeLine_size);
	
	%batteries capacity
	chargeData.maximum = (1000-tolerance)*rand(nr,1)+tolerance;
	
	%final charge vector (random value between 0 and the maximum)
	Q(:, end) = (chargeData.maximum-tolerance) .* betarnd(betadist.alpha,betadist.beta,nr,1) + tolerance;
	
	%initial charge vector (uniformelly distributed value between 0 and the maximum
	chargeData.initial = (chargeData.maximum-tolerance).*rand(nr,1) + tolerance;
	
	chargeData.minimum = (chargeData.initial-2*tolerance).*rand(nr,1) + tolerance;
	
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
		%guarantee the vectors are valid (the tolerance is used because Q>minimum charge)
		Q(:,time) = max(chargeData.minimum + tolerance, min(chargeData.maximum - tolerance, Q(:,time)));
	end
	
	if max(max(Q<=chargeData.minimum*ones(1,timeLine_size)))>0
		disp('Something is wrong with Q generation');
		success = false;
		return;
	end
	
	%decrease the minimum charge (relaxion)
	chargeData.minimum = chargeData.minimum .* ((1-tolerance)*betarnd(betadist.alpha,betadist.beta,nr,1)+tolerance);
	
	if max(max(Q<=chargeData.minimum*ones(1,timeLine_size)))>0
		disp('Something is wrong with Q-minimum relaxion');
		success=false;
	end
	
	%cloud as a discretization helper
	cloud = CloudHash(1, nSegments, chargeData.minimum, chargeData.maximum, 1, nt);
	
	%fixing the state sequency in order to include only centers of the hyper-cubes in the cloud
	D = cloud.discretize(chargeData.initial);
	[min_Q, max_Q] = cloud.dediscretize(D);
	chargeData.initial = (min_Q + max_Q) / 2;
	for time = 1:timeLine_size
		D = cloud.discretize(Q(:,time));
		[min_Q, max_Q] = cloud.dediscretize(D);
		%the center of the segment
		Q(:,time) = (min_Q + max_Q) / 2;	
	end
	
	if max(max(Q<=chargeData.minimum*ones(1,timeLine_size)))>0
		disp('Something is wrong with Q rectification');
		success=false;
		return;
	end
	
	%decrease the threshold charge (relaxion)
	chargeData.threshold = Q(:,end) .* ((1-tolerance)*betarnd(betadist.alpha,betadist.beta,nr,1)+tolerance);
	chargeData.threshold = max(chargeData.minimum + tolerance, chargeData.threshold);
	
	maxId = zeros(nr,1); %maximum discharge current
	%domain limits for the receiving current
	minIR = zeros(nr,1);
	maxIR = zeros(nr,1);
	%the domain of the inverse function
	minIn = zeros(nr,1);
	maxIn = zeros(nr,1);
	%calculating dt in order to every state to be feasible
	dt = 0;
	%for each device
	for r = 1:nr
	
		[minIR(r), maxIR(r)] = deviceData(r).domain_convACDC();
		[minIn(r),maxIn(r)] = deviceData(r).domain_iConvACDC();
		
		%maximum allowed discharge current
		[min_IC, ~] = deviceData(r).domain_effectiveChargeCurrent();
		maxId(r) = max(maxId(r), -min_IC - tolerance);
	
		%this is the charging/discharging range the device may provide
		[min_Ic, max_Ic] = deviceData(r).domain_iEffectiveChargeCurrent();
		
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
	dt = dt / ((1-tolerance)*betarnd(betadist.alpha,betadist.beta)+tolerance);
	
	%creating the solution
	for time = 1:timeLine_size
		disp(' ');
		disp(['Starting slot #',num2str(time)]);
	
		%required charging current
		if time > 1
			ic = (Q(:,time)-Q(:,time-1))/dt;
		else
			ic = (Q(:,1)-chargeData.initial)/dt;
		end
		%the receiving-current amplitude range whose resulting charging
		%current may be ic. so, y1 <= abs(ir) <= y2
		y1 = zeros(nr,1);
		y2 = zeros(nr,1);
		%the difference between input current and discharge current whose
		%resulting charging current is ic
		x1 = zeros(nr,1);
		x2 = zeros(nr,1);
		for r = 1:nr
			%the input range to generate ic(r)
			[x1(r),x2(r)] = deviceData(r).iEffectiveChargeCurrent(ic(r));
			
			%so, x1(r)<=deviceData(r).convACDC(abs(ir(r)))-Id(r)<=x2(r)
			
			%the minimum amplitude to generate ic(r)
			[y1(r),~] = deviceData(r).iConvACDC(max(x1(r),minIn(r)));
			
			%the maximum amplitude to generate ic(r)
			[~,y2(r)] = deviceData(r).iConvACDC(min(x2(r)+maxId(r),maxIn(r)));
		end
		
		%for safety reasons (eventual precision issues)
		y1 = max(minIR+tolerance, y1);
		y2 = min(maxIR-tolerance, y2);
		if max(y1>y2)==1
			disp('Failure!');
			success = false;
			return;
		end
		
		%load resistance for the chosen Q0
		RL = zeros(nr,1);
		if time == 1
			for r = 1:nr
				SOC = chargeData.initial(r)/chargeData.maximum(r);
				RL(r) = deviceData(r).getRLfromSOC(SOC);
			end
		else
			for r = 1:nr
				SOC = Q(r,time-1)/chargeData.maximum(r);
				RL(r) = deviceData(r).getRLfromSOC(SOC);
			end
		end
		
		chosen_V = zeros(nt,1);
		chosen_P = inf;
		chosen_Z = zeros(nt+nr);
		chosen_I = zeros(nt+nr,1);
		
		for s1 = 1:sample_size
			
			%getting a valid abs(ir)
			coeff = rand(nr,1);
			abs_ir = y1.*coeff + y2.*(1-coeff);
			
			%getting a valid ir
			theta = 2*pi*rand(nr,1);
			ir = cos(theta).*abs_ir + (1i)*sin(theta).*abs_ir;
			
			%generating a valid Z matrix for which ir is feasible and, therefore,
			%Zr*ir = M*pinv(M)*Zr*ir
			
			%wM is assumed to be in [0..1], so 
			wMR = betarnd(betadist_sim.alpha,betadist_sim.beta,nr);
			%M must be symmetrical and with main diagona equal to zero
			wMR = wMR-diag(diag(wMR)); wMR = wMR+wMR.';
			
			ZR = diag(RR + RL) - (1i)*wMR;
			
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
				M = -(1i)*betarnd(betadist_sim.alpha,betadist_sim.beta,nr,nt);
				while (max(max(abs(ZR*ir-M*pinv(M)*ZR*ir)))>tolerance)
					disp('Searching for a good M...');
					M = -(1i)*betarnd(betadist_sim.alpha,betadist_sim.beta,nr,nt);
				end
				
				it = -pinv(M)*ZR*ir + (eye(nt)-pinv(M)*M)*(normrnd(0,10,nt,1)+(1i)*normrnd(0,10,nt,1));
			else
				%One ir leads to one V, which MUST be real. So, M must be carefully built
				%in order to guarantee that.
				[M,it] = generateMMatrix(ZR, ir, nt);
			end
			
			%M is empty when the generator fails
			if ~isempty(M)
				rt = -ones(nt,1);
				iwMT = -(1i)*betarnd(betadist_sim.alpha,betadist_sim.beta,nt);
				iwMT = (iwMT + iwMT.') / 2;
				
				while true
				
					%calculate RT in order to imag(v)=0
					v = (1i)*iwMT*real(it)-imag(M.'*ir);
					rt = v./imag(it);
					
					%rt cannot be zero
					if max(rt<=0)==0
						break;
					end
					
					%fix rt_k
					[~,k] = min(rt);
					
					dm = (1i)*v(k)/real(it(k)) - (1i)*sign(imag(it(k)))*sign(real(it(k)))*gamrnd(1,2);
					
					iwMT(k,k) = iwMT(k,k) + dm;
				end
				
				ZT = diag(rt) + iwMT;
				
				%real() is just to delete some eventual imaginary noise
				V = real(ZT*it + M.'*ir);
				
				Z = [ZT, M.'; M, ZR];
				I_test = Z\[V;zeros(nr,1)];
				IT_test = I_test(1:nt);
				IR_test = I_test(nt+1:end);
				P = V.'*real(IT_test);
				
				if max(abs(abs(IR_test)-abs_ir))>tolerance
					disp('Invalid V');
				elseif max(max(abs(real(Z-diag(diag(Z))))))>tolerance
					disp('Invalid Z');
				elseif sum(diag(Z)<0)>0
					disp('Invalid diag(Z)');
				elseif max(minIR>abs(IR_test))==1 || max(maxIR<abs(IR_test))==1
					disp('IR out of domain');
				else
					%The solution is ok
					if P < chosen_P
						chosen_V = V;
						%the RL matrix is dynamic and then is desconsidered here
						chosen_Z = Z - diag([zeros(nt,1);RL]);
						chosen_I = I_test;
						chosen_P = P;
					end
				end
			end
			if mod(s1,round(sample_size/100))==0
				fprintf('|');
			end
		end
		
		solution = [solution, struct('Q',Q(:,time),'V',chosen_V)];
		slot.Z = chosen_Z;
		constraints.maxCurr = max(constraints.maxCurr, abs(chosen_I));
		constraints.maxPact = max(constraints.maxPact, chosen_P);
		slot.Id = zeros(nr,1);
		
		for r=1:nr
			%it is guaranteed to be inside the domain
			in = deviceData(r).convACDC(abs(chosen_I(nt+r)));
			%chose Id(r) so that x1<=in-Id(r)<=x2
			if in < x1
				disp('failure!!!');
				success = false;
				return;
			else
				%x1-in<=-Id(r)<=x2-in
				coeff = rand;
				slot.Id(r) = ((in-x1(r))*coeff + max(in-x2(r),0)*(1-coeff));
				if slot.Id(r) > maxId(r) || slot.Id(r) < 0
					disp('failure!!!');
					success = false;
					return;
				end
			end
			%verify if the the resulting effective charge is correct
			if abs(ic(r)-deviceData(r).effectiveChargeCurrent(in-slot.Id(r)))>tolerance
				disp('failure!!!!');
				success = false;
				return;
			end
		end
		
		timeLine = [timeLine; slot];
	end
	
	%increase the maximum energy spend
	constraints.maxCurr = constraints.maxCurr ./ ((1-tolerance)*betarnd(betadist.alpha,betadist.beta,nt+nr,1)+tolerance);
	constraints.maxPact = constraints.maxPact / ((1-tolerance)*betarnd(betadist.alpha,betadist.beta)+tolerance);
	
	%constraints.maxCurr = max(minIR, constraints.maxCurr);
	constraints.maxCurr(nt+1:end) = min(maxIR - tolerance, constraints.maxCurr(nt+1:end));
	
	%THIS WAS ADDED IN ORDER TO PREVENT WRONG ANSWERS CAUSED BY REMAINING BUGS
	disp(' ');
	ffModel = FeasibleFuture();  
	inst = NPortSourcingProblem(timeLine,dt,chargeData,deviceData,constraints,ffModel);
	[code,Q_list] = inst.verify([solution.V]);
	if(code~=0)
		%disp(Q_list);
		%disp([chargeData.initial,[solution.Q]]);
		disp(['ERROR CODE: ',num2str(code)]);
		success=false;
	elseif max(max(abs(Q_list-[chargeData.initial,[solution.Q]])>1e-6))>0
		disp('WARNING: THE SOLUTION IS VALID BUT THE STATE SEQUENCE IS NOT.');
	end
	disp('finished.')
end