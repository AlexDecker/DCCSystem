%this function generates the constraints based on provided data in order to
%guarantee there is at least one solution. It admits homogeneous networks
function [constraints,chargeData,sol,success] = randomConstraints(deviceData,...
	timeLine,dt,sampleSize,maxV,maxIn)
	
	success = true;
	
	%getting the dimensions of the problem
	nSlots = length(timeLine);
	nr = length(timeLine(1).Id);
	nt = length(timeLine(1).Z)-nr;
	
	chargeData.initial = rand(nr,1)*maxIn*dt*nSlots;
	chargeData.minimum = chargeData.initial - 1e-6; %to guarantee initial>min
	chargeData.maximum = chargeData.initial+rand(nr,1)*maxIn*dt*nSlots;
	chargeData.threshold = chargeData.initial;
	
	constraints.maxPapp = 0;
	constraints.maxPact = 0;
	constraints.maxCurr = zeros(nt+nr,1);
	
	%solution matrix
	sol = zeros(nt,0);
	
	%charge vector
	q = chargeData.initial;
	
	%integrating...
	for t=1:nSlots
		%calculating the load resistance of each receiving device
		Rl = [];
		for r = 1:nr
			Rl = [Rl; deviceData.getRLfromSOC(...
				q(r,end)/chargeData.maximum(r))];
		end
		%greedy solution. Generate a sample of currents and choose the element
		%which maximizes the efficiency
		eff = 0;
		V = zeros(nt,1);
		for i=1:sampleSize
			v = maxV*rand(nt,1);
			%calculating the phasor current vector
			current = (timeLine(t).Z+diag([zeros(nt,1);Rl]))\[v;zeros(nr,1)];
			%normalizing
			k = min((maxIn-1e-6)./abs(current));
			v = k*v;
			current = k*current;
			%calculating the active power
			pact = real(current(1:nt)'*v);
			%comparing the efficiency
			if sum(abs(current))/pact>eff
				eff = sum(abs(current))/pact;
				V = v;
			end
		end
		sol = [sol, V];
		I = (timeLine(t).Z+diag([zeros(nt,1);Rl]))\[V;zeros(nr,1)];
		%updating constraints
		constraints.maxCurr = max(constraints.maxCurr, abs(I));
		constraints.maxPapp = max(constraints.maxPapp, abs(I(1:nt)'*V));
		constraints.maxPact = max(constraints.maxPact, real(I(1:nt)'*V));
		
		%converting the currents
		chargeCurrent = [];
		for r = 1:nr
			%the input current less the discharge current
			curr = deviceData.convACDC(abs(current(nt+r)))-...
				timeLine(t).Id(r);
			%the charge/discharge current
			chargeCurrent = [chargeCurrent;...
				deviceData.effectiveChargeCurrent(curr)];
		end

		%updating the charge vector
		q = min(q + dt*chargeCurrent, chargeData.maximum);

		chargeData.minimum = min(q,chargeData.minimum);
		if sum(chargeData.minimum<0)>0
			success = false;
			return;
		end
	end
	
	chargeData.threshold = q;
end