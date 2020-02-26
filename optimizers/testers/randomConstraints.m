%this function generates the constraints based on provided data in order to
%guarantee there is at least one solution. It admits homogeneous networks.
% * deviceData: DeviceData object, which encapsulates the lookup tables involved 
%in the simulation process
% * timeLine: List with the impedance matrix and discharge vector at each moment
% * dt: Duration of one time slot (s)
% * sampleSize: This function chooses an efficient voltage vector for each timeslot.
%the search for the vector is via random sampling. So, the larger this value, the
%harder the problem instance.
% * maxV: maximum allowed voltage (absolute value) in volts.
function [constraints,chargeData,sol,success] = randomConstraints(deviceData,...
	timeLine,dt,sampleSize,maxV)
	
	success = true;
	
	%getting the dimensions of the problem
	nSlots = length(timeLine);
	nr = length(timeLine(1).Id);
	nt = length(timeLine(1).Z)-nr;

    %maximum receiving current (absolute value)
    [~,maxIr] = deviceData.domain_convACDC();
	
	chargeData.initial = rand(nr,1)*maxIr*dt*nSlots;
	chargeData.minimum = chargeData.initial - 1e-6; %to guarantee initial>min
	chargeData.maximum = chargeData.initial + rand(nr,1)*maxIr*dt*nSlots;
	chargeData.threshold = chargeData.initial;
	
	constraints.maxPapp = 0;
	constraints.maxPact = 0;
	constraints.maxCurr = zeros(nt+nr,1);
	
	%solution matrix
	sol = zeros(nt,0);
	
	%charge vector
	q = chargeData.initial;
    qLog = q;
	
	%integrating...
	for t=1:nSlots
		%calculating the load resistance of each receiving device
		Rl = [];
		for r = 1:nr
			Rl = [Rl; deviceData.getRLfromSOC(...
				q(r)/chargeData.maximum(r))];
		end
		%greedy solution. Generate a sample of currents and choose the element
		%which maximizes the efficiency
		eff = 0;
		V = zeros(nt,1);
		for i=1:sampleSize
			v = maxV*rand(nt,1);
			%calculating the phasor current vector
			current = (timeLine(t).Z+diag([zeros(nt,1);Rl]))\[v;zeros(nr,1)];
			receiving_current = current(nt+1:end);
			%normalizing
			k = min((maxIr-1e-6)./abs(receiving_current));
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
			curr = deviceData.convACDC(abs(I(nt+r)))-...
				timeLine(t).Id(r);
			%the charge/discharge current
			chargeCurrent = [chargeCurrent;...
				deviceData.effectiveChargeCurrent(curr)];
		end

		%updating the charge vector
		q = min(q + dt*chargeCurrent, chargeData.maximum);
        qLog = [qLog,q];

		chargeData.minimum = min(q-1e-6,chargeData.minimum);

		if sum(chargeData.minimum<0)>0
			success = false;
			return;
		end
	end

	%make the problem easier
    chargeData.minimum = 0.9*chargeData.minimum;
	chargeData.threshold = 0.9*(q-chargeData.minimum)+chargeData.minimum;
    constraints.maxPapp = (1.1)*constraints.maxPapp;
    constraints.maxPact = (1.1)*constraints.maxPact;
    constraints.maxCurr = (1.1)*constraints.maxCurr;
	%maxCurr in RX is limited by maxIr
	constraints.maxCurr = min(constraints.maxCurr, [inf*ones(nt,1);maxIr*ones(nr,1)]);
end
