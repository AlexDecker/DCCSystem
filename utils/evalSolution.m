%This function calculates the number of discrete moments (with [timeSkip] seconds of
%interval) needed to finish the charging of a WPT system defined by WPTManager (which
%is an envListManagerBAT object) without consumption, considering an optimized solu-
%tion [sol]. [sol] is a matrix in which each column corresponds to the transmitting
%voltage vector at a given moment. The function also adjusts the values of each
%vector in order to spend [Pmax] of active power given the state of the system at the
%corresponding moment. We admit a MagMIMO setup (and therefore groups with single
%coils. We also admit that the devices have a magMIMOLinearBattery.
function [t,sol] = evalSolution(WPTManager,timeSkip,Pmax,sol)
	t = inf;%default
	%extracting some parameters in order to spare function calls
	nt = WPTManager.nt_groups;
	s = size(sol);
	tmax = s(2);%number of moments=number of columns
	if(nt~=s(1))%simple compatibility test
		disp('nt must be equal to the number of rows of sol');
		return;
	end
	n = length(WPTManager.ENV.R_group);
	nr = n-nt;

	%coefficients for the nr receivers (load resistance)
	coeff = [];

	%more data
	Q = [];%charge
	Qmax = [];%maximum charge
	I0 = [];%minimal current
	I1 = [];%maximal current
	
	%getting the data from WPTManager (and preserving the original structs)
	for i = 1:nr
		%some more data, such as limit currents and initial charge
		Q = [Q; WPTManager.deviceList(i).obj.bat.Q];
		Qmax = [Qmax; WPTManager.deviceList(i).obj.bat.Qmax];
		I0 = [I0; WPTManager.deviceList(i).obj.bat.constantCurrent_min];
		I1 = [I1; WPTManager.deviceList(i).obj.bat.constantCurrent_max];

		%linear regression
		data = WPTManager.deviceList(i).obj.bat.rlLookupTable.table;
		X = [ones(length(data(:,1)),1),Qmax(end)*data(:,1)];%the original data is in
		%terms of SOC
		Y = data(:,2);
		B = ((X'*X)\X')*Y;
		coeff = [coeff; B'];%RL_i(q=q')=coeff(i,1)+coeff(i,2)*q'
	end

	t=0;%starting moment
	RL = coeff(:,1)+coeff(:,2).*Q;%initial 
	RL = (RL>0).*RL;%RL always more than 0
	while(t<tmax)
		%getting the n-port impedance matrix
		Z = getZ(WPTManager.ENV,t*timeSkip)+diag([zeros(nt,1);RL]);
		
		V = [sol(:,t+1);zeros(nr,1)];%composing V vector 
		I = Z\V;%obtaining the current

		P = real(I'*V);
		k = sqrt(Pmax/P);%correction factor
		I = k*I;%fixing the current
		sol(:,t+1) = k*sol(:,t+1);%fixing the solution

		Ir = I(nt+1:end);	
		Ir = (Ir>I0).*Ir;%cutting by minimal current
		Ir = Ir + (Ir>I1).*(I1-Ir);%cutting by maximal current

		Q1 = Q+timeSkip*Ir;%new charge
		
		if(sum(Q1<Qmax)==0)%finished. Q can be more then Qmax for simplification
			break
		end
		%updating RL
		RL = coeff(:,2).*(Q1-Q);

		Q = Q1;	
		t = t+timeSkip;
	end
end
