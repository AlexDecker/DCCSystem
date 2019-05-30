%This application implements the system MultiSpot proposed by Shi at al.
%("Wireless Power Hotspot that Charges All of Your Devices") but using a
%totally accurated Y matrix instead of estimations based on current
%measurements
classdef omniMultiSpot_tx < powerTXApplication
    properties
		%for the regular operation
        timeSkip %step between the time events (seconds)
		Pmax %max active power to spend
		Vt %last transmitting voltages
    end
    methods
        function obj = omniMultiSpot_tx(timeSkip,Pmax)
            obj@powerTXApplication();%building superclass structure
            obj.APPLICATION_LOG.DATA = zeros(2,0);
            obj.timeSkip = timeSkip;
			obj.Pmax = Pmax;
        end

        function [obj,netManager,WPTManager] = init(obj,netManager,WPTManager)
			%Getting the parammeters (omniscient agent)
			[Zt,~,Rt,~,~,Y] = getParammeters(obj,WPTManager);

			%The first iteraction (expected bad results)
			[obj,Itbf] = calculateBeamFormingCurrents(obj,Zt,Rt,Y);
			
			obj.Vt = (Zt + Y)*Itbf;
			
			%apply the corresponding voltages
			WPTManager = setSourceVoltages(obj,WPTManager,obj.Vt,0);

			%start the timer
            netManager = setTimer(obj,netManager,0,obj.timeSkip);
        end

        function [obj,netManager,WPTManager] = handleMessage(obj,data,GlobalTime,netManager,WPTManager)
			%Not needed
		end

        function [obj,netManager,WPTManager] = handleTimer(obj,GlobalTime,netManager,WPTManager)
			%Getting the parammeters (omniscient agent)
			[Zt,Zr,Rt,Rr,wM,Y] = getParammeters(obj,WPTManager);

			%prints omniscient information about the charging procedure
			obj = evaluateEstimations(obj,WPTManager,GlobalTime,Zt,Zr,Rt,Rr,wM,Y);

			%calculating the beamforming currents again
			[obj,Itbf] = calculateBeamFormingCurrents(obj,Zt,Rt,Y);
			
			obj.Vt = (Zt + Y)*Itbf;

            WPTManager = setSourceVoltages(obj,WPTManager,obj.Vt,GlobalTime);

			%scheduling the next event	
            netManager = setTimer(obj,netManager,GlobalTime,obj.timeSkip);
        end
		
		function [obj,Itbf] = calculateBeamFormingCurrents(obj,Zt,Rt,Y)
			[Vectors,Values] = eig(real(Y));
			%the receiving power is maximized when the transmitting currents are proportional to the 
			%eigenvector of greatest eigenvalue of real(Y) and all power is spent
			[~,i] = max(diag(Values));

			maxEigenvector = Vectors(:,i);
			
			%calculating c (see the paper, specially Apendix A and equation 14)
			c = sqrt(obj.Pmax/(maxEigenvector'*(Rt + real(Y))*maxEigenvector));

			Itbf = c*maxEigenvector;
		end
			
		function [Zt,Zr,Rt,Rr,wM,Y] = getParammeters(obj,WPTManager)
			%requires omniscience, but can be understood as a stage of pre-processing.
			Z = getCompleteLastZMatrix(WPTManager);
			Zt = Z(1:WPTManager.nt_groups,1:WPTManager.nt_groups);
			Rt = diag(diag(Zt));
			wM = imag(Z(WPTManager.nt_groups+1:end,1:WPTManager.nt_groups));
			Zr = Z(WPTManager.nt_groups+1:end,WPTManager.nt_groups+1:end);
			Y = (wM.')*inv(Zr)*wM;
			Rr = diag(diag(Zr));%receivers' resistance
		end
		
		function obj = evaluateEstimations(obj,WPTManager,GlobalTime,Zt,Zr,Rt,Rr,wM,Y)
			disp('ESTIMATIONS:---------------------------------X');
			disp(['Global time (h): ', num2str(GlobalTime/3600)]);
			
			
			disp('State of the devices:');	
			%omniscient information again, this time for getting the latest current measurements
			I = WPTManager.latestCI;
			It = I(1:WPTManager.nt_groups);
			Ir = I(WPTManager.nt_groups+1:end);%receivers' current
			charging = false;
			for i=1:length(Ir)
				SOC = getSOC(WPTManager.deviceList(i).obj.bat);
				if SOC~=1
					charging = true;
					disp(['Max current for device ',num2str(i),': ', num2str(abs(Ir(i))),...
						'A; SOC(0-1): ',num2str(SOC)]);
				else
					disp(['Device ',num2str(i),' ended charging']);
				end
			end
			%used to avoid waiting too long simulations
			if ~charging
				disp('All devices ended charging');
				obj = endSimulation(obj);
			end


			%the assumptions of the model: the currents and voltages are related according to equation 6
			disp(['Equation 6 error: ', num2str(mean(abs(obj.Vt-(Zt+Y)*It)))]);

			%The accurate delivered power
			rxPowerRr = Ir'*Rr*Ir;

			%The delivered power using Y
			rxPowerY = It'*real(Y)*It;%If Y is good, this value will approach rxPower

			%rx power + tx power
			totalPower = It'*Rt*It + rxPowerRr;

			%x vector (see Apendix A)
			x = sqrt(Rt)*It;
			A = sqrt(inv(Rt))*real(Y)*sqrt(inv(Rt));
			%this value must be equal to rxPower
			rxPowerA = x'*A*x;

			[Q,Lambda] = eig(A);
			%Q must be orthogonal
			disp(['Q orthogonality error: ', num2str(mean(mean(abs(eye(length(Q))-Q'*Q))))]);

			%y corresponds to x' from the paper (see Apendix A). Y must have only the first position other 
			%than 0
			y = Q'*x;

			%This value must be equal to rxPower
			rxPowery = y'*Lambda*y;
			
			%getting the maximal eigenvalue (the max eigenvalue of real(Y) cannot be used for this
			%estimation, because some constants are omitted for simplicity by the authors)
			lambda = max(max(Lambda));
			
			%The power expected to be delivered according to theorem 4.1
			delivered = lambda*obj.Pmax/(lambda+1);
			
			disp(['rx Power values in W (they all must be equal): ',num2str(rxPowerRr),...
				';',num2str(rxPowerY),';',num2str(rxPowerA),';',num2str(rxPowery),';',...
				num2str(delivered)]);
			disp(['Total power spent (W): ', num2str(totalPower)]);
			disp('X--------------------------------------------X');
		end	
    end
end
