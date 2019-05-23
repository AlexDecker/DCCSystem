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
		lambda %last high eigenvalue
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
			[Zt,Rt,~,Y] = getParammeters(obj,WPTManager);

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
			[Zt,Rt,~,Y] = getParammeters(obj,WPTManager);

			%calculating the beamforming currents again
			[obj,Itbf] = calculateBeamFormingCurrents(obj,Zt,Rt,Y);
			
			obj.Vt = (Zt + Y)*Itbf;

            WPTManager = setSourceVoltages(obj,WPTManager,obj.Vt,GlobalTime);

			evaluateEstimations(obj,WPTManager,GlobalTime)

			%scheduling the next event	
            netManager = setTimer(obj,netManager,GlobalTime,obj.timeSkip);
        end
		
		function [obj,Itbf] = calculateBeamFormingCurrents(obj,Zt,Rt,Y)
			[Vectors,Values] = eig(real(Y));

			%the receiving power is maximized when the transmitting currents are proportional to the 
			%eigenvector of greatest eigenvalue of real(Y) and all power is spent
			[obj.lambda,i] = max(diag(Values));

			maxEigenvector = Vectors(:,i);
			
			%calculating c (see the paper, specially Apendix A and equation 14)
			c = sqrt(obj.Pmax/(maxEigenvector'*(Rt + real(Y))*maxEigenvector));

			Itbf = c*maxEigenvector;
		end
			
		function [Zt,Rt,Rr,Y] = getParammeters(obj,WPTManager)
			%requires omniscience, but can be understood as a stage of pre-processing.
			Z = getZ(WPTManager.ENV,0);
			Zt = Z(1:WPTManager.nt_groups,1:WPTManager.nt_groups);
			Rt = diag(diag(Zt));
			wM = imag(Z(WPTManager.nt_groups+1:end,1:WPTManager.nt_groups));
			Zr = Z(WPTManager.nt_groups+1:end,WPTManager.nt_groups+1:end);
			Y = (wM.')*inv(Zr)*wM;
			Rr = diag(diag(Zr));%receivers' resistance
		end
		
		function evaluateEstimations(obj,WPTManager,GlobalTime)
			%Getting the parammeters (omniscient agent)
			[Zt,~,Rr,Y] = getParammeters(obj,WPTManager);

			disp(['Global time (min): ', num2str(GlobalTime/60)]);
			
			%omniscient information again, this time for getting the latest current measurements
			I = WPTManager.latestCI;
			It = I(1:WPTManager.nt_groups);
			Ir = I(WPTManager.nt_groups+1:end);%receivers' current
			charging = false;
			for i=1:length(Ir)
				SOC = getSOC(WPTManager.deviceList(i).obj.bat);
				if SOC~=1
					charging = true;
				end
				disp(['Max current for device ',num2str(i),': ', num2str(abs(Ir(i))),...
					'A; SOC(0-1): ',num2str(SOC)]);
			end
			%used to avoid waiting too long simulations
			if ~charging
				obj = endSimulation(obj);
			end

			disp(['Equation 6 error: ', num2str(mean(abs(obj.Vt-(Zt+Y)*It)))]);

			%power evaluations
			delivered = obj.lambda*obj.Pmax/(obj.lambda+1);
			rxPower = Ir'*Rr*Ir;
			rxPower2 = It'*Y*It;%If Y is good, this value will approach rxPower

			totalPower = It'*diag(diag(Zt))*It + rxPower;
			disp(['Power expected to be delivered: ',num2str(delivered),'W; power delivered: ',...
				num2str(rxPower), 'W; Power spent: ',num2str(totalPower),'W; Lambda: ', num2str(obj.lambda),...
				'; Powe calc using the accurate Y: ',num2str(rxPower2),'W']);

			disp(['Equation 9 error: ', num2str(abs(delivered-rxPower)/delivered),' (0-1)']);
		end	
    end
end
