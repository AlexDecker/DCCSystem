%This application implements the system MultiSpot proposed by Shi at al.
%("Wireless Power Hotspot that Charges All of Your Devices")
classdef multiSpot_tx < powerTXApplication
    properties
		%for the regular operation
        timeSkip %step between the time events (seconds)
		Pmax %max active power to spend
		Zt %transmitters' impedance submatrix
		Y %receiver-receiver couplings and impedance

		%for evaluation
		Itbf %beamforming currents (transmitter)
		Vt %last voltage applied
		It %last measured transmitting current
    end
    methods
        function obj = multiSpot_tx(timeSkip,Pmax)
            obj@powerTXApplication();%building superclass structure
            obj.APPLICATION_LOG.DATA = zeros(2,0);
            obj.timeSkip = timeSkip;
			obj.Pmax = Pmax;
        end

        function [obj,netManager,WPTManager] = init(obj,netManager,WPTManager)
			%This is pre-parametrized
			obj.Zt = getZt(obj,WPTManager);

			%This can be initialized using any matrix, according to the paper.
			%But actually initializing it with a matrix so much different from the
			%actual one cause divergence.
			obj.Y = zeros(WPTManager.nt_groups);

			%The first iteraction (expected bad results)
			obj = calculateBeamFormingCurrents(obj);
			
			obj.Vt = (obj.Zt + obj.Y)*obj.Itbf;
			
			%apply the corresponding voltages
			WPTManager = setSourceVoltages(obj,WPTManager,obj.Vt,0);

			%start the timer
            netManager = setTimer(obj,netManager,0,obj.timeSkip);
        end

        function [obj,netManager,WPTManager] = handleMessage(obj,data,GlobalTime,netManager,WPTManager)
			%Not needed
		end

        function [obj,netManager,WPTManager] = handleTimer(obj,GlobalTime,netManager,WPTManager)
			%getting the current measures
            [obj.It,WPTManager] = getCurrents(obj,WPTManager,GlobalTime);
			
			%computing the error
			deltaIt = obj.Itbf-obj.It;

			%if the channel changes (deltaIt is other than the null vector)
			if (sum(deltaIt~=0)>0)
				deltaVt = (obj.Zt + obj.Y)*deltaIt;
				obj.Y = obj.Y + deltaVt*(deltaVt.')/(deltaVt.'*obj.It);
				evaluateEstimations(obj,WPTManager,GlobalTime);
			end
			
			%calculating the beamforming currents again
			obj = calculateBeamFormingCurrents(obj);
			
			obj.Vt = (obj.Zt + obj.Y)*obj.Itbf;
            WPTManager = setSourceVoltages(obj,WPTManager,obj.Vt,GlobalTime);

			%scheduling the next event	
            netManager = setTimer(obj,netManager,GlobalTime,obj.timeSkip);
        end
		
		function obj = calculateBeamFormingCurrents(obj)
			[Vectors,Values] = eig(real(obj.Y));

			%the receiving power is maximized when the transmitting currents are proportional to the 
			%eigenvector of greatest eigenvalue of real(Y) and all power is spent
			[~,i] = max(abs(diag(Values)));%abs just in case the error for Y is still high

			maxEigenvector = Vectors(:,i);
			%as we consider all the internal impedance of each transmitter being the resistance,
			%the transmitting resistance vector is numerically equal to the main diagonal of the Zt matrix
			Rt = diag(diag(obj.Zt));
			
			%calculating c (see the paper, specially Apendix A and equation 14)
			c = sqrt(obj.Pmax/(maxEigenvector'*(Rt + real(obj.Y))*maxEigenvector));

			obj.Itbf = c*maxEigenvector;
		end
			
		function Zt = getZt(obj,WPTManager)
			%requires omniscience, but can be understood as a stage of pre-processing.
			Z = getCompleteLastZMatrix(WPTManager);
			Zt = Z(1:WPTManager.nt_groups,1:WPTManager.nt_groups);
		end
		
		function obj = evaluateEstimations(obj,WPTManager,GlobalTime)
			disp('ESTIMATIONS:---------------------------------X');
			disp(['Global time (h): ', num2str(GlobalTime/60)]);

			%this function uses omniscient information in order to evaluate Y estimation accuracy
			Z = getCompleteLastZMatrix(WPTManager);

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
			
			%first verification: is this I valid toward obj.It?
			disp(['Error for It measurements: ',...
				num2str(sum(abs(It-obj.It))/sum(abs(It))),' (0-1)']);
			
			%getting some useful identities
			Zt = Z(1:WPTManager.nt_groups,1:WPTManager.nt_groups);
			Zr = Z(WPTManager.nt_groups+1:end,WPTManager.nt_groups+1:end);
			Rt = diag(diag(Zt));%transmitters' resistance
			Rr = diag(diag(Zr));%receivers' resistance
			wM = -imag(Z(WPTManager.nt_groups+1:end,1:WPTManager.nt_groups));
			iZr = eye(length(Zr))/Zr;%inversion of Zr
			Y = (wM.')*iZr*wM;%the accurate value of Y
			A = sqrt(inv(Rt))*real(Y)*sqrt(inv(Rt));

			%testing equations 5 and 6, which refers to the adopted WPT modeling and must return a very
			%small error in order to guarantee that the assumptions made by Shi et al will be ok
			disp(['Error for eq 5: ',num2str(sum(abs((1i)*iZr*wM*obj.It-Ir))/sum(abs(Ir))), ' (0-1)']);
			disp(['Error for eq 6: ', num2str(mean(abs(obj.Vt-(obj.Zt+(wM')*iZr*(wM))*obj.It)))]);
			
			%The accurate delivered power
			rxPowerRr = Ir'*Rr*Ir;

			%The delivered power using Y
			rxPowerY = It'*real(Y)*It;%If Y is good, this value will approach rxPower

			%The delivered power according to the estimated Y
			rxPowereY = It'*real(obj.Y)*It;

			%tx power + rx power
			totalPower = It'*Rt*It + rxPowerRr;

			%x vector (see Apendix A)
		    x = sqrt(Rt)*It;
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
			disp(['rx Power according to the estimated Y (W) ', num2str(rxPowereY)]);
			disp(['Total power spent (W): ', num2str(totalPower)]);
			disp('X--------------------------------------------X');
		end	
    end
end
