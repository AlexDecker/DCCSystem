%This application implements the system MultiSpot proposed by Shi at al.
%("Wireless Power Hotspot that Charges All of Your Devices")
classdef multiSpot_tx < powerTXApplication
    properties
		%for the regular operation
        timeSkip %step between the time events (seconds)
		Pmax %max active power to spend
		Zt %transmitters' impedance submatrixi
		Y %receiver-receiver couplings and impedance

		%for evaluation
		Itbf %beamforming currents (transmitter)
		Vt %last voltage applied
		It %last measured transmitting current
		lambda %max eigenvalue of real(Y) matrix

		%for Y estimation (extended form)
		VtWindow
		ItWindow
		Y_safe
		rY %most accurate value of Y
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

			%This can be initialized using any matrix, according to the paper
			obj.Y = rand(WPTManager.nt_groups);
			obj.Y_safe = obj.Y;
			obj.rY = obj.Y;

			%The first iteraction (expected bad results)
			obj = calculateBeamFormingCurrents(obj);
			
			obj.Vt = (obj.Zt + obj.rY)*obj.Itbf;
			
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
			
			%For the extended way of estimating Y
			obj = calculateYSafely(obj);

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
			
			obj.Vt = (obj.Zt + obj.rY)*obj.Itbf;
            WPTManager = setSourceVoltages(obj,WPTManager,obj.Vt,GlobalTime);
			%scheduling the next event	
            netManager = setTimer(obj,netManager,GlobalTime,obj.timeSkip);
        end
		
		function obj = calculateBeamFormingCurrents(obj)
			[Vectors,Values] = eig(real(obj.rY));
			%[Vectors,Values] = eig(real(obj.Y_safe));
			%[Vectors,Values] = eig(real(obj.Y));

			%the receiving power is maximized when the transmitting currents are proportional to the 
			%eigenvector of greatest eigenvalue of real(Y) and all power is spent
			[obj.lambda,i] = max(abs(diag(Values)));%abs just in case the error for Y is still high

			maxEigenvector = Vectors(:,i);
			%as we consider all the internal impedance of each transmitter being the resistance,
			%the transmitting resistance vector is numerically equal to the main diagonal of the Zt matrix
			Rt = diag(diag(obj.Zt));
			
			%calculating c (see the paper, specially Apendix A and equation 14)
			%c = sqrt(obj.Pmax/(maxEigenvector'*(Rt + real(obj.Y))*maxEigenvector));
			%c = sqrt(obj.Pmax/(maxEigenvector'*(Rt + real(obj.Y_safe))*maxEigenvector));
			c = sqrt(obj.Pmax/(maxEigenvector'*(Rt + real(obj.rY))*maxEigenvector));

			obj.Itbf = c*maxEigenvector;
		end
			
		function Zt = getZt(obj,WPTManager)
			%requires omniscience, but can be understood as a stage of pre-processing.
			Z = getZ(WPTManager.ENV,0);
			Zt = Z(1:WPTManager.nt_groups,1:WPTManager.nt_groups);
		end
		
		function obj = calculateYSafely(obj)
			obj.VtWindow = [obj.VtWindow, obj.Vt];
			obj.ItWindow = [obj.ItWindow, obj.It];
			s = size(obj.VtWindow);
			if(s(1)<s(2))%more than nt columns
				obj.VtWindow = obj.VtWindow(:,2:end);
				obj.ItWindow = obj.ItWindow(:,2:end);
				obj.Y_safe = (obj.VtWindow/obj.ItWindow) - obj.Zt;
			else
				obj.Y_safe = rand(s(1));%not ready for real estimations
			end
		end

		function obj = evaluateEstimations(obj,WPTManager,GlobalTime)
			disp(['Global time (min): ', num2str(GlobalTime/60)]);
			%this function uses omniscient information in order to evaluate Y estimation accuracy
			Z = getCompleteLastZMatrix(WPTManager);
			%omniscient information again, this time for getting the latest current measurements
			I = WPTManager.latestCI;
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
			
			%first verification: is this I valid toward obj.It?
			disp(['Error for It measurements: ',...
				num2str(sum(abs(I(1:WPTManager.nt_groups)-obj.It))/sum(abs(obj.It))),' (0-1)']);
			
			%receivers' impedance sub-matrix
			Zr = Z(WPTManager.nt_groups+1:end,WPTManager.nt_groups+1:end);
			Rr = diag(diag(Zr));%receivers' resistance

			imgZ = imag(Z);
			wM = -imgZ(WPTManager.nt_groups+1:end,1:WPTManager.nt_groups);
			iZr = eye(length(Zr))/Zr;%inversion of Zr
			
			%testing equation 5:
			disp(['Error for eq 5: ',num2str(sum(abs((1i)*iZr*wM*obj.It-Ir))/sum(abs(Ir))), ' (0-1)']);
		
			disp(['Equation 6 error: ', num2str(mean(abs(obj.Vt-(obj.Zt+(wM')*iZr*(wM))*obj.It)))]);
			
			%the accurate value of Y
			obj.rY = (wM.')*iZr*wM;

			%power evaluations
			delivered = obj.lambda*obj.Pmax/(obj.lambda+1);
			rxPower = Ir'*Rr*Ir;
			rxPower2 = obj.It'*real(obj.Y)*obj.It;%If Y is good, this value will approach rxPower
			rxPower3 = obj.It'*real(obj.Y_safe)*obj.It;
			rxPower4 = obj.It'*real(obj.rY)*obj.It;

			totalPower = obj.It'*diag(diag(obj.Zt))*obj.It + rxPower;
			disp(['Power expected to be delivered: ',num2str(delivered),'W; power delivered: ',...
				num2str(rxPower), 'W; Power spent: ',num2str(totalPower),'W; Lambda: ', num2str(obj.lambda),...
				'; Power estimation using Y: ',num2str(rxPower2),'W; Using Y (extended calculation): ',...
				num2str(rxPower3),'W; Using the accurate Y: ',num2str(rxPower4),'W']);

			disp(['Equation 9 error: ', num2str(abs(delivered-rxPower)/delivered),' (0-1)']);


			%error between the estimated value and the real, normalized
			disp(['Error for the estimation of Y: ',...
				num2str(mean(mean(abs(obj.rY-obj.Y)))/mean(mean(abs(obj.rY)))),' (0-1)']);
			disp(['Error for the extended calculation of Y: ',...
				num2str(mean(mean(abs(obj.rY-obj.Y_safe)))/mean(mean(abs(obj.rY)))),' (0-1)']);

		end	
    end
end
