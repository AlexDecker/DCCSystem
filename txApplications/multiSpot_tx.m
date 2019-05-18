%This application implements the system MultiSpot proposed by Shi at al.
%("Wireless Power Hotspot that Charges All of Your Devices")
classdef multiSpot_tx < powerTXApplication
    properties
        timeSkip %step between the time events (seconds)
		Pmax %max active power to spend
		Zt %transmitters' impedance submatrixi
		Y %receiver-receiver couplings and impedance
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

			%This can be initialized using any matrix, according to the paper
			obj.Y = rand(WPTManager.nt_groups);

			%The first iteraction (expected bad results)
			obj.Itbf = calculateBeamFormingCurrents(obj);
			
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
			obj.Itbf = calculateBeamFormingCurrents(obj);
			obj.Vt = (obj.Zt + obj.Y)*obj.Itbf;
            WPTManager = setSourceVoltages(obj,WPTManager,obj.Vt,GlobalTime);
			
			%scheduling the next event	
            netManager = setTimer(obj,netManager,GlobalTime,obj.timeSkip);
        end
		
		function Itbf = calculateBeamFormingCurrents(obj)
			[Vectors,Values] = eig(real(obj.Y));
			%the receiving power is maximized when the transmitting currents are proportional to the 
			%eigenvector of greatest eigenvalue of real(Y) and all power is spent
			[~,i] = max(diag(Values));
			maxEigenvector = Vectors(:,i);
			%as we consider all the internal impedance of each transmitter being the resistance,
			%the transmitting resistance vector is numerically equal to the main diagonal of the Zt matrix
			Rt = diag(obj.Zt);
			%calculating c (see the paper, specially Apendix A and equation 14)
			c = sqrt(obj.Pmax/(maxEigenvector'*(Rt + real(obj.Y))*maxEigenvector));
			Itbf = c*maxEigenvector;
		end
			
		function Zt = getZt(obj,WPTManager)
			%requires omniscience, but can be understood as a stage of pre-processing.
			Z = getZ(WPTManager.ENV,0);
			Zt = Z(1:WPTManager.nt_groups,1:WPTManager.nt_groups);
		end
		
		function evaluateEstimations(obj,WPTManager,GlobalTime)
			%this function uses omniscient information in order to evaluate Y estimation accuracy
			Z = getCompleteLastZMatrix(WPTManager);
			%receivers' impedance sub-matrix
			Zr = Z(WPTManager.nt_groups+1:end,WPTManager.nt_groups+1:end);
			Rr = diag(diag(Zr));%receivers' resistance

			imgZ = imag(Z);
			wM = imgZ(WPTManager.nt_groups+1:end,1:WPTManager.nt_groups);

			iZr = eye(length(Zr))/Zr;%inversion of Zr
			rY = -(wM')*(iZr')*Rr*iZr*wM;

			%rY

			%real(obj.Y)

			%error between the estimated value and the real, normalized
			disp(['Error for the estimation of Y: ',...
				num2str(100*mean(mean(abs(rY-real(obj.Y))))/mean(mean(abs(rY)))),' percent']);

			disp(['Equation 6 error: ', num2str(mean(obj.Vt-(obj.Zt+(wM')*iZr*(wM))*obj.It))]);
		end	
    end
end
