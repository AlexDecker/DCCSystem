classdef WPTSystem
	properties
		systemName
		nt % number of transmitters
		nr % number of receivers
		
		% topological distribution of the area
		hierarchy
		
		% components
        transmitters
        receivers
        coupler
		
	end
	
	% setup methods
	methods
		% name: name of the system. Must be unique.
		% nt: number of transmitters. Must be greater or equal than 1.
		% nr: number of transmitters. Must be greater or equal than 1.
		% top_rect: simulation area. If empty the default information will be used.
		% root: the path to the Mutual Inductance references. If empty the pwd will be used.
		function obj = WPTSystem(name, nt, nr, top_rect, root)
			obj.nt = nt;
			obj.nr = nr;
			
			if obj.nt < 1 || obj.nr < 1
				error('nt and nr sum must be at least 1.');
			end
			
			% environment setup
			obj.systemName = name;
            obj.destroy(); % to guarantee there is no other system with the same name
			new_system(name);
			open_system(name);
			set_param(name,'StartTime', '0.0');
			set_param(name,'StopTime', '0.5');
			
			% solver definitions: default
			add_block('powerlib/powergui',[name, '/powergui']);
			set_param([name, '/powergui'],'SimulationMode','Discrete');
			set_param([name, '/powergui'],'SampleTime','1e-7');
			%set_param('wpt', 'MaxStep', '1');
			
			% topological hierarchy
			if isempty(top_rect) || top_rect(1) <= top_rect(3) || top_rect(2) <= top_rect(4)
				disp('Using default area information...');
				top_rect = [0, 0, 1200, 300 * max(nt, nr)];
			end
			obj.hierarchy = Hierarchy(top_rect);
			
			% spliting the area into the TX-side circuits (left), the coupling abstraction
            % (center), and the RX-side circuits (right)
            obj.hierarchy = obj.hierarchy.horizontalCut([0.25, 0.15, 0.6]);
            [obj.hierarchy.children{1}, obj] = obj.setupTXSide(obj.hierarchy.children{1});
            [obj.hierarchy.children{2}, obj] = obj.setupCouplingAbstraction(obj.hierarchy.children{2}, root);
            [obj.hierarchy.children{3}, obj] = obj.setupRXSide(obj.hierarchy.children{3});
			
			obj.connectDevices();
		end
		
		function [hierarchy, obj] = setupTXSide(obj, hierarchy)
            obj.transmitters = cell(0);
            
			% cut into equal sub areas
            max_n = max(obj.nt, obj.nr);
			distribution = ones(1, max_n) / max_n;
			hierarchy = hierarchy.verticalCut(distribution);
			
			% build each circuit
			for i = 1 : obj.nt
                % place the circuits next to the center of the area so that the wires from the coupling
                % abstraction are shorter on average
                slot = floor((max_n - obj.nt)/2) + i;
				hierarchy.children{slot} = addPadding(hierarchy.children{slot}, 0.1);
				obj.transmitters{end + 1} = Transmitter(obj.systemName, hierarchy.children{slot}.children{1}, i);
			end
		end
		
		function [hierarchy, obj] = setupCouplingAbstraction(obj, hierarchy, root)
			hierarchy = addPadding(hierarchy, 0.6);
            obj.coupler = Coupler(obj.systemName, hierarchy.children{1}, obj.nt, obj.nr, root);
		end
		
		function [hierarchy, obj] = setupRXSide(obj, hierarchy)
            obj.receivers = cell(0);
		
			% cut into equal sub areas
            max_n = max(obj.nt, obj.nr);
			distribution = ones(1, max_n) / max_n;
			hierarchy = verticalCut(hierarchy, distribution);
			
			% build each circuit
			for i = 1 : obj.nr
                % place the circuits next to the center of the area so that the wires from the coupling
                % abstraction are shorter on average
                slot = floor((max_n - obj.nt)/2) + i;
				hierarchy.children{slot} = addPadding(hierarchy.children{slot}, 0.1);
				obj.receivers{end + 1} = Receiver(obj.systemName, hierarchy.children{slot}.children{1}, i);
			end
		end
		
		function destroy(obj)
			% close inconditionally
			bdclose(obj.systemName);
		end
	end
	
	% Runtime settings
	methods
		% Set the transmitting voltages
		function obj = setVoltages(obj, voltage_vector)
			assert(length(voltage_vector) == obj.nt);
			for (i = 1:obj.nt)
				obj.transmitters{i} = obj.transmitters{i}.setVoltage(...
					voltage_vector(i)...
				);
			end
		end
		
		function obj = setFrequency(obj, frequency)
			assert(length(frequency) == 1);
			assert(frequency > 0);
			for (i = 1:obj.nt)
				obj.transmitters{i} = obj.transmitters{i}.setFrequency(...
					frequency...
				);
			end
		end
		
		% Set fixed resistances of the rlc rings
		function obj = setResistances(obj, resistance_vector)
			assert(length(resistance_vector) == obj.nt + obj.nr);
			assert(sum(resistance_vector <= 0) == 0);
			for (i = 1:obj.nt)
				obj.transmitters{i} = obj.transmitters{i}.setResistance(...
					resistance_vector(i)...
				);
			end
			for (i = 1:obj.nr)
				obj.receivers{i} = obj.receivers{i}.setResistance(...
					resistance_vector(obj.nt + i)...
				);
			end
		end
		
		function obj = setConsummerResistances(obj, resistance_vector)
			assert(length(resistance_vector) == obj.nr);
			assert(sum(resistance_vector <= 0) == 0);
			for (i = 1:obj.nr)
				obj.receivers{i} = obj.receivers{i}.setConsummerResistance(...
					resistance_vector(i)...
				);
			end
		end
		
		% States-of-charge (%)
		function obj = setSOC(obj, SOC_vector)
			assert(length(SOC_vector) == obj.nr);
			assert(sum(SOC_vector < 0 | SOC_vector > 100) == 0);
			for (i = 1:obj.nr)
				obj.receivers{i} = obj.receivers{i}.setSOC(SOC_vector(i));
			end
		end
		
		function obj = changeCouplings(obj)
			obj.coupler = obj.coupler.changeCouplings();
			% self inductances
			L = abs(diag(obj.coupler.getInductanceMatrix()));
			% operational frequency
			f = obj.transmitters{1}.getFrequency();
			% resonance capacitances (2*pi*f*L = 1./(2*pi*f*C))
			C = 1./(4*pi^2*f^2*L);
			
			for (i = 1:obj.nt)
				obj.transmitters{i} = obj.transmitters{i}.setCapacitance(C(i));
			end
			for (i = 1:obj.nr)
				obj.receivers{i} = obj.receivers{i}.setCapacitance(C(obj.nt+i));
			end
		end
		
		function settings = getSettings(obj)
			% creating the structure to store the circuit data
			settings.voltages = zeros(obj.nt, 1);
			settings.internalImpedance = zeros(obj.nt + obj.nt, 1);
		    settings.SOC = zeros(obj.nt, 1);
		    settings.M = obj.coupler.getInductanceMatrix();
			settings.f = obj.transmitters{1}.getFrequency();
		    settings.RL = zeros(obj.nt, 1);
			
			% transmitter-wise data
			for (i = 1:obj.nt)
				settings.voltages(i) = obj.transmitters{i}.getVoltage();
				settings.internalImpedance(i) = obj.transmitters{i}.getResistance();
			end
			% receiver-wise data
			for (i = 1:obj.nr)
				settings.internalImpedance(obj.nt + i) = obj.receivers{i}.getResistance();
				settings.SOC(i) = obj.receivers{i}.getSOC();
				settings.RL(i) = obj.receivers{i}.getConsummerResistance();
			end
		end
	end
	
	% connections
	methods
		function connectDevices(obj)
			for (i = 1:obj.nt)
				add_line(obj.systemName,...
					obj.transmitters{i}.positiveHandler(),...
					obj.coupler.txNegativeHandler(i),...
					'Autorouting', 'on'...
				);
				add_line(obj.systemName,...
					obj.coupler.txPositiveHandler(i),...
					obj.transmitters{i}.negativeHandler(),...
					'Autorouting', 'on'...
				);
			end
			for (i = 1:obj.nr)
				add_line(obj.systemName,...
					obj.receivers{i}.positiveHandler(),...
					obj.coupler.rxNegativeHandler(i),...
					'Autorouting', 'on'...
				);
				add_line(obj.systemName,...
					obj.coupler.rxPositiveHandler(i),...
					obj.receivers{i}.negativeHandler(),...
					'Autorouting', 'on'...
				);
			end
		end
	end
	
	% operational methods
	methods
		function [result, t] = run(obj, plotData)
			tic;
			sim(obj.systemName);
			t = toc;
			
			result.txData = cell(0);
			result.rxData = cell(0);
			
			for i = 1:obj.nr
				acqName = obj.receivers{i}.getLoadCurrentAcquisition();
				acq = eval(acqName);
				data.loadCurrent = [acq.Time, acq.Data];
				clear acqName;
				
				acqName = obj.receivers{i}.getReceivingCurrentAcquisition();
				acq = eval(acqName);
				data.current = [acq.Time, acq.Data];
				clear acqName;
				
				acqName = obj.receivers{i}.getReceivingVoltageAcquisition();
				acq = eval(acqName);
				data.voltage = [acq.Time, acq.Data];
				clear acqName;
				
				acqName = obj.receivers{i}.getBatteryAcquisition();
				acq = eval(acqName);
				data.soc = [acq.SOC____.Time, acq.SOC____.Data];
				data.batVoltage = [acq.Voltage__V_.Time, acq.Voltage__V_.Data];
				data.batCurrent = [acq.Current__A_.Time, acq.Current__A_.Data];
				clear acq;
				
				result.rxData{end+1} = data;
			end
			
			for i = 1:obj.nt
				acqName = obj.transmitters{i}.getTransmittingVoltageAcquisition();
				acq = eval(acqName);
				data.voltage = [acq.Time, acq.Data];
				clear acqName;

				acqName = obj.transmitters{i}.getTransmittingCurrentAcquisition();
				acq = eval(acqName);
				data.current = [acq.Time, acq.Data];
				clear acqName;
				
				result.txData{end+1} = data;
			end
			
			if plotData
				for i = 1:obj.nr
					figure;
					plot(result.rxData{i}.loadCurrent(:,1), result.rxData{i}.loadCurrent(:,2));
					title(['Load Currents (',num2str(i),')']);
					figure;
					plot(result.rxData{i}.current(:,1), result.rxData{i}.current(:,2));
					title(['RX Currents (',num2str(i),')']);
					figure;
					plot(result.rxData{i}.voltage(:,1), result.rxData{i}.voltage(:,2));
					title(['RX RC Voltages (',num2str(i),')']);
					figure;
					plot(result.rxData{i}.soc(:,1), result.rxData{i}.soc(:,2));
					title(['SOC (',num2str(i),')']);
					figure;
					plot(result.rxData{i}.batVoltage(:,1), result.rxData{i}.batVoltage(:,2));
					title(['Bat Voltages (',num2str(i),')']);
					figure;
					plot(result.rxData{i}.batCurrent(:,1), result.rxData{i}.batCurrent(:,2));
					title(['Bat Currents (',num2str(i),')']);
				end
				
				for i = 1:obj.nt
					figure;
					plot(result.txData{i}.voltage(:,1), result.txData{i}.voltage(:,2));
					title(['TX RC Voltages (',num2str(i),')']);
					figure;
					plot(result.txData{i}.current(:,1), result.txData{i}.current(:,2));
					title(['TX Currents (',num2str(i),')']);
				end
			end
			
		end
	end
end