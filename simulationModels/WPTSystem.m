classdef WPTSystem
	properties
		name
		nt % number of transmitters
		nr % number of receivers
		
		% topological distribution of the area
		hierarchy
		
		% tx-side components
		sources
		tx_rlc_branches
		
		% rx-side components
		rx_rlc_branches
		acdc_converters
		consumers
		batteries
		
		% coupling model
		mutual_coupler % block
		coupling_helper % mutual inductance generator
		
		operating_frequency
	end
	
	% setup methods
	methods
		function obj = WPTSystem(name, nt, nr, top_rect)
		
			% local setup (get all coupling references in the pwd)
			obj.coupling_helper = CouplingHelper([], false);
			obj.nt = nt;
			obj.nr = nr;
			obj.operating_frequency = 1e+6; % 1 MHz
			
			if obj.nt < 1 || obj.nr < 1
				error('nt and nr sum must be at least 1.');
			end
			
			% environment setup
			obj.name = name;
			new_system(name);
			open_system(name);
			
			% solver definitions: default
			add_block('powerlib/powergui',[name, '/powergui']);
			%set_param([name, '/powergui'],'SimulationMode','Continuous');
			set_param([name, '/powergui'],'SimulationMode','Discrete');
			%set_param([name, '/powergui'],'SimulationMode','Phasor');
			set_param([name, '/powergui'],'SampleTime','1e-6');
			set_param('wpt', 'MaxStep', '1');
			
			% topological hierarchy
			if isempty(top_rect) || top_rect(1) <= top_rect(3) || top_rect(2) <= top_rect(4)
				disp('Using default area information...');
				top_rect = [0, 0, 1200, 300 * max(nt, nr)];
			end
			obj.hierarchy = Hierarchy(top_rect);
			
			% block containers setup
			obj.sources.count = 0;
			obj.sources.elements = cell(obj.nt,1);
			
			obj.tx_rlc_branches.count = 0;
			obj.tx_rlc_branches.components = cell(obj.nt,1);
			obj.tx_rlc_branches.voltmeters = cell(obj.nt,1);
			obj.tx_rlc_branches.acquisitions = cell(obj.nt,1);
			
			obj.rx_rlc_branches.count = 0;
			obj.rx_rlc_branches.components = cell(obj.nr,1);
			obj.rx_rlc_branches.voltmeters = cell(obj.nr,1);
			obj.rx_rlc_branches.acquisitions = cell(obj.nr,1);
			
			obj.acdc_converters.count = 0;
			obj.acdc_converters.bridges = cell(obj.nr,1);
			obj.acdc_converters.filters = cell(obj.nr,1);
			obj.acdc_converters.diodes = cell(obj.nr,1);
			
			obj.consumers.count = 0;
			obj.consumers.elements = cell(obj.nr,1);
			
			obj.batteries.count = 0;
			obj.batteries.components = cell(obj.nr,1);
			obj.batteries.acquisitions = cell(obj.nr,1); % data acquisition
			
			% dividing the area between tx, coupling, and rx panels
			obj.hierarchy = horizontalCut(obj.hierarchy, [0.14, 0.2, 0.66]);
			[obj.hierarchy.children{1}, obj] = obj.setupTXSide(obj.hierarchy.children{1});
			[obj.hierarchy.children{2}, obj] = obj.setupCouplingAbstraction(obj.hierarchy.children{2});
			[obj.hierarchy.children{3}, obj] = obj.setupRXSide(obj.hierarchy.children{3});
			
			% add all connections
			obj = obj.connectTXComponents();
			obj = obj.connectRXComponents();
		end
		
		function [hierarchy, obj] = setupTXSide(obj, hierarchy)
			
			% cut into equal sub areas
			distribution = ones(1, obj.nt) / obj.nt;
			hierarchy = verticalCut(hierarchy, distribution);
			
			% build each circuit
			for i = 1 : length(hierarchy.children)
				hierarchy.children{i} = addPadding(hierarchy.children{i}, 0.1);
				[hierarchy.children{i}.children{1}, obj] = obj.buildTXCircuit(hierarchy.children{i}.children{1});
			end
			
		end
		
		function [hierarchy, obj] = setupCouplingAbstraction(obj, hierarchy)
			hierarchy = addPadding(hierarchy, 0.6);
			obj = obj.newCoupledInductors(hierarchy.children{1});
		end
		
		function [hierarchy, obj] = setupRXSide(obj, hierarchy)
		
			% cut into equal sub areas
			distribution = ones(1, obj.nr) / obj.nr;
			hierarchy = verticalCut(hierarchy, distribution);
			
			% build each circuit
			for i = 1 : length(hierarchy.children)
				hierarchy.children{i} = addPadding(hierarchy.children{i}, 0.1);
				[hierarchy.children{i}.children{1}, obj] = obj.buildRXCircuit(hierarchy.children{i}.children{1});
			end
			
		end
		
		function [hierarchy, obj] = buildTXCircuit(obj, hierarchy)
			
			hierarchy = verticalCut(hierarchy, [0.2, 0.23, 0.02, 0.37, 0.18]);
			
			hierarchy.children{2} = horizontalCut(hierarchy.children{2}, [0.3, 0.7]);
			obj = obj.newSource(hierarchy.children{2}.children{1});
			
			hierarchy.children{4} = horizontalCut(hierarchy.children{4}, [0.2, 0.8]);
			[hierarchy.children{4}.children{2}, obj] = obj.newRC(true, hierarchy.children{4}.children{2});
			
		end
		
		function [hierarchy, obj] = buildRXCircuit(obj, hierarchy)
			% 2 vertical spacers
			hierarchy = verticalCut(hierarchy, [0.16,0.61,0.23]);
			
			% Panels: rc, acdc converter, consumer, battery and 3 spacers
			hierarchy.children{2} = horizontalCut(hierarchy.children{2}, [0.17,0.04,0.37,0.04,0.03,0.04,0.32]);
			
			hierarchy.children{2}.children{1} = verticalCut(hierarchy.children{2}.children{1}, [0.2, 0.65, 0.15]);
			[hierarchy.children{2}.children{1}.children{2}, obj] = obj.newRC(false, hierarchy.children{2}.children{1}.children{2});
			
			[hierarchy.children{2}.children{3}, obj] = obj.newACDCConverter(hierarchy.children{2}.children{3});
			
			hierarchy.children{2}.children{5} = verticalCut(hierarchy.children{2}.children{5}, [0.35,0.3,0.35]);
			obj = obj.newPoweredDevice(hierarchy.children{2}.children{5}.children{2});
			
			[hierarchy.children{2}.children{7}, obj] = obj.newBattery(0, hierarchy.children{2}.children{7});
			
		end
				
		function destroy(obj)
			% close inconditionally
			bdclose(obj.name);
		end
	end
	
	% Block creation
	methods
		function obj = newSource(obj, hierarchy)
		
			obj.sources.count = obj.sources.count + 1;
			
			element.name = [obj.name,'/src_', num2str(obj.sources.count)];
			
			add_block('powerlib/Electrical Sources/AC Voltage Source',...
				element.name,...
				'Position', hierarchy.rect ...
			);
				
			% for later connection of blocks
			element.hnd = get_param(element.name, 'PortHandles');
			
			obj.sources.elements{obj.sources.count} = element;
			
		end
		
		% used in both TX and RX (one for each circuit)
		function [hierarchy, obj] = newRC(obj, isTX, hierarchy)
			
			hierarchy = verticalCut(hierarchy, [0.16, 0.37, 0.47]);
			hierarchy.children{1} = horizontalCut(hierarchy.children{1}, [0.53, 0.47]);
			hierarchy.children{3} = horizontalCut(hierarchy.children{3}, [0.11, 0.32, 0.25, 0.32]);
			
			if isTX
				obj.tx_rlc_branches.count = obj.tx_rlc_branches.count + 1;
				count = obj.tx_rlc_branches.count;
				rlc_element.name = [obj.name,'/tx_rlc_', num2str(count)];
				voltmeter_element.name = [obj.name,'/tx_volt_', num2str(count)];
				acquisition_element.name = [obj.name,'/tx_acq_', num2str(count)];
				acquisition_element.variable = ['tx_acq_', num2str(count)];
			else
				obj.rx_rlc_branches.count = obj.rx_rlc_branches.count + 1;
				count = obj.rx_rlc_branches.count;
				rlc_element.name = [obj.name,'/rx_rlc_', num2str(count)];
				voltmeter_element.name = [obj.name,'/rx_volt_', num2str(count)];
				acquisition_element.name = [obj.name,'/rx_acq_', num2str(count)];
				acquisition_element.variable = ['rx_acq_', num2str(count)];
			end

			add_block('powerlib/Elements/Series RLC Branch',...
				rlc_element.name,...
				'Position', hierarchy.children{1}.children{1}.rect ...
			);
			
			add_block('powerlib/Measurements/Voltage Measurement',...
				voltmeter_element.name,...
				'Position', hierarchy.children{3}.children{2}.rect ...
			);
			
			add_block('simulink/Sinks/To Workspace',...
				acquisition_element.name,...
				'Position', hierarchy.children{3}.children{4}.rect ...
			);
			
			set_param(rlc_element.name, 'BranchType', 'RC');
			set_param(rlc_element.name, 'Measurements', 'Branch voltage');
			set_param(acquisition_element.name, 'VariableName', acquisition_element.variable);
			
			% for the connection of the blocks
			rlc_element.hnd = get_param(rlc_element.name,'PortHandles');
			voltmeter_element.hnd = get_param(voltmeter_element.name,'PortHandles');
			acquisition_element.hnd = get_param(acquisition_element.name,'PortHandles');
			
			add_line(obj.name,...
				rlc_element.hnd.RConn,...
				voltmeter_element.hnd.LConn(1),...
				'Autorouting', 'on'...
			);
			
			add_line(obj.name,...
				rlc_element.hnd.LConn,...
				voltmeter_element.hnd.LConn(2),...
				'Autorouting', 'on'...
			);
			
			add_line(obj.name,...
				voltmeter_element.hnd.Outport,...
				acquisition_element.hnd.Inport,...
				'Autorouting', 'on'...
			);
			
			if isTX
				obj.tx_rlc_branches.components{count} = rlc_element;
				obj.tx_rlc_branches.voltmeters{count} = voltmeter_element;
				obj.tx_rlc_branches.acquisitions{count} = acquisition_element;
			else
				obj.rx_rlc_branches.components{count} = rlc_element;
				obj.rx_rlc_branches.voltmeters{count} = voltmeter_element;
				obj.rx_rlc_branches.acquisitions{count} = acquisition_element;
			end
			
		end
		
		% https://www.mathworks.com/help/physmod/sps/powersys/ref/mutualinductance.html
		function obj = newCoupledInductors(obj, hierarchy)

			element.name = [obj.name,'/coupler'];

			add_block('powerlib/Elements/Mutual Inductance',...
				element.name, ...
				'Position', hierarchy.rect ...
			);
			
			set_param(element.name, 'TypeOfMutual', 'Generalized mutual inductance');
			set_param(element.name, 'NumberOfWindings', num2str(obj.nt + obj.nr));
			
			% for later connection of blocks
			element.hnd = get_param(element.name,'PortHandles');
			
			obj.mutual_coupler = element;
			
			% first set of internal parameters
			obj = obj.changeCouplings();
		end
		
		% Chargeable battery which supplies the powered device. 
		function [hierarchy, obj] = newBattery(obj, SOC, hierarchy)
			
			if SOC < 0 || SOC > 100
				error('Invalid state-of-charge.');
			end
			
			hierarchy = horizontalCut(hierarchy, [0.5, 0.1, 0.4]);
			
			obj.batteries.count = obj.batteries.count + 1;
			
			bat_element.name = [obj.name,'/bat_', num2str(obj.batteries.count)];
			acquisition_element.name = [obj.name,'/bat_acq_', num2str(obj.batteries.count)];
			acquisition_element.variable = ['bat_acq_', num2str(obj.batteries.count)];
			
			add_block('electricdrivelib/Extra Sources/Battery',...
				bat_element.name, ...
				'Position', hierarchy.children{1}.rect ...
			);
			
			add_block('simulink/Sinks/To Workspace',...
				acquisition_element.name,...
				'Position', hierarchy.children{3}.rect ...
			);
			
			set_param(bat_element.name, 'Orientation', 'right');
			set_param(bat_element.name, 'SOC', num2str(SOC));
			set_param(acquisition_element.name, 'VariableName', acquisition_element.variable);
			
			% LConn [1..2], Outport
			bat_element.hnd = get_param(bat_element.name,'PortHandles');
			% Inport
			acquisition_element.hnd = get_param(acquisition_element.name,'PortHandles');
			
			% internal connections
			add_line(obj.name,...
				bat_element.hnd.Outport,...
				acquisition_element.hnd.Inport,...
				'Autorouting', 'on'...
			);
			
			obj.batteries.components{obj.batteries.count} = bat_element;
			obj.batteries.acquisitions{obj.batteries.count} = acquisition_element;
		end
		
		% Provides interface between the rx rlc ring and the DC internal circuit
		function [hierarchy, obj] = newACDCConverter(obj, hierarchy)
			
			% 3 panels: one for the bridge, the second one for the filter capacitor
			% and the third one for the isolation diode (here implemented using another
			% bridge). 2 other panels used for spacing
			hierarchy = horizontalCut(hierarchy, [0.34, 0.13, 0.06, 0.13, 0.34]);
			hierarchy.children{3} = verticalCut(hierarchy.children{3}, [0.35,0.3,0.35]);
		
			obj.acdc_converters.count = obj.acdc_converters.count + 1;
			
			bridge.name = [obj.name,'/bridge_', num2str(obj.acdc_converters.count)];
			filter.name = [obj.name,'/filter_', num2str(obj.acdc_converters.count)];
			diode.name = [obj.name,'/diode_', num2str(obj.acdc_converters.count)];
			
			add_block('powerlib/Power Electronics/Universal Bridge',...
				bridge.name, ...
				'Position', hierarchy.children{1}.rect ...
			);
			add_block('powerlib/Elements/Series RLC Branch',...
				filter.name, ...
				'Position', hierarchy.children{3}.children{2}.rect ...
			);
			add_block('powerlib/Power Electronics/Universal Bridge',...
				diode.name, ...
				'Position', hierarchy.children{5}.rect ...
			);
			
			set_param(bridge.name, 'Arms', '2');
			set_param(bridge.name, 'Device', 'Diodes');
			
			set_param(filter.name, 'BranchType', 'C');
			set_param(filter.name, 'Capacitance', '0.1'); % Farad
			set_param(filter.name, 'Orientation', 'down');
			
			set_param(diode.name, 'Arms', '2');
			set_param(diode.name, 'Device', 'Diodes');
			
			% LConn[1..2], RConn[1..2]
			bridge.hnd = get_param(bridge.name,'PortHandles');
			diode.hnd = get_param(diode.name,'PortHandles');
			% LConn (+), RConn (-)
			filter.hnd = get_param(filter.name,'PortHandles');
			
			% internal connections
			add_line(obj.name,...
				bridge.hnd.RConn(1),...
				filter.hnd.LConn,...
				'Autorouting', 'on'...
			);
			add_line(obj.name,...
				bridge.hnd.RConn(2),...
				filter.hnd.RConn,...
				'Autorouting', 'on'...
			);
			add_line(obj.name,...
				filter.hnd.LConn,...
				diode.hnd.LConn(1),...
				'Autorouting', 'on'...
			);
			add_line(obj.name,...
				filter.hnd.RConn,...
				diode.hnd.LConn(2),...
				'Autorouting', 'on'...
			);
			
			obj.acdc_converters.bridges{obj.acdc_converters.count} = bridge;
			obj.acdc_converters.filters{obj.acdc_converters.count} = filter;
			obj.acdc_converters.diodes{obj.acdc_converters.count} = diode;
		end
		
		% For consistency sake between TX and RX readings.
		function obj = newGround(obj, hierarchy)
		end
		
		function obj = newPoweredDevice(obj, hierarchy)
			obj.consumers.count = obj.consumers.count + 1;
		
			element.name = [obj.name,'/consumer_', num2str(obj.consumers.count)];

			add_block('powerlib/Elements/Series RLC Branch',...
				element.name,...
				'Position', hierarchy.rect ...
			);
			
			set_param(element.name, 'BranchType', 'R');
			set_param(element.name, 'Resistance', num2str(1000)); %ohms
			set_param(element.name, 'Orientation', 'down');
			
			% for later connection of blocks
			element.hnd = get_param(element.name,'PortHandles');
			
			obj.consumers.elements{obj.consumers.count} = element;
		end
	end
	
	% Runtime settings
	methods
		function obj = setVoltages(obj, voltage_vector)
			if length(voltage_vector) ~= obj.nt
				error('Inconsistent voltage vector');
			end
			
			for i = 1 : obj.sources.count
				set_param(obj.sources.elements{i}.name, 'Amplitude', num2str(voltage_vector(i)));
				set_param(obj.sources.elements{i}.name, 'Phase', '0');
				set_param(obj.sources.elements{i}.name, 'Frequency', num2str(obj.operating_frequency));
			end
			
		end
		
		% This function generates and applies a new set of mutual inductances. It also generates
		% the self-inductances and the capacitances in a way to keep the system resonant.
		function obj = changeCouplings(obj)
			
			% the mutual induction in a simulated homogeneous system of windings
			M = 4*pi*1e-6 * obj.coupling_helper.generateMutualInductionMatrix(obj.nt + obj.nr);
			
			% the self-induction of each coil
			L = CouplingHelper.referenceSelfInductance();
			
			% the complete inductance matrix
			inductances = -M + L * eye(obj.nt + obj.nr);
			
			% setting up the ideal inductors
			set_param(obj.mutual_coupler.name, 'InductanceMatrix', mat2str(inductances));
			set_param(obj.mutual_coupler.name, 'ResistanceMatrix', mat2str(zeros(obj.nt + obj.nr)));
			
			% achieving resonance using the capacitive reactance
			C = 1 / ( ( 2 * pi * obj.operating_frequency )^2 * L );
			for i = 1 : obj.tx_rlc_branches.count
				set_param(obj.tx_rlc_branches.components{i}.name, 'Capacitance', num2str(C));
			end
			for i = 1 : obj.rx_rlc_branches.count
				set_param(obj.rx_rlc_branches.components{i}.name, 'Capacitance', num2str(C));
			end
		end
		
		function obj = setResistances(obj, resistance_vector)
			if length(resistance_vector) ~= obj.tx_rlc_branches.count + obj.rx_rlc_branches.count
				error('Inconsistent resistance vector');
			end
			
			for i = 1 : obj.tx_rlc_branches.count
				set_param(obj.tx_rlc_branches.components{i}.name, 'Resistance', num2str(resistance_vector(i)));
			end
			
			for i = 1 : obj.rx_rlc_branches.count
				set_param(obj.rx_rlc_branches.components{i}.name, 'Resistance', num2str(resistance_vector(obj.nt + i)));
			end
		end
	end
	
	% connections
	methods
		function obj = connectTXComponents(obj)
			for i = 1 : obj.sources.count
				% connecting the source to the rlc
				add_line(obj.name,...
					obj.sources.elements{i}.hnd.LConn,...
					obj.tx_rlc_branches.components{i}.hnd.LConn,...
					'Autorouting', 'on'...
				);
				
				% connecting the source to the coil
				add_line(obj.name,...
					obj.sources.elements{i}.hnd.RConn,...
					obj.mutual_coupler.hnd.RConn(i),...
					'Autorouting', 'on'...
				);
				
				% connecting the coil to the rlc
				add_line(obj.name,...
					obj.tx_rlc_branches.components{i}.hnd.RConn,...
					obj.mutual_coupler.hnd.LConn(i),...
					'Autorouting', 'on'...
				);
			end
		end
		
		function obj = connectRXComponents(obj)
			for i = 1 : obj.nr
				
				add_line(obj.name,...
					obj.rx_rlc_branches.components{i}.hnd.RConn,...
					obj.acdc_converters.bridges{i}.hnd.LConn(2),...
					'Autorouting', 'on'...
				);
				
				add_line(obj.name,...
					obj.mutual_coupler.hnd.LConn(obj.nt + i),...
					obj.acdc_converters.bridges{i}.hnd.LConn(1),...
					'Autorouting', 'on'...
				);
				
				add_line(obj.name,...
					obj.mutual_coupler.hnd.RConn(obj.nt + i),...
					obj.rx_rlc_branches.components{i}.hnd.LConn,...
					'Autorouting', 'on'...
				);
				
				add_line(obj.name,...
					obj.acdc_converters.diodes{i}.hnd.RConn(1),...
					obj.consumers.elements{i}.hnd.LConn,...
					'Autorouting', 'on'...
				);
				
				add_line(obj.name,...
					obj.acdc_converters.diodes{i}.hnd.RConn(2),...
					obj.consumers.elements{i}.hnd.RConn,...
					'Autorouting', 'on'...
				);
				
				add_line(obj.name,...
					obj.consumers.elements{i}.hnd.LConn,...
					obj.batteries.components{i}.hnd.LConn(1),...
					'Autorouting', 'on'...
				);
				
				add_line(obj.name,...
					obj.consumers.elements{i}.hnd.RConn,...
					obj.batteries.components{i}.hnd.LConn(2),...
					'Autorouting', 'on'...
				);
				
			end
		end
	end
	
	% operational methods
	methods
		function [result, t] = run(obj)
			tic;
			result = sim(obj.name);
			t = toc;
		end
	end
end