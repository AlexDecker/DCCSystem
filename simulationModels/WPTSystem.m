classdef WPTSystem
	properties
		name
		nt % number of transmitters
		nr % number of receivers
		padding % proportion of circuit-area padding (0-1)
		
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
			% local setup
			obj.coupling_helper = CouplingHelper([], false);
			obj.nt = nt;
			obj.nr = nr;
			obj.padding = 0.1; % 10%
			obj.operating_frequency = 1e+6; % 1 MHz
			
			if obj.coupling_helper.n_coils ~= obj.nt + obj.nr
				error(['nt and nr sum must be equal to ', num2str(obj.coupling_helper.n_coils)]);
			end
			
			% environment setup
			obj.name = name;
			new_system(name);
			open_system(name);
			
			% topological hierarchy
			if isempty(top_rect) || top_rect(1) <= top_rect(3) || top_rect(2) <= top_rect(4)
				disp('Using default area information...');
				top_rect = [0, 0, 1000, 175 * max(nt, nr)];
			end
			
			% block containers setup
			obj.sources.count = 0;
			obj.sources.elements = cell(obj.nt,1);
			
			obj.tx_rlc_branches.count = 0;
			obj.tx_rlc_branches.elements = cell(obj.nt,1);
			
			obj.rx_rlc_branches.count = 0;
			obj.rx_rlc_branches.elements = cell(obj.nr,1);
			
			obj.acdc_converters.count = 0;
			obj.acdc_converters.bridges = cell(obj.nr,1);
			obj.acdc_converters.filters = cell(obj.nr,1);
			obj.acdc_converters.diodes = cell(obj.nr,1);
			
			obj.consumers.count = 0;
			obj.consumers.elements = cell(obj.nr,1);
			
			obj.batteries.count = 0;
			obj.batteries.elements = cell(obj.nr,1);
			
			obj.hierarchy.rect = top_rect;
			obj.hierarchy.children = cell(0);
			
			% dividing the area between tx, coupling, and rx panels
			obj.hierarchy = obj.horizontalCut(obj.hierarchy, [0.15, 0.05, 0.8]);
			[obj.hierarchy.children{1}, obj] = obj.setupTXSide(obj.hierarchy.children{1});
			[obj.hierarchy.children{2}, obj] = obj.setupCouplingAbstraction(obj.hierarchy.children{2});
			[obj.hierarchy.children{3}, obj] = obj.setupRXSide(obj.hierarchy.children{3});
			
			% add all connections
			obj = obj.connectTXComponents()
			obj = obj.connectRXComponents()
		end
		
		function [hierarchy, obj] = setupTXSide(obj, hierarchy)
			
			% cut into equal sub areas
			distribution = ones(1, obj.nt) / obj.nt;
			hierarchy = obj.verticalCut(hierarchy, distribution);
			
			% build each circuit
			for i = 1 : length(hierarchy.children)
				hierarchy.children{i} = obj.addPadding(hierarchy.children{i});
				[hierarchy.children{i}.children{1}, obj] = obj.buildTXCircuit(hierarchy.children{i}.children{1});
			end
			
		end
		
		function [hierarchy, obj] = setupCouplingAbstraction(obj, hierarchy)
			hierarchy = obj.addPadding(hierarchy);
			obj = obj.newCoupledInductors(hierarchy.children{1});
		end
		
		function [hierarchy, obj] = setupRXSide(obj, hierarchy)
		
			% cut into equal sub areas
			distribution = ones(1, obj.nr) / obj.nr;
			hierarchy = obj.verticalCut(hierarchy, distribution);
			
			% build each circuit
			for i = 1 : length(hierarchy.children)
				hierarchy.children{i} = obj.addPadding(hierarchy.children{i});
				[hierarchy.children{i}.children{1}, obj] = obj.buildRXCircuit(hierarchy.children{i}.children{1});
			end
			
		end
		
		function [hierarchy, obj] = buildTXCircuit(obj, hierarchy)
			
			hierarchy = obj.horizontalCut(hierarchy, [0.25, 0.75]);
			
			hierarchy.children{1} = obj.verticalCut(hierarchy.children{1}, [0.2, 0.4, 0.4]);
			obj = obj.newSource(hierarchy.children{1}.children{2});
			
			hierarchy.children{2} = obj.verticalCut(hierarchy.children{2}, [0.6, 0.4]);
			obj = obj.newRC(true, hierarchy.children{2}.children{2});
			
		end
		
		function [hierarchy, obj] = buildRXCircuit(obj, hierarchy)
			% Panels: rc, acdc converter, consumer, battery
			hierarchy = obj.horizontalCut(hierarchy, [0.1,0.65,0.05,0.2]);
			
			for i = 1 : length(hierarchy.children)
				hierarchy.children{i} = obj.addPadding(hierarchy.children{i});
			end
			
			hierarchy.children{1}.children{1} = obj.verticalCut(hierarchy.children{1}.children{1}, [0.8, 0.2]);
			obj = obj.newRC(false, hierarchy.children{1}.children{1}.children{2});
			
			[hierarchy.children{2}.children{1}, obj] = obj.newACDCConverter(hierarchy.children{2}.children{1});
			
			hierarchy.children{3}.children{1} = obj.verticalCut(hierarchy.children{3}.children{1}, [0.4, 0.2, 0.4]);
			obj = obj.newPoweredDevice(hierarchy.children{3}.children{1}.children{2});
			
			obj = obj.newBattery(0, hierarchy.children{4}.children{1});
			
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
		function obj = newRC(obj, isTX, hierarchy)

			if isTX
				obj.tx_rlc_branches.count = obj.tx_rlc_branches.count + 1;
				count = obj.tx_rlc_branches.count;
				element.name = [obj.name,'/tx_rlc_', num2str(count)];
			else
				obj.rx_rlc_branches.count = obj.rx_rlc_branches.count + 1;
				count = obj.rx_rlc_branches.count;
				element.name = [obj.name,'/rx_rlc_', num2str(count)];
			end

			add_block('powerlib/Elements/Series RLC Branch',...
				element.name,...
				'Position', hierarchy.rect ...
			);
			
			set_param(element.name, 'BranchType', 'RC');
			
			% for later connection of blocks
			element.hnd = get_param(element.name,'PortHandles');
			
			if isTX
				obj.tx_rlc_branches.elements{count} = element;
			else
				obj.rx_rlc_branches.elements{count} = element;
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
		function obj = newBattery(obj, SOC, hierarchy)
			
			if SOC < 0 || SOC > 100
				error('Invalid state-of-charge.');
			end
			
			obj.batteries.count = obj.batteries.count + 1;
			
			element.name = [obj.name,'/bat_', num2str(obj.batteries.count)];
			
			add_block('electricdrivelib/Extra Sources/Battery',...
				element.name, ...
				'Position', hierarchy.rect ...
			);
			
			set_param(element.name, 'Orientation', 'right');
			set_param(element.name, 'SOC', num2str(SOC));
			
			% LConn [1..2]
			element.hnd = get_param(element.name,'PortHandles');
			
			obj.batteries.elements{obj.batteries.count} = element;
		end
		
		% Provides interface between the rx rlc ring and the DC internal circuit
		function [hierarchy, obj] = newACDCConverter(obj, hierarchy)
			
			% 3 panels: one for the bridge, the second one for the filter capacitor
			% and the third one for the isolation diode (here implemented using another
			% bridge)
			hierarchy = obj.horizontalCut(hierarchy, [0.4, 0.2, 0.4]);
			hierarchy.children{1} = obj.addPadding(hierarchy.children{1});
			hierarchy.children{2} = obj.addPadding(hierarchy.children{2});
			hierarchy.children{3} = obj.addPadding(hierarchy.children{3});
		
			obj.acdc_converters.count = obj.acdc_converters.count + 1;
			
			bridge.name = [obj.name,'/bridge_', num2str(obj.acdc_converters.count)];
			filter.name = [obj.name,'/filter_', num2str(obj.acdc_converters.count)];
			diode.name = [obj.name,'/diode_', num2str(obj.acdc_converters.count)];
			
			add_block('powerlib/Power Electronics/Universal Bridge',...
				bridge.name, ...
				'Position', hierarchy.children{1}.children{1}.rect ...
			);
			hierarchy.children{2}.children{1} = obj.verticalCut(hierarchy.children{2}.children{1}, [0.4,0.2,0.4]);
			add_block('powerlib/Elements/Series RLC Branch',...
				filter.name, ...
				'Position', hierarchy.children{2}.children{1}.children{2}.rect ...
			);
			add_block('powerlib/Power Electronics/Universal Bridge',...
				diode.name, ...
				'Position', hierarchy.children{3}.children{1}.rect ...
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
			M = 4*pi*1e-6 * obj.coupling_helper.generateMutualInductionMatrix();
			
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
				set_param(obj.tx_rlc_branches.elements{i}.name, 'Capacitance', num2str(C));
			end
			for i = 1 : obj.rx_rlc_branches.count
				set_param(obj.rx_rlc_branches.elements{i}.name, 'Capacitance', num2str(C));
			end
		end
		
		function obj = setResistances(obj, resistance_vector)
			if length(resistance_vector) ~= obj.tx_rlc_branches.count + obj.rx_rlc_branches.count
				error('Inconsistent resistance vector');
			end
			
			for i = 1 : obj.tx_rlc_branches.count
				set_param(obj.tx_rlc_branches.elements{i}.name, 'Resistance', num2str(R(i)));
			end
			
			for i = 1 : obj.rx_rlc_branches.count
				set_param(obj.rx_rlc_branches.elements{i}.name, 'Resistance', num2str(R(obj.nt + i)));
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
					obj.tx_rlc_branches.elements{i}.hnd.LConn,...
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
					obj.tx_rlc_branches.elements{i}.hnd.RConn,...
					obj.mutual_coupler.hnd.LConn(i),...
					'Autorouting', 'on'...
				);
			end
		end
		
		function obj = connectRXComponents(obj)
			for i = 1 : obj.sources.count
				
			end
		end
	end
	
	% hierarchical division of the circuit area
	methods
		function hierarchy = horizontalCut(obj, hierarchy, distribution)
			
			if sum(distribution) <= 0 || ~isempty(hierarchy.children)
				error('Invalid cut');
			end
			
			distribution = distribution / sum(distribution);
			
			% top walking dimension
			w = hierarchy.rect(3) - hierarchy.rect(1);
			
			% horizontal cut point
			w_start = hierarchy.rect(1);
			
			for i = 1 : length(distribution)
				
				% horizontal cut point (end)
				w_end = w_start + distribution(i) * w;
				
				% sub-area
				child.rect = [w_start, hierarchy.rect(2), w_end, hierarchy.rect(4)];
				child.children = cell(0);
				
				hierarchy.children{end + 1} = child;
				
				w_start = w_end;
			end
		end
		
		function hierarchy = verticalCut(obj, hierarchy, distribution)
			
			if sum(distribution) <= 0 || ~isempty(hierarchy.children)
				error('Invalid cut');
			end
			
			distribution = distribution / sum(distribution);
			
			% top walking dimension
			h = hierarchy.rect(4) - hierarchy.rect(2);
			
			% vertical cut point
			h_start = hierarchy.rect(2);
			
			for i = 1 : length(distribution)
				
				% vertical cut point (end)
				h_end = h_start + distribution(i) * h;
				
				% sub-area
				child.rect = [hierarchy.rect(1), h_start, hierarchy.rect(3), h_end];
				child.children = cell(0); % no children so far
				
				% new child-area
				hierarchy.children{end + 1} = child;
				
				h_start = h_end;
			end
		end
		
		function hierarchy = addPadding(obj, hierarchy)
			
			if ~isempty(hierarchy.children)
				error('Invalid padding insertion');
			end
			
			% top dimensions
			h = hierarchy.rect(3) - hierarchy.rect(1);
			w = hierarchy.rect(4) - hierarchy.rect(2);
					
			% sub-area
			child.rect = hierarchy.rect + [h * obj.padding / 2, w * obj.padding / 2,...
				-h * obj.padding / 2, -w * obj.padding / 2];
			child.children = cell(0); % no children so far
			
			hierarchy.children{end + 1} = child;
		end
	end
end