classdef WPTSystem
	properties
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
		diodes
		acdc_converters
		batteries
		
		% coupling model
		mutual_coupler % block
		coupling_helper % mutual inductance generator
		
		operating_frequency
	end
	
	% setup methods
	methods
		function obj = WPTSystem(nt, nr, top_rect)
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
			new_system('wpt_system');
			open_system('wpt_system');
			
			% topological hierarchy
			if isempty(top_rect) || top_rect(1) <= top_rect(3) || top_rect(2) <= top_rect(4)
				disp('Using default area information...');
				top_rect = [0, 0, 100 * max(nt, nr), 1000];
			end
			
			% block containers setup
			obj.sources.count = 0; obj.sources.elements = cell(obj.nt,1);
			obj.tx_rlc_branches.count = 0; obj.tx_rlc_branches.elements = cell(obj.nt,1);
			obj.rx_rlc_branches.count = 0; obj.rx_rlc_branches.elements = cell(obj.nr,1);
			obj.diodes.count = 0; obj.diodes.elements = cell(obj.nr,1);
			obj.acdc_converters.count = 0; obj.acdc_converters.elements = cell(obj.nr,1);
			obj.batteries.count = 0; obj.batteries.elements = cell(obj.nr,1);
			
			obj.hierarchy.rect = top_rect;
			obj.hierarchy.children = cell(0);
			
			% dividing the area between tx, coupling, and rx panels
			obj.hierarchy = obj.horizontalCut(obj.hierarchy, [0.3, 0.2, 0.5]);
			[obj.hierarchy.children{1}, obj] = obj.setupTXSide(obj.hierarchy.children{1});
			[obj.hierarchy.children{2}, obj] = obj.setupCouplingAbstraction(obj.hierarchy.children{2});
			[obj.hierarchy.children{3}, obj] = obj.setupRXSide(obj.hierarchy.children{3});
			
			% add all connections
			obj = obj.connectTXComponents()
			obj = obj.connectRXComponents()
			obj = obj.connectCoupler()
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
			obj = obj.newRC(hierarchy.children{2}.children{2});
		end
		
		function [hierarchy, obj] = buildRXCircuit(obj, hierarchy)
			
		end
				
		function obj = destroy(obj)
			% close inconditionally
			bdclose('wpt_system');
		end
	end
	
	% Block creation
	methods
		function obj = newSource(obj, hierarchy)
		
			obj.sources.count = obj.sources.count + 1;
			
			element.name = ['wpt_system/src_', num2str(obj.sources.count)];
			
			add_block('powerlib/Electrical Sources/AC Voltage Source',...
				element.name,...
				'Position', hierarchy.rect ...
			);
				
			% for later connection of blocks
			element.hnd = get_param(element.name, 'PortHandles');
			
			obj.sources.elements{obj.sources.count} = element;
			
		end
		
		% used in both TX and RX (one for each circuit)
		function obj = newRC(obj, hierarchy)

			obj.tx_rlc_branches.count = obj.tx_rlc_branches.count + 1;
		
			element.name = ['wpt_system/tx_rlc_', obj.tx_rlc_branches.count];

			add_block('powerlib/Elements/Series RLC Branch',...
				element.name,...
				'Position', hierarchy.rect ...
			);
			
			set_param(element.name, 'BranchType', 'RC');
			
			% for later connection of blocks
			element.hnd = get_param(element.name,'PortHandles');
			
			obj.tx_rlc_branches.elements{obj.tx_rlc_branches.count} = element;
			
		end
		
		% https://www.mathworks.com/help/physmod/sps/powersys/ref/mutualinductance.html
		function obj = newCoupledInductors(obj, hierarchy)

			element.name = ['wpt_system/coupler'];

			add_block('powerlib/Elements/Mutual Inductance',...
				element.name, ...
				'Position', hierarchy.rect ...
			);
			
			set_param(element.name, 'TypeOfMutual', 'Generalized mutual inductance');
			set_param(element.name, 'NumberOfWindings', obj.nt + obj.nr);
			
			% first set of internal parameters
			obj = obj.changeCouplings();
			
			% for later connection of blocks
			element.hnd = get_param(element.name,'PortHandles');
			
			obj.mutual_coupler = element;
		end
		
		% Chargeable battery which supplies the powered device. 
		function obj = newBattery(obj, hierarchy)
		end
		
		% Provides interface between the rx rlc ring and the DC internal circuit
		function obj = newACDCConverter(obj, hierarchy)
		end
		
		% Avoids the battery power to interfere in the rlc behavior.
		function obj = newDiode(obj, hierarchy)
		end
		
		% For consistency sake between TX and RX readings.
		function obj = newGround(obj, hierarchy)
		end
		
		function obj = newPoweredDevice(obj, hierarchy)
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
			M = obj.CouplingHelper.generateMutualInductionMatrix();
			
			% the self-induction of each coil
			L = obj.CouplingHelper.referenceSelfInductance();
			
			% the complete inductance matrix
			inductances = -M + L * eye(obj.nt + obj.nr);
			
			% setting up the ideal inductors
			set_param(obj.coupler.name, 'InductanceMatrix', mat2str(M));
			set_param(obj.coupler.name, 'ResistanceMatrix', mat2str(zeros(obj.nt + obj.nr)));
			
			% achieving resonance using the capacitive reactance
			C = 1 / ( ( 2 * pi * obj.operating_frequency )^2 * L );
			for i = 1 : length(obj.tx_rlc)
				set_param(obj.tx_rlc{i}.name, 'Capacitance', num2str(C));
			end
			for i = 1 : length(obj.rx_rlc)
				set_param(obj.rx_rlc{i}.name, 'Capacitance', num2str(C));
			end
		end
		
		function obj = setResistances(obj, resistance_vector)
			if length(resistance_vector) ~= obj.nt + obj.nr
				error('Inconsistent resistance vector');
			end
			
			for i = 1 : obj.nt
				set_param(obj.tx_rlc{i}.name, 'Resistance', num2str(R(i)));
			end
			
			for i = 1 : obj.nr
				set_param(obj.rx_rlc{i}.name, 'Resistance', num2str(R(obj.nt + i)));
			end
		end
	end
	
	% connections
	methods
		function obj = connectTXComponents(obj)
			for i = 1 : obj.sources.count
				% connecting the source to the rlc
				add_line('wpt_system',...
					obj.sources.elements{i}.hnd.LConn,...
					obj.tx_rlc_branches.elements{i}.hnd.LConn,...
					'Autorouting', 'on'...
				);
				
				% connecting the source to the coil
				
				% connecting the coil to the rlc
			end
		end
		
		function obj = connectRXComponents(obj)
			% connecting the coil to the rlc
			% connecting the rlc to the ground
			% connecting the coil to the ACDC converter
			% connecting the ACDC converter to the diode
			% connecting the ACDC converter to the ground
			% connecting the diode to the battery
			% connecting the battery to the ground
			% connecting the battery to the device
			% connecting the device to the ground
		end
		
		function obj = connectCoupler(obj)
			% connecting the tx rlc to the coil
			% connecting the coil to the source
			% connecting the rx rlc to the coil
			% connecting the coil to the source
		end
	end
	
	% hierarchical division of the circuit area
	methods
		function hierarchy = horizontalCut(obj, hierarchy, distribution)
			
			if sum(distribution) ~= 1 || ~isempty(hierarchy.children)
				error('Invalid cut');
			end
			
			% top walking dimension
			w = hierarchy.rect(4) - hierarchy.rect(2);
			
			% horizontal cut point
			w_start = hierarchy.rect(2);
			
			for i = 1 : length(distribution)
				
				% horizontal cut point (end)
				w_end = w_start + distribution(i) * w;
				
				% sub-area
				child.rect = [hierarchy.rect(1), w_start, hierarchy.rect(3), w_end];
				child.children = cell(0);
				
				hierarchy.children{end + 1} = child;
				
				w_start = w_end;
			end
		end
		
		function hierarchy = verticalCut(obj, hierarchy, distribution)
			
			if sum(distribution) ~= 1 || ~isempty(hierarchy.children)
				error('Invalid cut');
			end
			
			% top walking dimension
			h = hierarchy.rect(3) - hierarchy.rect(1);
			
			% vertical cut point
			h_start = hierarchy.rect(1);
			
			for i = 1 : length(distribution)
				
				% vertical cut point (end)
				h_end = h_start + distribution(i) * h;
				
				% sub-area
				child.rect = [h_start, hierarchy.rect(2), h_end, hierarchy.rect(4)];
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