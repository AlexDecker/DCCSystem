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
			
			% solver definitions: default
			add_block('powerlib/powergui',[name, '/powergui']);
			%set_param([name, '/powergui'],'SimulationMode','Continuous');
			set_param([name, '/powergui'],'SimulationMode','Discrete');
			set_param([name, '/powergui'],'SampleTime','1e-6');
			set_param('wpt', 'MaxStep', '1');
			
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
		function obj = setVoltages(obj, voltage_vector)
		end
		
		function obj = setResistances(obj, resistance_vector)
		end
		
		function obj = updateParameters(obj)
		end
		
		function settings = getSettings(obj)
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
		function [result, t] = run(obj)
			tic;
			result = sim(obj.systemName);
			t = toc;
		end
	end
end