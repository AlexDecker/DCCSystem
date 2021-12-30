classdef Rectifier
	properties (Access = private)
        name
        hnd
	end
	methods
		function obj = Rectifier(systemName, hierarchy, index)
			obj.name = [systemName,'/bridge_', num2str(index)];
            add_block('powerlib/Power Electronics/Universal Bridge',...
				obj.name, ...
				'Position', hierarchy.rect ...
			);
            set_param(obj.name, 'Arms', '2');
			set_param(obj.name, 'Device', 'Diodes');
            obj.hnd = get_param(obj.name,'PortHandles');
		end
		
		function hnd = inPositiveHandler(obj)
			hnd = obj.hnd.RConn(1);
		end
		
		function hnd = inNegativeHandler(obj)
			hnd = obj.hnd.RConn(2);
		end
		
		function hnd = outNegativeHandler(obj)
			hnd = obj.hnd.LConn(1);
		end
		
        function hnd = outPositiveHandler(obj)
			hnd = obj.hnd.LConn(2);
		end
	end
end