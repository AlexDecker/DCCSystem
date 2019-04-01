%this application generates the optimum voltages in order to minimize the time required to finish
%all the charging processes given a 'timeSkip' discretization interval
classdef optimum_tx < powerTXApplication
    properties
        timeSkip
		count
        voltage_progression%matriz in which each column corresponds to the vt vector at each moment
    end
    methods
        function obj = optimum_tx(timeSkip)
            obj@powerTXApplication();%building superclass structure
            obj.timeSkip = timeSkip;
			obj.count = 1;
        end

        function [obj,netManager,WPTManager] = init(obj,netManager,WPTManager)
			%calculate the optimization voltages
			obj.voltage_progression = calc_optimum_voltages(WPTManager,obj.timeSkip);
			%apply the calculated voltages
            [obj,WPTManager] = applyVoltages(obj,0,WPTManager);
			netManager = setTimer(obj,netManager,0,obj.timeSkip);%begin the time events
        end

        function [obj,netManager,WPTManager] = handleMessage(obj,data,GlobalTime,netManager,WPTManager)          
        end

        function [obj,netManager,WPTManager] = handleTimer(obj,GlobalTime,netManager,WPTManager)
            [obj,WPTManager] = applyVoltages(obj,GlobalTime,WPTManager);
			netManager = setTimer(obj,netManager,GlobalTime,obj.timeSkip);%schedule next event
        end

		function [obj,WPTManager] = applyVoltages(obj,GlobalTime,WPTManager)
			s = size(obj.voltage_progression);
			if(s(2)>=obj.count)%if all the voltages have been applyied
				WPTManager = setSourceVoltages(obj,WPTManager,obj.voltage_progression(:,obj.count),GlobalTime);
				obj.count = obj.count+1;
			else
				obj = endSimulation(obj);
			end
		end
    end
end
