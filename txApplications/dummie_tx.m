%this application only generates random voltages limited by the min-max-vector 'voltage'
%using intervals of 'timeskip' seconds
classdef dummie_tx < powerTXApplication
    properties
        timeSkip
        voltage
    end
    methods
        function obj = dummie_tx(voltage, timeSkip)
            obj@powerTXApplication();%building superclass structure
            obj.voltage = voltage;
            obj.timeSkip = timeSkip;
        end

        function [obj,netManager,WPTManager] = init(obj,netManager,WPTManager)
            netManager = setTimer(obj,netManager,0,obj.timeSkip);%begin the time events
        end

        function [obj,netManager,WPTManager] = handleMessage(obj,data,GlobalTime,netManager,WPTManager)          
        end

        function [obj,netManager,WPTManager] = handleTimer(obj,GlobalTime,netManager,WPTManager)
            voltVector = zeros(WPTManager.nt,1);
            for i=1:WPTManager.nt %random voltages
                voltVector(i) = obj.voltage(1) + rand*(obj.voltage(2)-obj.voltage(1));
            end
            WPTManager = setSourceVoltages(obj,WPTManager,voltVector,GlobalTime);
            netManager = setTimer(obj,netManager,GlobalTime,obj.timeSkip);%schedule next event
        end
    end
end
