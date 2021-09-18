classdef ImpedanceManager
	properties 
        file
	end
	
	methods
		function obj = ImpedanceManager(file)
            obj.file = file;
		end
		
        function obj = addReadings(obj, charge, internalResistance, equivalentResistance)
            
        end
	end
end
