classdef InstanceGenerator
	properties
		% outputs
		solution
		chargeData
		constraints
		timeLine
		dt
		
		% inputs
		deviceData
		nt
		nSegments
		timeLine_size
		
	end
	
	methods
		
		function obj = InstanceGenerator(deviceData, nt, nSegments, timeLine_size)
			obj.deviceData = deviceData;
			obj.nt = nt;
			obj.nSegments = nSegments;
			obj.timeLine_size = timeLine_size;
			obj = obj.generateSimpleInstance();
		end
		
		function  obj = generateSimpleInstance(obj)
		end
		
		function obj = generateSampledInstance(obj, nSamples)
		end
		
	end
end