%Builds both rx and tx applications based on configurations from an informed file (paramFile) and
%using appName as model

%appName: model filename (sample_rx.m,sample_tx.m -> appName='sample')
%paramFile: filename of the file containing the parameters of both rx and tx (path/param1.m -> 'param1')

function [powerTX, powerRX] = appBuilder(appName, paramFile, NRX)
	param = eval(paramFile);%execute the function in paramFile and get the parammeters
	switch(appName)
		case 'dummie'
			%this application only generates random voltages limited by the min-max-vector 'voltage'
			%using intervals of 'timeskip' seconds
			powerTX = dummie_tx(param.voltage, param.timeSkip);
			powerRX = [];
			for i=1:NRX
				powerRX = [powerRX, struct('obj',powerRXApplication(i))];
			end
		case 'optimum'
			%this application calculates and applies the optimum voltages for a given timeSkip
			powerTX = optimum_tx(param.Pmax,param.timeSkip);
			powerRX = [];
			for i=1:NRX
				powerRX = [powerRX, struct('obj',powerRXApplication(i))];
			end
		otherwise
			disp('Unknown application');
	end
end
