%Collects the data from nRepetitions of the specified simulation, considering the specified
%instances 
function generateData(appName, paramFile, configFile, fileInstances, nRepetitions)
	output = ['data_',appName,'_',paramFile,'_',configFile,...
		'_nRepeat_',num2str(nRepetitions)];
	
	conf = eval(configFile);%getting the number of devices
	%initializing the structure for each device
	devData = [];
	for d=1:conf.NRX
		devData = [devData,struct('SOC',[],'VB',[],'CC',[])];
	end

	for i=fileInstances
		for n=1:nRepetitions
			LOG = simulate(2018,appName,paramFile,configFile,i,false);
			for d=1:length(LOG)
				if(isempty(devData(d).CC))%first data set
					devData(d).CC = LOG(d).CC;
					plot(LOG(d).CC(2,:),LOG(d).CC(1,:));
					figure;
					devData(d).VB = LOG(d).VB;
					plot(LOG(d).VB(2,:),LOG(d).VB(1,:));
					figure;
					devData(d).SOC = LOG(d).SOC;
					plot(LOG(d).SOC(2,:),LOG(d).SOC(1,:));
				else
					%adjust the temporal series and then sum the new data and the accumulator	
					CC = LOG(d).CC;
					devData(d).CC = adjustAndSum(CC,devData(d).CC);
				
					VB = LOG(d).VB;
					devData(d).VB = adjustAndSum(VB,devData(d).VB);
				
					SOC = LOG(d).SOC;
					devData(d).SOC = adjustAndSum(SOC,devData(d).SOC);
				end
			end
		end
	end
	
	num = length(fileInstances)*nRepetitions;

	for d=1:conf.NRX
		devData(d).CC(:,1) = devData(d).CC(:,1)./num;
		devData(d).VB(:,1) = devData(d).VB(:,1)./num;
		devData(d).SOC(:,1) = devData(d).SOC(:,1)./num;
	end

	save(output, 'devData');
end
