%Collects the data from nRepetitions of the specified simulation, considering the specified
%instances 
function generateData(appName, paramFile, configFile, fileInstances, nRepetitions)
    %updating index.txt with the information regarding 
    [cod, line] = system('tail -n 1 out/index.txt');
    if cod ~= 0 %success of tail command
        error('Could not read index.txt');
    end
    lineCell = strsplit(line, ';');
    if length(lineCell)~=6 %6 fields
        error('Unexpected line in index.txt');
    end
    fileNumber = str2num(lineCell{1});%get the first field
    if isempty(fileNumber)%not numeric
        fileNumber = 0;
    end
    fileNumber = fileNumber + 1;%this file will use the next available file number
    %open in append mode
    index_fp = fopen('out/index.txt','a+');
	fprintf(index_fp, [num2str(fileNumber),'; ',appName,'; ',paramFile,'; ',...
        configFile,'; |',sprintf('%d|',fileInstances),'; ',num2str(nRepetitions),'\n']);
    fclose(index_fp);
    
	%Getting stuff to execute the simulation
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
					devData(d).VB = LOG(d).VB;
					devData(d).SOC = LOG(d).SOC;
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
		devData(d).CC(1,:) = devData(d).CC(1,:)./num;
		devData(d).VB(1,:) = devData(d).VB(1,:)./num;
		devData(d).SOC(1,:) = devData(d).SOC(1,:)./num;
	end
	
    save(['out/file',num2str(fileNumber),'.mat', 'devData');
    
end
