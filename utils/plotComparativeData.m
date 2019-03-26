%given a cell array with the devData filenames (see generateData.m), this function plots
%a comparative chart for the SOC progression, for the battery voltage progression and for
%the charging current progression.
function plotComparativeData(fileStructArray)
	%loading all data
	nRX = 0;
	devDataList = {};
	for f=1:length(fileStructArray)
		devData=[];load(fileStructArray{f});%load each file
		devDataList{end+1} = devData;%append to the general list
		nRX = max(nRX,length(devData));%the maximum number of RX
	end

	%SOC progression
	for d=1:nRX
		figure; hold on;
		for f=1:length(fileStructArray)
			devData = devDataList{f};
			if(length(devData)>=d)
				plot(devData(d).SOC(2,:)/60,100*devData(d).SOC(1,:));
			end
		end
		%legend(fileStructArray);
		xlabel('Time (min)');
		ylabel('SOC (%)');
		title(['Device ',num2str(d)]);
	end
	
	%Voltage progression
	for d=1:nRX
		figure; hold on;
		for f=1:length(fileStructArray)
			devData = devDataList{f};
			if(length(devData)>=d)
				plot(devData(d).VB(2,:)/60,devData(d).VB(1,:));
			end
		end
		%legend(fileStructArray);
		xlabel('Time (min)');
		ylabel('Battery Voltage (V)');
		title(['Device ',num2str(d)]);
	end
	
	%Current progression
	for d=1:nRX
		figure; hold on;
		for f=1:length(fileStructArray)
			devData = devDataList{f};
			if(length(devData)>=d)
				plot(devData(d).CC(2,:)/60,1000*devData(d).CC(1,:));
			end
		end
		%legend(fileStructArray);
		xlabel('Time (min)');
		ylabel('Charge Current (mA)');
		title(['Device ',num2str(d)]);
	end	
end
