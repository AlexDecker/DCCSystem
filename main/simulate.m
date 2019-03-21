%-appName, paramFile: see appBuilder
%-configFile: string contaning the name of the function in configs/ folder with
%the desired set of parameters
%-fileInstance: number of the correct envList file inside the folder defined in 
%configFile
%version: year of your matlab version (numeric)
function simulate(version, appName, paramFile, configFile, fileInstance)
    disp('Reminding: Please be sure that the workspace is clean (use clear all)');
    
    conf = eval(configFile);

    %GENERAL ASPECTS
    NTX = conf.NTX; %number of transmitting coils
    NRX = conf.NRX; %number of receiving coils
    W = conf.W;
    R = conf.R;%RLCs resistance
    C = conf.C;%RLC's capacitance (-1=use data from .mat)
    MAX_ACT_POWER = conf.MAX_ACT_POWER;%W
    MAX_APP_POWER = conf.MAX_APP_POWER;%W
    TOTAL_TIME = conf.TOTAL_TIME;%seconds of simulation (virtual time)

    %BATTERY (iPhone 4s battery. We didn't find some values, so LIR18650 was used
    %as refernce for these)
    fase1Limit = 0.9;          % (90%, according to the disturb point at Rl curve)
    limitToBegin = 0.93;       % (93%)
    constantCurrent_min = 0.1; % (A)
    constantCurrent_max = 1.5;   % (A)
    constantVoltage = 4.2;     % (V)
    Rc = -1;      % (ohm. -1=calculate automatically)
    Rd = -1;       % (ohm. -1=calculate automatically)
    R_MAX = 1e7;   % (ohm. Used in oder to make matrix inversion easier)
    Q0 = 0;       % initial charge (As)
    Qmax = 5148;  % (As), equivalent to 1430 mAh (iPhone 4s battery)

    bat = magMIMOLinearBattery('magMIMOLinearBattery_data.txt','Li_Ion_Battery_LIR18650.txt',...
                      Rc,Rd,Q0,Qmax,R_MAX,fase1Limit,constantCurrent_min,constantCurrent_max,...
                      constantVoltage,limitToBegin,false);
    %DEVICE
    power_m = conf.power_m; % (W)
    power_sd = 0;
    minV = conf.minV;     % (V)
    minVTO = conf.minVTO;   % (V)
    err = 0.05;     % (5%)
    efficiency = conf.efficiency; % (95% efficiency for AC/DC conversion)

    STEP=0.2;     % (s)

    dev = genericDeviceWithBattery(bat,power_m,power_sd,minV,minVTO,err,efficiency);
	DEVICE_LIST = [];
	for i=1:NRX
    	DEVICE_LIST = [DEVICE_LIST, struct('obj',dev)];
	end

    %APPLICATIONS
    [powerTX, powerRX] = appBuilder(appName, paramFile, NRX);

    %SIMULATOR

    IFACTOR=1.5;
    DFACTOR=2;
    INIT_VEL=0.01;
    MAX_ERR = 0.005;

    SHOW_PROGRESS = true;

    B_SWIPT = 0.5;%minimum SINR for the message to be undertood
    B_RF = 0.5;%minimum SINR for the message to be undertood
    A_RF = 2;%expoent for free-space path loss (RF only)
    N_SWIPT = 0.1;%Noise for SWIPT (W)
    N_RF = 0.1;%Noise for RF (W)

    [~,LOG_dev_list,LOG_app_list] = Simulate([conf.envFolder,'/',num2str(fileInstance),'.mat'],...
		NTX,R,C,W,TOTAL_TIME,MAX_ERR,R_MAX,IFACTOR,DFACTOR,INIT_VEL,MAX_ACT_POWER,...
		MAX_APP_POWER,DEVICE_LIST,STEP,SHOW_PROGRESS,powerTX,powerRX,B_SWIPT,B_RF,...
		A_RF,N_SWIPT,N_RF);

    %VISUALIZATION
        
    for i=1:length(LOG_dev_list)
        LOG = endDataAquisition(LOG_dev_list(i));
        if(version <= 2010)
            plotBatteryChart2010(LOG);
        else
            plotBatteryChart(LOG);
        end
    end

end
