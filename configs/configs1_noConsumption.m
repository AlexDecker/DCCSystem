function genSettings = configs1()
    genSettings.NTX = 6;%number of transmitters
    genSettings.NRX = 3;%number of receivers
    genSettings.W = 1e6;%operational frequency
    genSettings.R = 0.5*ones(9,1);%RLCs resistance
    genSettings.C = -1*ones(9,1);%RLC's capacitance (-1=use data from .mat)
    genSettings.MAX_ACT_POWER = 200;%W
    genSettings.MAX_APP_POWER = 20000;%W
    genSettings.TOTAL_TIME = 18000;%seconds of simulation (virtual time)=5h
    genSettings.power_m = 0;%consumption (W)
    genSettings.minV = 2.3;%minimal voltage to operate
    genSettings.minVTO = 3.3;%minimal voltage to turn on
    genSettings.efficiency = 0.95;%ACDC conversion
    genSettings.envFolder = 'mobData1';%the folder containing the envList files
end
