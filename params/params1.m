%Parammeter function for dummie application
function params = params1()
	params.voltage = [5,12];%voltage range (V)
	params.timeSkip = 60;%time interval between voltage updates (s)
	params.Pmax = 500;%max active power to be spent
end
