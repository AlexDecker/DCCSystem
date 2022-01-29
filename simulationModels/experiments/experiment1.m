% This script aims on gather data from a 1TX1RX wpt circuit, in order to model
% the set converter/battery/consumer in terms only of the state-of-charge and the
% consuming current (considering the other parameters of the battery are constant)

wpt = WPTSystem('wpt', 1, 1, []);

% Change voltages, resistances, couplings and SOC
wpt = wpt.updateParameters();

wpt.run();

wpt.getSettings();

wpt.destroy();