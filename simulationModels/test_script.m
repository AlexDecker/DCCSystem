clc;
clear all;

nt = 1;
nr = 1;

model = WPTSystem('wpt', nt, nr, []);

model = model.setVoltages(30*ones(nt,1));
model = model.setResistances(ones(nt + nr,1));

%[result, t] = model.run();
%tic;
%sim('wpt');
%toc

%who

%model.destroy();