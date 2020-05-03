boxplots = false;
%variable_of_interest = 'nr';
%variable_of_interest = 'nt';
%variable_of_interest = 'timeLine_size';
%variable_of_interest = 'nSegments';
variable_of_interest = 'sample_size';

figure;
hold on;
disp('FFPareto')
parser_charging(100, 'charging', 'exp_charging_', 3, variable_of_interest,'-+')
disp('FFDoubleDummie')
parser_charging(100, 'charging', 'exp_charging_', 2, variable_of_interest,'-o')
disp('FFGreedy')
parser_charging(100, 'charging', 'exp_charging_', 4, variable_of_interest,'--o')
disp('FFRobin')
parser_charging(100, 'charging', 'exp_charging_', 5, variable_of_interest,':o')
disp('FFDummie')
parser_charging(100, 'charging', 'exp_charging_', 1, variable_of_interest,'-.o')
ylabel('Normalized charging time (0-1)')
%xlabel('Number of passive circuits')
%xlabel('Number of active circuits')
%xlabel('Number of time-slots')
%xlabel('Number of intervals for charge discretization')
xlabel('Sample size (for instance generation)')
legend('Pareto','Simple', 'Max Sum', 'Max Min', 'Fly-weight')
