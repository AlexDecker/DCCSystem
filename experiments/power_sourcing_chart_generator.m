boxplots = false;
%variable_of_interest = 'nr';
%variable_of_interest = 'nt';
%variable_of_interest = 'timeLine_size';
%variable_of_interest = 'nSegments';
variable_of_interest = 'sample_size';
measurement = 'success';
%measurement = 'time';
%measurement = 'size_MB';

figure;
hold on;
disp('FFPareto')
parser(100, 'power_sourcing', 'exp_', 3, measurement,variable_of_interest,boxplots,'-+')
disp('FFDoubleDummie')
parser(100, 'power_sourcing', 'exp_', 2, measurement,variable_of_interest,boxplots,'-o')
disp('FFGreedy')
parser(100, 'power_sourcing', 'exp_', 4, measurement,variable_of_interest,boxplots,'--o')
disp('FFRobin')
parser(100, 'power_sourcing', 'exp_', 5, measurement,variable_of_interest,boxplots,':o')
disp('FFDummie')
parser(100, 'power_sourcing', 'exp_', 1, measurement,variable_of_interest,boxplots,'-.o')
ylabel('Success (0-1)')
%ylabel('Success (absolute)')
%ylabel('Time (s)')
%ylabel('Space (MB)')
%xlabel('Number of passive circuits')
%xlabel('Number of active circuits')
%xlabel('Number of time-slots')
%xlabel('Number of intervals for charge discretization')
xlabel('Sample size (for instance generation)')
legend('Pareto','Simple', 'Max Sum', 'Max Min', 'Fly-weight')
