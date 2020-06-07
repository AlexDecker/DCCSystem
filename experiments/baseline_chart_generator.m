boxplots = false;

measurement = 'success';

disp('FFMaxPower')
[x_nr4,y_nr4,err_nr4, x_nt4,y_nt4,err_nt4, x_t4,y_t4,err_t4] = solve_and_plot(351, 'baselines', 'exp_', 4);
disp('FFMultiSpot')
[x_nr5,y_nr5,err_nr5, x_nt5,y_nt5,err_nt5, x_t5,y_t5,err_t5] = solve_and_plot(351, 'baselines', 'exp_', 5);

variable_of_interest = 'nr';
figure;
hold on;
disp('FFGreedyCurrents')
parser(351, 'baselines', 'exp_', 1, measurement,variable_of_interest,boxplots,'-+')
disp('FFGreedy')
parser(351, 'baselines', 'exp_', 2, measurement,variable_of_interest,boxplots,'--o')
disp('FFRobin')
parser(351, 'baselines', 'exp_', 3, measurement,variable_of_interest,boxplots,':o')
errorbar(x_nr4,y_nr4,err_nr4,'-o', 'linewidth', 3);
errorbar(x_nr5,y_nr5,err_nr5,'-.o', 'linewidth', 3);
legend('Max Sum (currents)', 'Max Sum', 'Max Min','Max Power', 'MultiSpot')
xlabel('Number of passive circuits')
ylabel('Success (0-1)')

variable_of_interest = 'nt';
figure;
hold on;
disp('FFGreedyCurrents')
parser(351, 'baselines', 'exp_', 1, measurement,variable_of_interest,boxplots,'-+')
disp('FFGreedy')
parser(351, 'baselines', 'exp_', 2, measurement,variable_of_interest,boxplots,'--o')
disp('FFRobin')
parser(351, 'baselines', 'exp_', 3, measurement,variable_of_interest,boxplots,':o')
errorbar(x_nt4,y_nt4,err_nt4,'-o', 'linewidth', 3);
errorbar(x_nt5,y_nt5,err_nt5,'-.o', 'linewidth', 3);
legend('Max Sum (currents)', 'Max Sum', 'Max Min','Max Power', 'MultiSpot')
xlabel('Number of active circuits')
ylabel('Success (0-1)')

variable_of_interest = 'timeLine_size';
figure;
hold on;
disp('FFGreedyCurrents')
parser(351, 'baselines', 'exp_', 1, measurement,variable_of_interest,boxplots,'-+')
disp('FFGreedy')
parser(351, 'baselines', 'exp_', 2, measurement,variable_of_interest,boxplots,'--o')
disp('FFRobin')
parser(351, 'baselines', 'exp_', 3, measurement,variable_of_interest,boxplots,':o')
errorbar(x_t4,y_t4,err_t4,'-o', 'linewidth', 3);
errorbar(x_t5,y_t5,err_t5,'-.o', 'linewidth', 3);
legend('Max Sum (currents)', 'Max Sum', 'Max Min','Max Power', 'MultiSpot')
xlabel('Number of time-slots')
ylabel('Success (0-1)')