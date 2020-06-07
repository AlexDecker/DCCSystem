boxplots = false;

[x_nr_MS,y_nr_MS,err_nr_MS, x_nt_MS,y_nt_MS,err_nt_MS, x_t_MS,y_t_MS,err_t_MS,success_vector] = solve_and_plot_charging(100, 'charging', 'exp_charging_', 5, ones(100,1));

use_vector = success_vector;
%use_vector = ones(100,1);

disp('FFGreedyCurrents');
[x_nr_GC,y_nr_GC,err_nr_GC, x_nt_GC,y_nt_GC,err_nt_GC, x_t_GC,y_t_GC,err_t_GC,~] = solve_and_plot_charging(100, 'charging', 'exp_charging_', 3, use_vector);
disp('FFMaxPower');
[x_nr_MP,y_nr_MP,err_nr_MP, x_nt_MP,y_nt_MP,err_nt_MP, x_t_MP,y_t_MP,err_t_MP,~] = solve_and_plot_charging(100, 'charging', 'exp_charging_', 4, use_vector);
disp('FFMultSpot');


variable_of_interest = 'nr';
figure;
hold on;
errorbar(x_nr_GC,y_nr_GC,err_nr_GC,'-+','linewidth',3);
disp('FFGreedy')
parser_charging(100, 'charging', 'exp_charging_', 4, variable_of_interest,'--o',use_vector)
disp('FFRobin')
parser_charging(100, 'charging', 'exp_charging_', 5, variable_of_interest,':o',use_vector)
errorbar(x_nr_MP,y_nr_MP,err_nr_MP,'-o','linewidth',3);
errorbar(x_nr_MS,y_nr_MS,err_nr_MS,'-.o','linewidth',3);
disp('FFPareto')
parser_charging(100, 'charging', 'exp_charging_', 3, variable_of_interest,'--+',use_vector)
ylabel('Normalized charging time (0-1)')
xlabel('Number of passive circuits')
legend('Max-Sum-Of-Currents','Max Sum', 'Max Min', 'MaxPower', 'MultiSpot', 'Pareto')

variable_of_interest = 'nt';
figure;
hold on;
errorbar(x_nt_GC,y_nt_GC,err_nt_GC,'-+','linewidth',3);
disp('FFGreedy')
parser_charging(100, 'charging', 'exp_charging_', 4, variable_of_interest,'--o',use_vector)
disp('FFRobin')
parser_charging(100, 'charging', 'exp_charging_', 5, variable_of_interest,':o',use_vector)
errorbar(x_nt_MP,y_nt_MP,err_nt_MP,'-o','linewidth',3);
errorbar(x_nt_MS,y_nt_MS,err_nt_MS,'-.o','linewidth',3);
disp('FFPareto')
parser_charging(100, 'charging', 'exp_charging_', 3, variable_of_interest,'--+',use_vector)
ylabel('Normalized charging time (0-1)')
xlabel('Number of active circuits')
legend('Max-Sum-Of-Currents','Max Sum', 'Max Min', 'MaxPower', 'MultiSpot', 'Pareto')

variable_of_interest = 'timeLine_size';
figure;
hold on;
errorbar(x_t_GC,y_t_GC,err_t_GC,'-+','linewidth',3);
disp('FFGreedy')
parser_charging(100, 'charging', 'exp_charging_', 4, variable_of_interest,'--o',use_vector)
disp('FFRobin')
parser_charging(100, 'charging', 'exp_charging_', 5, variable_of_interest,':o',use_vector)
errorbar(x_t_MP,y_t_MP,err_t_MP,'-o','linewidth',3);
errorbar(x_t_MS,y_t_MS,err_t_MS,'-.o','linewidth',3);
disp('FFPareto')
parser_charging(100, 'charging', 'exp_charging_', 3, variable_of_interest,'--+',use_vector)
ylabel('Normalized charging time (0-1)')
xlabel('Number of time-slots')
legend('Max-Sum-Of-Currents','Max Sum', 'Max Min', 'MaxPower', 'MultiSpot', 'Pareto')