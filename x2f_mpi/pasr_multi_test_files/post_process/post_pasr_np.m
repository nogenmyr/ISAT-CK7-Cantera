% script to plot basic statistics for pasr_multi from the files pasr_#.out
%
% specify: np       (number of processes)
%          run_path (location of output files isat_1_#.op)

clear all
np = 32;
run_path = '../test_case_output/';
xscale_type = 'linear'; % 'log' or 'linear'
yscale_type = 'linear'; % 'log' or 'linear'

figure();
hold on;
box on;

cmap = jet(64);
for i = 1:np
    clr = cmap(1+floor(63*(np-i)/np),:);
    x = load([run_path '/pasr_' num2str(i) '.out']);
    
    istep = x(:,1);
    t = x(:,2);
    temp_mean = x(:,3);
    x1_mean = x(:,6);
    x3_mean = x(:,7);
    x6_mean = x(:,8);
    max_t = max(t);
    
    % Temperature
    subplot(2,2,1); box on; hold on;
    plot(t,temp_mean,'b','linewidth',2,'color',clr);
    set(gca,'xscale',xscale_type,'yscale',yscale_type);
    xlim([0 max_t]);
    xlabel('Time (s)');
    ylabel('<T> (K)');
    
    % Species 1
    subplot(2,2,2); box on; hold on;
    plot(t,x1_mean,'b','linewidth',2,'color',clr);
    set(gca,'xscale',xscale_type,'yscale',yscale_type);
    xlim([0 max_t]);
    xlabel('Time (s)');
    ylabel('<z_1>');
    
    % Species 3
    subplot(2,2,3); box on; hold on;
    plot(t,x3_mean,'b','linewidth',2,'color',clr);
    set(gca,'xscale',xscale_type,'yscale',yscale_type);
    xlim([0 max_t]);
    xlabel('Time (s)');
    ylabel('<z_3>');
    
    % Species 6
    subplot(2,2,4); box on; hold on;
    plot(t,x6_mean,'b','linewidth',2,'color',clr);
    set(gca,'xscale',xscale_type,'yscale',yscale_type);
    xlim([0 max_t]);
    xlabel('Time (s)');
    ylabel('<z_6>');
end
legend([num2str(transpose([1:np]))],-1);
