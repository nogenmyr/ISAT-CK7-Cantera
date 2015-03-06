% For two runs of pasr_multi, calculate statistics based on the files
% pasr_#.out and plot the differences
%
% specify: np         (number of processes)
%          run_path   (location of output files isat_1_#.op for run 1)
%          run_2_path (location of output files isat_1_#.op for run 2)

%clear all
np = 32;
run_path = '../test_case_output/';
run_2_path = '../test_case_output2/';
xscale_type = 'linear'; % 'log' or 'linear'
yscale_type = 'linear'; % 'log' or 'linear'
x3max = 0;

figure();
hold on;
box on;

cmap = jet(64);
for i = 1:np
    clr = cmap(1+floor(63*(np-i)/np),:);
    x = load([run_path '/pasr_' num2str(i) '.out']);
    y = load([run_2_path '/pasr_' num2str(i) '.out']);

    istep = x(:,1);
    t = x(:,2);
    temp_mean = x(:,3);
    x1_mean = x(:,6);
    x3_mean = x(:,7);
    x6_mean = x(:,8);
    ytemp_mean = y(:,3);
    y1_mean = y(:,6);
    y3_mean = y(:,7);
    y6_mean = y(:,8);
    max_t = max(t);
    
    % Temperature
    subplot(2,2,1); box on; hold on;
    plot(t,temp_mean-ytemp_mean,'b','linewidth',2,'color',clr);
    set(gca,'xscale',xscale_type,'yscale',yscale_type);
    xlim([0 max_t]);
    xlabel('Time (s)');
    ylabel('Delta <T> (K)');
    
    % Species 1
    subplot(2,2,2); box on; hold on;
    plot(t,x1_mean-y1_mean,'b','linewidth',2,'color',clr);
    set(gca,'xscale',xscale_type,'yscale',yscale_type);
    xlim([0 max_t]);
    xlabel('Time (s)');
    ylabel('Delta <z_1>');
    
    % Species 3
    subplot(2,2,3); box on; hold on;
    plot(t,x3_mean-y3_mean,'b','linewidth',2,'color',clr);
    set(gca,'xscale',xscale_type,'yscale',yscale_type);
    xlim([0 max_t]);
    xlabel('Time (s)');
    ylabel('Delta <z_3>');
    
    % Species 6
    subplot(2,2,4); box on; hold on;
    plot(t,x6_mean-y6_mean,'b','linewidth',2,'color',clr);
    set(gca,'xscale',xscale_type,'yscale',yscale_type);
    xlim([0 max_t]);
    xlabel('Time (s)');
    ylabel('Delta <z_6>');
    if max(abs(x3_mean-y3_mean))>x3max
        [x3max,maxloc]=max(abs(x3_mean-y3_mean));
        vals=[x3_mean y3_mean];
    end
end

