% script to plot basic ISAT information for pasr_multi from the files
% isat_1_#.op
%
% specify: np       (number of processes)
%          run_path (location of output files isat_1_#.op)

clear all
np = 32;
run_path = '../test_case_output/';
xscale_type = 'log'; % 'log' or 'linear'
yscale_type = 'log'; % 'log' or 'linear'

figure();
hold on;
box on;

cmap = jet(64);
for i = 1:np
    clr = cmap(1+floor(63*(np-i)/np),:);
    x = load([run_path '/isat_1_' num2str(i-1) '.op']);
    
    qu = x(:,1);
    pr = x(:,2);
    sr = x(:,3);
    gr = x(:,4);
    ad = x(:,5);
    de = x(:,7);
    max_qu = max(qu);
    
    % Primary and secondary retrieves
    subplot(2,2,1); box on; hold on;
    plot(qu,pr+sr,'b','linewidth',2,'color',clr);
    set(gca,'xscale',xscale_type,'yscale',yscale_type);
    xlim([0 max_qu]);
    xlabel('Queries');
    ylabel('Retrieves');
    
    % Grows
    subplot(2,2,2); box on; hold on;
    plot(qu,gr,'b','linewidth',2,'color',clr);
    set(gca,'xscale',xscale_type,'yscale',yscale_type);
    xlim([0 max_qu]);
    xlabel('Queries');
    ylabel('Grows');
    
    % Adds
    subplot(2,2,3); box on; hold on;
    plot(qu,ad,'b','linewidth',2,'color',clr);
    set(gca,'xscale',xscale_type,'yscale',yscale_type);
    xlim([0 max_qu]);
    xlabel('Queries');
    ylabel('Adds');
    
    % Direct evaluations
    subplot(2,2,4); box on; hold on;
    plot(qu,de,'b','linewidth',2,'color',clr);
    set(gca,'xscale',xscale_type,'yscale',yscale_type);
    xlim([0 max_qu]);
    xlabel('Queries');
    ylabel('Direct evaluations');
end
legend([num2str(transpose([1:np]))],-1);
