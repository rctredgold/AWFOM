clear

format_results('hornsrev_all_wd_patternsearchandfmincon.mat')
load("formatted_results.mat", "axesData")
for j = 1:2
for i = 1:37

    wd(i,j) = axesData.dataArray(j,i,16);
    power_optim(i,j) = axesData.dataArray(j,i,5);
    power_init(i,j) =  axesData.dataArray(j,i,4);
    iter_tot(i,j) = axesData.dataArray(j,i,7);
    iter_count(i,j) = axesData.dataArray(j,i,13);
    farm_power_iter(i,j) = axesData.dataArray(j,i,14);
    
end
end
   wd = cell2mat(wd);
   power_optim = cell2mat(power_optim);
   power_init = cell2mat(power_init);
   iter_tot = cell2mat(iter_tot);
%    iter_count_ = cell2mat(iter_count);
%    farm_power_iter = cell2mat(farm_power_iter);

   power_increase = power_optim./power_init;
   
   wd = wd.*(pi/180);

   figure(1)
   polarplot(wd(:,1), power_init(:,1));
   hold on
   polarplot(wd(:,1), power_optim(:,1),'-.', LineWidth=2);
   polarplot(wd(:,2), power_optim(:,2),'--', LineWidth=2);
   polarplot(wd(:,1), iter_tot(:,1)./10, 'o');
   polarplot(wd(:,2), iter_tot(:,2)./10, 'o');
   legend('Initial Farm Power','Optimal Farm Power returned by fmincon', ...
          'Optimal Farm Power returned by patternsearch', ...
          'Normalised Solver Iterations of fmincon',...
          'Normalised Solver Iterations of patternsearch');
   title({'Polar Plot showing the total Power Output of ',...
          'Horns Rev 1 over all wind directions'});
   ax = gca;
   ax.ThetaZeroLocation = 'top';
   ax.ThetaDir = 'clockwise';
  

    figure(2)
    plot(cell2mat(iter_count(27,1)), cell2mat(farm_power_iter(27,1)), 'r')
    hold on
    plot(cell2mat(iter_count(27,2)), cell2mat(farm_power_iter(27,2)), 'b')
    plot(cell2mat(iter_count(15,1)), cell2mat(farm_power_iter(15,1)), 'r-.')
    plot(cell2mat(iter_count(15,2)), cell2mat(farm_power_iter(15,2)), 'b-.')
    plot(cell2mat(iter_count(1,1)), cell2mat(farm_power_iter(1,1)), 'r--')
    plot(cell2mat(iter_count(1,2)), cell2mat(farm_power_iter(1,2)), 'b--')

    legend('fmincon \(260^{\circ}\)','patternsearch \(260^{\circ}\)','fmincon \(140^{\circ}\)','patternsearch \(140^{\circ}\)',...
            'fmincon \(0^{\circ}\)','patternsearch \(0^{\circ}\)', 'Interpreter', 'latex')
    grid on
    xlim([0 100])
    xlabel('Solver Iterations')
    ylabel('Farm Power (MW)')
   