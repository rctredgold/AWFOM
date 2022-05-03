load('results.mat', 'params')
%oad horns rev
horns_rev = params.farm.turbine_centres;
n_turb = params.farm.n_turb;
density = params.env.density;
diameter = params.turb.diameters;
power_curve = params.turb.power_curve;
location = [0 0 0];



%% pd_8a_hansen
pd_8a_hansen = readmatrix('pwr_deficit_8a_hansen2012.csv');
yaw_angles = zeros(n_turb,1);
wind_direction = 250;

for i = 1:41
[power_8a_hansen(:,i),~] = floris(8,density,wind_direction,...
    horns_rev,yaw_angles,diameter,power_curve,location);
    wind_direction = wind_direction + 1;
end

pwr_deficit = 1- power_8a_hansen(12,:)./power_8a_hansen(11,:);

figure(1)
plot(pd_8a_hansen(:,1), pd_8a_hansen(:,2), 'o','LineWidth',2);
hold on
plot(linspace(-20,20, numel(pwr_deficit)), pwr_deficit, 'LineWidth',2)
xlabel('Normalised Wind Direction (degrees)')
ylabel('Power Deficit [-]')
legend('Observed Data','\texttt{floris.m} performance', 'Interpreter', 'latex')
title({'Performance of \texttt{floris.m} compared to',' experimental data for Horns Rev 1'}, 'Interpreter', 'latex')

%% row_pd_10d_hansen
rowpd_10d_hansen = readmatrix('row_pwr_deficit_10d_hansen2012.csv');
yaw_angles = zeros(n_turb,1);
wind_direction = 312;

[power_10d_hansen,~] = floris(8,density,wind_direction,...
    horns_rev,yaw_angles,diameter,power_curve,location);

row1 = [51 42 33 24 15 6];
row2 = [61 52 43 34 25 16 7];
row3 = [71 62 53 44 35 26 17 8];
row4 = [72 63 54 45 36 27 18 9];
row5 = [73 64 55 46 37 28 19 10];

rd1 = power_10d_hansen(row1);
rd1 = 1-(rd1./rd1(1));

rd2 = power_10d_hansen(row2);
rd2 = 1-(rd2./rd2(1));

rd3 = power_10d_hansen(row3);
rd3 = 1-(rd3./rd3(1));

rd4 = power_10d_hansen(row4);
rd4 = 1-(rd4./rd4(1));

rd5 = power_10d_hansen(row5);
rd5 = 1-(rd5./rd5(1));

figure(2)
plot(rowpd_10d_hansen(:,3), rowpd_10d_hansen(:,4), 'o', 'LineWidth', 1.5);
hold on
plot(rowpd_10d_hansen(:,1), rowpd_10d_hansen(:,2), 'o', 'LineWidth', 1.5);
plot(rowpd_10d_hansen(:,5), rowpd_10d_hansen(:,6), 'o', 'LineWidth', 1.5);
plot(linspace(0,10.4*(numel(rd1)-1), numel(rd1)), rd1, 'LineWidth', 1.5)
plot(linspace(0,10.4*(numel(rd2)-1), numel(rd2)), rd2, 'LineWidth', 1.5)
plot(linspace(0,10.4*(numel(rd3)-1), numel(rd3)), rd3, 'LineWidth', 1.5)
plot(linspace(0,10.4*(numel(rd4)-1), numel(rd4)), rd4, 'LineWidth', 1.5)
plot(linspace(0,10.4*(numel(rd5)-1), numel(rd5)), rd5, 'LineWidth', 1.5)
legend('cL = very stable', 'cL = stable', 'cL = other',...
       'Cloumn 1', 'Row 2', 'Row 3', 'Row 4', 'Row 5', 'NumColumns', 3, 'Color', 'none')
ylim([0 0.6])
xlim([0 70])
xlabel('Downstream distance (m)')
ylabel('Power deficit [-]')
title('Turbine Spacing \(10.4D\), Wind direction 312\(^{\circ}\)', 'Interpreter', 'latex')
grid on


%% Plot and label horns rev
figure (3)
horns_rev = horns_rev./1000;
plot(horns_rev(:,1),horns_rev(:,2), 'o', 'LineWidth',1.5)
hold on
plot(horns_rev(11,1),horns_rev(11,2), 'x', 'MarkerSize',15, 'LineWidth',3)
plot(horns_rev(12,1),horns_rev(12,2), 'x', 'MarkerSize',15,'LineWidth',3)
legend('','wt07', 'wt17')
xlim([-2.5 6])
ylim([-0.5 5 ])
grid on
arrowlength = 0.1;
arrow270x = [0.2 0.2-arrowlength*sind(270)];
arrow270y = [0.5 0.5-arrowlength*cosd(270)];
annotation('textarrow', arrow270x,arrow270y, 'String', '270\(^{\circ}\)', 'Interpreter', 'latex')

arrow312x = [0.2 0.2-arrowlength*sind(312)];
arrow312y = [0.85 0.85-arrowlength*cosd(312)];
annotation('textarrow', arrow312x,arrow312y, 'String', '312\(^{\circ}\)', 'Interpreter', 'latex')

title('Horns Rev 1 Layout')
xlabel('Relative x-coordinates (km)')
ylabel('Relative y-coordinates (km)')
