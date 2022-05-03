
power_curve = readmatrix('V80_powercurve.csv');
power_curve=sortrows(power_curve);
ws = linspace(min(power_curve(:,1)),max(power_curve(:,1)),1000);
power=interp1(power_curve(:,1),power_curve(:,2),ws,'linear',0);
cp = power./(0.5*1.225*pi*((80/2)^2)*(ws.^3));

plot(ws,cp)