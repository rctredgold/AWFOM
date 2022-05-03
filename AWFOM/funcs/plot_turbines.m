function tip_locs = plot_turbines(yaw_current, n_turb, wind_direction,...
                                    diameters, turbine_centres)

radii = diameters./2;

%yaw is relative to wind direction therefore to plot yaw relative to axes
%we need to add the wind direction

yaw_relative = yaw_current(:,1) + (wind_direction.*ones(n_turb,1));

tipx_relative = zeros(n_turb,2);
tipy_relative = zeros(n_turb,2);

for i = 1:n_turb
    tipx_relative(i,:) = [-radii(i)*cosd(yaw_relative(i)), radii(i)*cosd(yaw_relative(i))];
    tipy_relative(i,:) = [radii(i)*sind(yaw_relative(i)), -radii(i)*sind(yaw_relative(i))];
end

tipx = [turbine_centres(:,1)+tipx_relative(:,1), turbine_centres(:,1)+tipx_relative(:,2)];
tipy = [turbine_centres(:,2)+tipy_relative(:,1), turbine_centres(:,2)+tipy_relative(:,2)];


tip_locs = [tipx, tipy];
end

