function wind_spd_array = flow_field_eval(farm_size, wind_direction, turbine_centres,...
                                            yaw, diameters)

load(".\results.mat",'params');
ff_res = params.farm.ff_res;

%create mesh and plot evaluation of FLORIS for wind speeds at optimised yaw
%positions


%create max values for mesh from existing turbine centres 
%with a buffer around the system of 3 blade diameters
%diameter_max = abs(max(DIAMETER));

x_min = farm_size(1);
y_min = farm_size(2);
x_max = farm_size(3);
y_max = farm_size(4);

%set resolution of flow field to be evaluated


%create x and y values for mesh at locations in a resolution specified by
%ff_res
x = x_min:((x_max-x_min)/(ff_res-1)):x_max;
y = y_min:((y_max-y_min)/(ff_res-1)):y_max;

%create meshes of gridpoints
[X,Y] = meshgrid(x,y);

coords(:,1) = reshape(X,1,numel(X));
coords(:,2) = reshape(Y,1,numel(Y));
coords(:,3) = zeros;
%This returns an array ff_res x ff_res,

[~,spd] = floris(params.env.wind_speed,params.env.density,wind_direction,...
         turbine_centres,yaw, diameters,params.turb.power_curve,coords);


wind_spd_array = reshape(spd,numel(y),numel(x));

end