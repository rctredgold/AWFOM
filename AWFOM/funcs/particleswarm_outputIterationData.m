function [stop] = particleswarm_outputIterationData(OptimValues, state)

% optimvalues has the following fields
% funccount = evaluations of floris pwr
% bestfval = power outut of the farm
% iteration = iteration number
% bestx = yaw_angles
% swarm/swarmfvals is x and fval for each particle at current iteration
stop = false;

switch state

    case 'init'

        % set output of state values to the current field , and set counter
        % for table row to the solver iteration value
        load("results.mat", "params");
        wind_speed = params.env.wind_speed;
        density = params.env.density;
        wind_direction = params.env.wind_direction;
        diameter = params.turb.diameters;
        turbine_centres = params.farm.turbine_centres;
        power_curve = params.turb.power_curve;
        case_id = params.results.caseID;

        locations = [0 0 0];

        i = OptimValues.iteration;

        params.results.(case_id(1)).(case_id(2)).yaw_initial = OptimValues.bestx;
        params.results.(case_id(1)).(case_id(2)).pwr_per_turb_initial = floris(wind_speed,density,wind_direction,...
                     turbine_centres,OptimValues.bestx,diameter,power_curve,locations)';
        params.results.(case_id(1)).(case_id(2)).power_initial = floris_pwr(wind_speed,density,wind_direction,...
                      turbine_centres,OptimValues.bestx,diameter,power_curve,locations);
                
        params.results.(case_id(1)).(case_id(2)).state_vals.iteration(i+1) = i;
        params.results.(case_id(1)).(case_id(2)).state_vals.fval(i+1) = OptimValues.bestfval;
        params.results.(case_id(1)).(case_id(2)).state_vals.funccount(i+1) = OptimValues.funccount;
        params.results.(case_id(1)).(case_id(2)).state_vals.swarm(:,:,i+1) = OptimValues.swarm;
        params.results.(case_id(1)).(case_id(2)).state_vals.yaw_angles(:,i+1) = OptimValues.bestx;
        params.results.(case_id(1)).(case_id(2)).state_vals.swarmfvals(:,i+1) = OptimValues.swarmfvals;
        save("results.mat", "params")

    case 'iter'
        load("results.mat", "params");
        case_id = params.results.caseID;
        i = OptimValues.iteration;
                
        params.results.(case_id(1)).(case_id(2)).state_vals.iteration(i+1) = i;
        params.results.(case_id(1)).(case_id(2)).state_vals.fval(i+1) = OptimValues.bestfval;
        params.results.(case_id(1)).(case_id(2)).state_vals.funccount(i+1) = OptimValues.funccount;
        params.results.(case_id(1)).(case_id(2)).state_vals.swarm(:,:,i+1) = OptimValues.swarm;
        params.results.(case_id(1)).(case_id(2)).state_vals.yaw_angles(:,i+1) = OptimValues.bestx;
        params.results.(case_id(1)).(case_id(2)).state_vals.swarmfvals(:,i+1) = OptimValues.swarmfvals;
        save("results.mat", "params")

    case 'done'
        load("results.mat", "params");
        case_id = params.results.caseID;
        i = OptimValues.iteration;
        
        
        params.results.(case_id(1)).(case_id(2)).state_vals.iteration(i+1) = i;
        params.results.(case_id(1)).(case_id(2)).state_vals.fval(i+1) = OptimValues.bestfval;
        params.results.(case_id(1)).(case_id(2)).state_vals.funccount(i+1) = OptimValues.funccount;
        params.results.(case_id(1)).(case_id(2)).state_vals.swarm(:,:,i+1) = OptimValues.swarm;
        params.results.(case_id(1)).(case_id(2)).state_vals.yaw_angles(:,i+1) = OptimValues.bestx;
        params.results.(case_id(1)).(case_id(2)).state_vals.swarmfvals(:,i+1) = OptimValues.swarmfvals;
        save("results.mat", "params")       
    otherwise
end

end

