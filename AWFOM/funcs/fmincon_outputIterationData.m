function [stop] = fmincon_outputIterationData(x,OptimValues,state)

% optimvalues has the following fields
% funccount = evaluations of floris pwr
% fval = power outut of the farm
% iteration = iteration number
% procedure = value
% x is yaw at current iteration
stop = false;

switch state

    case 'init'

        % set output of state values to the current field , and set counter
        % for table row to the solver iteration value
        load("params.mat", "params");
        case_id = params.results.caseID;
        i = OptimValues.iteration;
                
        params.results.(case_id(1)).(case_id(2)).state_vals.iteration(i+1) = i;
        params.results.(case_id(1)).(case_id(2)).state_vals.fval(i+1) = OptimValues.fval;
        params.results.(case_id(1)).(case_id(2)).state_vals.funccount(i+1) = OptimValues.funccount;
        params.results.(case_id(1)).(case_id(2)).state_vals.yaw_angles(:,i+1) = x;
        save("params.mat", "params")

    case 'iter'
        load("params.mat", "params");
        case_id = params.results.caseID;
        i = OptimValues.iteration;
                
        params.results.(case_id(1)).(case_id(2)).state_vals.iteration(i+1) = i;
        params.results.(case_id(1)).(case_id(2)).state_vals.fval(i+1) = OptimValues.fval;
        params.results.(case_id(1)).(case_id(2)).state_vals.funccount(i+1) = OptimValues.funccount;
        params.results.(case_id(1)).(case_id(2)).state_vals.yaw_angles(:,i+1) = x;
        save("params.mat", "params")

    case 'done'
        load("params.mat", "params");
        case_id = params.results.caseID;
        i = OptimValues.iteration;
        
        
        params.results.(case_id(1)).(case_id(2)).state_vals.iteration(i+1) = i;
        params.results.(case_id(1)).(case_id(2)).state_vals.fval(i+1) = OptimValues.fval;
        params.results.(case_id(1)).(case_id(2)).state_vals.funccount(i+1) = OptimValues.funccount;
        params.results.(case_id(1)).(case_id(2)).state_vals.yaw_angles(:,i+1) = x;
        save("params.mat", "params")       
    otherwise
end

end

