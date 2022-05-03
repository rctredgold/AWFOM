function [stop, options, optchanged] = patternsearch_outputIterationData(optimvalues, options, flag)

% optimvalues has the following fields
% funccount = evaluations of floris pwr
% fval = power outut of the farm
% iteration = iteration number
% procedure = value
% x is yaw at current iteration


stop = false;
optchanged = false;
switch flag

    case 'init'

        % set output of state values to the current field , and set counter
        % for table row to the solver iteration value
        load("results.mat", "params");
        case_id = params.results.caseID;
        i = optimvalues.iteration;
                
        params.results.(case_id(1)).(case_id(2)).state_vals.iteration(i+1) = i;
        params.results.(case_id(1)).(case_id(2)).state_vals.fval(i+1) = optimvalues.fval;
        params.results.(case_id(1)).(case_id(2)).state_vals.funccount(i+1) = optimvalues.funccount;
        params.results.(case_id(1)).(case_id(2)).state_vals.yaw_angles(:,i+1) = optimvalues.x;
        save("results.mat", "params")

    case 'iter'
        load("results.mat", "params");
        case_id = params.results.caseID;
        i = optimvalues.iteration;
                
        params.results.(case_id(1)).(case_id(2)).state_vals.iteration(i+1) = i;
        params.results.(case_id(1)).(case_id(2)).state_vals.fval(i+1) = optimvalues.fval;
        params.results.(case_id(1)).(case_id(2)).state_vals.funccount(i+1) = optimvalues.funccount;
        params.results.(case_id(1)).(case_id(2)).state_vals.yaw_angles(:,i+1) = optimvalues.x;
        save("results.mat", "params")
    
    case 'interrupt'
    case 'done'
        load("results.mat", "params");
        case_id = params.results.caseID;
        i = optimvalues.iteration;
        
        
        params.results.(case_id(1)).(case_id(2)).state_vals.iteration(i+1) = i;
        params.results.(case_id(1)).(case_id(2)).state_vals.fval(i+1) = optimvalues.fval;
        params.results.(case_id(1)).(case_id(2)).state_vals.funccount(i+1) = optimvalues.funccount;
        params.results.(case_id(1)).(case_id(2)).state_vals.yaw_angles(:,i+1) = optimvalues.x;
        save("results.mat", "params")       
    otherwise
end

end

