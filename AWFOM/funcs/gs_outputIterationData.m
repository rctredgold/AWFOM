function [stop] = gs_outputIterationData(optimValues, state)


stop = false;

switch state

    case 'init'

        % set output of state values to the current field , and set counter
        % for table row to the solver iteration value
        load("results.mat", "params");
        case_id = params.results.caseID;
        i = optimValues.localrunindex;
                    
        params.results.(case_id(1)).(case_id(2)).state_vals.iteration(i+1) = i;
        params.results.(case_id(1)).(case_id(2)).state_vals.fval(i+1) = ...
                                                                    -params.results.(case_id(1)).(case_id(2)).power_initial;
        params.results.(case_id(1)).(case_id(2)).state_vals.funccount(i+1) = 0;        
        params.results.(case_id(1)).(case_id(2)).state_vals.yaw_angles(:,i+1) = params.results.(case_id(1)).(case_id(2)).yaw_initial;        
        save("results.mat", "params")

    case 'iter'
        load("results.mat", "params");
        case_id = params.results.caseID;
        i = optimValues.localrunindex;
                
        params.results.(case_id(1)).(case_id(2)).state_vals.iteration(i+1) = i;
        params.results.(case_id(1)).(case_id(2)).state_vals.fval(i+1) = optimValues.bestfval;
        params.results.(case_id(1)).(case_id(2)).state_vals.funccount(i+1) = optimValues.funccount;        
        params.results.(case_id(1)).(case_id(2)).state_vals.yaw_angles(:,i+1) = optimValues.bestx;        
        save("results.mat", "params")        
 
    case 'done'
        load("results.mat", "params");
        case_id = params.results.caseID;
        i = optimValues.localrunindex;        
        
        params.results.(case_id(1)).(case_id(2)).state_vals.iteration(i+1) = i;
        params.results.(case_id(1)).(case_id(2)).state_vals.fval(i+1) = optimValues.bestfval;
        params.results.(case_id(1)).(case_id(2)).state_vals.funccount(i+1) = optimValues.funccount;        
        params.results.(case_id(1)).(case_id(2)).state_vals.yaw_angles(:,i+1) = optimValues.bestx;  
        
        save("results.mat", "params")       
    otherwise
end

end

