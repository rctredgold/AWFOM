function format_results(filename)

load(filename, 'params');

% Delete caseID counter from results field so that it isn't read as a
% solver name
if isfield(params.results, 'caseID') == 1
    params.results = rmfield(params.results, "caseID");
end

if isfield(params.results, 'wind_speed_array') == 1
    params.results = rmfield(params.results, 'wind_speed_array');
end


%set number of cases to iterate over
case_iterations = max([params.globcon.rnd_n, params.globcon.eql_n, params.globcon.farmsz_n,...
                      params.globcon.wdsweep_n]);

%% INITIALISE UNIVERSAL PLOT DATA 
axesData.solverNames = fieldnames(params.results);


%% CREATE ARRAY WITH ALL IMPORTANT DATA FOR THE PLOTTING APP

%preallocate cell array

solver_n = length(axesData.solverNames);
dataArray = cell(solver_n, case_iterations, 18);
%fill cell array with applicable plot data, loop through cases for each
%solver
for j = 1:solver_n

    s_n = char(axesData.solverNames(j));

    for i = 1:case_iterations
    
        it_n = num2str(i);
        x = strcat('r_',it_n); 
        % set case names
        axesData.caseNames{i,1} = x;
        
        % entry 1 in array is initial yaw settings
        dataArray{j,i,1} = params.results.(s_n).(x).yaw_initial;
        % entry 2 in array is optimal yaw settings
        dataArray{j,i,2} = params.results.(s_n).(x).yaw_optimal;
        % entry 3 in array is delta yaw
        dataArray{j,i,3} = dataArray{j,i,2} - dataArray{j,i,1};
        

        % POWER IN MW FOR WHOLE FARM
        % entry 4 in array is initial farm power
        dataArray{j,i,4} = params.results.(s_n).(x).power_initial/1000000;
        % entry 5 in array is optimised farm power
        dataArray{j,i,5} = -params.results.(s_n).(x).power_optimal/1000000;
        % entry 6 in array is delta farm power
        dataArray{j,i,6} = dataArray{j,i,5} - dataArray{j,i,4};

        % entry 7 in array is total solver iterations
        dataArray{j,i,7} = params.results.(s_n).(x).output.iterations;
        
        % POWER IN MW FOR INDIVIDUAL TURBINES
        % entry 8 in array is initial power per turbine
        dataArray{j,i,8} = params.results.(s_n).(x).pwr_per_turb_initial/1000000;
        % entry 9 in array is optimised power per turbine
        dataArray{j,i,9} = params.results.(s_n).(x).pwr_per_turb/1000000;
        % entry 10 in array is delta power per turbine
        dataArray{j,i,10} = dataArray{j,i,9} - dataArray{j,i,8};

        % entry 11 in array is total function evaluations
        dataArray{j,i,11} = params.results.(s_n).(x).output.funccount;
        

        % FOR EACH ITERATION
        % entry 12 in array is the function value at each iteration
        dataArray{j,i,12} = params.results.(s_n).(x).state_vals.yaw_angles;
        % entry 13 in array is the iteration count
        dataArray{j,i,13} = params.results.(s_n).(x).state_vals.iteration;        
        % entry 14 in array is the farm power in MW at each iteration
        dataArray{j,i,14} = -params.results.(s_n).(x).state_vals.fval/1000000;
        % entry 15 in array is the function count at each iteration
        dataArray{j,i,15} = params.results.(s_n).(x).state_vals.funccount;
        
        % entry 16 in array is the wind_direction
        dataArray{j,i,16} = params.results.(s_n).(x).wind_direction;
        % entry 17 in array is the turbine_centres
        dataArray{j,i,17} = params.results.(s_n).(x).turbine_centres;
        % entry 18 in array is the diameters
        dataArray{j,i,18} = params.results.(s_n).(x).diameters;

                        
    end
end
axesData.dataArray = dataArray;
save('formatted_results.mat', 'axesData');


%% FUNCTIONS

% Function to assign turbine positions to polar coordinates
function polarTurbPositions = turbTheta(n)
    turb_theta = (2*pi/n).*(0:n-1);
    %turb_theta(1:2:end) = -turb_theta(1:2:end);
    polarTurbPositions = turb_theta;
end
end