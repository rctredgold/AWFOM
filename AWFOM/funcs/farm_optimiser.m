function farm_optimiser(solver_selection, tx)
%% VARIABLE SETUP
tx.Value = 'Initialising variables and running first optimisation';
drawnow()
%assign values for running certain solvers
use_fmincon = solver_selection(1);
use_patternsearch = solver_selection(2);
use_particleswarm = solver_selection(3);
use_globalsearch = solver_selection(4);
use_simulannealbnd = solver_selection(5);
use_surrogateopt = solver_selection(6);

%load farm settings and clear results folder
load("params.mat",'params');
if isfield(params,'results') == 1
params = rmfield(params,"results");
save("params.mat", 'params');
end



n_turb = params.farm.n_turb;
wind_speed = params.env.wind_speed;
density = params.env.density;
wind_direction = params.env.wind_direction;
diameter = params.turb.diameters;
turbine_centres = params.farm.turbine_centres;
power_curve = params.turb.power_curve;

%set locations to zero as having multiple locations reduces the floris run
%time and they are irrelevant for optimising turbine yaw
locations = [0 0 0];



%set global constraints from params file
ub = params.globcon.ub;
lb = params.globcon.lb;
x0 = params.globcon.yaw_init;

wdsweep_n = params.globcon.wdsweep_n;
farmsz_n = params.globcon.farmsz_n;

case_iterations = size(x0,2);

%set objective function



%% OPTIMISE USING FMINCON

if use_fmincon == 1
    load(".\optimisation\fmincon_options.mat", 'fmincon_options');
    fmincon_options.OutputFcn = @fmincon_outputIterationData;  
    
    
    for p = 1:case_iterations
    
    % before running optimisation create the results field and add the initial 
    % power value and case identifier to it. the case identifier is so that
    % the iteration ouput function knows where to send the results for the
    % current iteration
    it_n = num2str(p);
    r_n = strcat('r_',it_n);    
    params.results.caseID = string({'fmincon',r_n});
    
    if wdsweep_n > 0
    wind_direction = params.globcon.wd_range(p);
    end

    if farmsz_n > 0
        turbine_centres = cell2mat(params.globcon.turbine_centres(p));
        if params.globcon.sortlocs == 1
            turbine_centres = sortlocs(turbine_centres, wind_direction);
            params.globcon.turbine_centres{p} = turbine_centres;
        end
    else
        if params.globcon.sortlocs == 1
            turbine_centres = sortlocs(turbine_centres, wind_direction);
            params.farm.turbine_centres = turbine_centres;
        end
    end
    
    yaw_angles = cell2mat(params.globcon.yaw_init(p));
    ub = cell2mat(params.globcon.ub(p));
    lb = cell2mat(params.globcon.lb(p));
    diameter = cell2mat(params.turb.diameters(p));

    objective = @(yaw_angles)-floris_pwr(wind_speed,density,wind_direction,...
    turbine_centres,yaw_angles,diameter,power_curve,locations);

    params.results.fmincon.(r_n).yaw_initial = yaw_angles;
    params.results.fmincon.(r_n).pwr_per_turb_initial = floris(wind_speed,density,wind_direction,...
    turbine_centres,yaw_angles,diameter,power_curve,locations)';
    params.results.fmincon.(r_n).power_initial = floris_pwr(wind_speed,density,wind_direction,...
    turbine_centres,yaw_angles,diameter,power_curve,locations);
    
    save("params.mat", "params");
    [yaw_optimal_fmincon, power_optimal, exit_flag, output] = fmincon(objective, yaw_angles,[],[],[],[], lb, ub,[],fmincon_options);
    load("params.mat", "params");
    
    params.results.fmincon.(r_n).yaw_optimal = yaw_optimal_fmincon;
    params.results.fmincon.(r_n).power_optimal = power_optimal;
    params.results.fmincon.(r_n).exit_flag = exit_flag;
    params.results.fmincon.(r_n).output = output;
    params.results.fmincon.(r_n).pwr_per_turb = floris(wind_speed,density,wind_direction,...
    turbine_centres,yaw_optimal_fmincon,diameter,power_curve,locations)';

    %decapitalise funccount for fmincon to be consistent with other optimiser ouptuts 
    params.results.fmincon.(r_n).output.funccount = params.results.fmincon.(r_n).output.funcCount; 
    params.results.fmincon.(r_n).output = rmfield(params.results.fmincon.(r_n).output, 'funcCount');

    params.results.fmincon.(r_n).wind_direction = wind_direction;
    params.results.fmincon.(r_n).turbine_centres = turbine_centres;
    params.results.fmincon.(r_n).diameters = diameter;
    
    save("params.mat", "params");
    

    msg = strcat(params.results.caseID{1},{' '},'reached stopping constraints for case',{' '},...
                 it_n,'/',string(case_iterations),{' '},'in',{' '},...
                 string(output.iterations),{' '},'iterations');
    tx.Value = msg;
    drawnow()
    end

end

%% OPTIMISE USING PATTERNSEARCH

if use_patternsearch == 1
    load(".\optimisation\patternsearch_options.mat", 'patternsearch_options');
    patternsearch_options.OutputFcn = @patternsearch_outputIterationData; 

    for p = 1:case_iterations
    it_n = num2str(p);
    r_n = strcat('r_',it_n);
    params.results.caseID = string({'patternsearch',r_n});
    
    if wdsweep_n > 0
    wind_direction = params.globcon.wd_range(p);
    end

    if farmsz_n > 0
        turbine_centres = cell2mat(params.globcon.turbine_centres(p));
        if params.globcon.sortlocs == 1
            turbine_centres = sortlocs(turbine_centres, wind_direction);
            params.globcon.turbine_centres{p} = turbine_centres;
        end
    else
        if params.globcon.sortlocs == 1
            turbine_centres = sortlocs(turbine_centres, wind_direction);
            params.farm.turbine_centres = turbine_centres;
        end
    end

    yaw_angles = cell2mat(params.globcon.yaw_init(p));
    ub = cell2mat(params.globcon.ub(p));
    lb = cell2mat(params.globcon.lb(p));
    diameter = cell2mat(params.turb.diameters(p));

    objective = @(yaw_angles)-floris_pwr(wind_speed,density,wind_direction,...
    turbine_centres,yaw_angles,diameter,power_curve,locations);

    params.results.patternsearch.(r_n).yaw_initial = yaw_angles;
    params.results.patternsearch.(r_n).pwr_per_turb_initial = floris(wind_speed,density,wind_direction,...
    turbine_centres,yaw_angles,diameter,power_curve,locations)';
    params.results.patternsearch.(r_n).power_initial = floris_pwr(wind_speed,density,wind_direction,...
    turbine_centres,yaw_angles,diameter,power_curve,locations);

    save("params.mat", "params");
    [yaw_optimal_patternsearch, power_optimal, exit_flag, output] = patternsearch(objective, yaw_angles,[],[],[],[], lb, ub,[], patternsearch_options);
    load("params.mat", "params");

    params.results.patternsearch.(r_n).yaw_optimal = yaw_optimal_patternsearch;
    params.results.patternsearch.(r_n).power_optimal = power_optimal;
    params.results.patternsearch.(r_n).exit_flag = exit_flag;
    params.results.patternsearch.(r_n).output = output;
    params.results.patternsearch.(r_n).output = output;
    params.results.patternsearch.(r_n).pwr_per_turb = floris(wind_speed,density,wind_direction,...
    turbine_centres,yaw_optimal_patternsearch,diameter,power_curve,locations)';

    params.results.patternsearch.(r_n).wind_direction = wind_direction;
    params.results.patternsearch.(r_n).turbine_centres = turbine_centres;
    params.results.patternsearch.(r_n).diameters = diameter;

    save("params.mat", "params");


   msg = strcat(params.results.caseID{1},{' '},'reached stopping constraints for case',{' '},...
                 it_n,'/',string(case_iterations),{' '},'in',{' '},...
                 string(output.iterations),{' '},'iterations');
    tx.Value = msg;
    drawnow()
    end

end



%% OPTIMISE USING PARTICLE SWARM

if use_particleswarm == 1
    load(".\optimisation\particleswarm_options.mat", 'particleswarm_options');
    particleswarm_options.OutputFcn = @particleswarm_outputIterationData; 

    %if the user has selected a hybrid function then add it to the ps
    %options using the options the user has already set for the hybrid
    %function
    if  isempty(particleswarm_options.HybridFcn) == 0            
        load(".\optimisation\patternsearch_options.mat", 'patternsearch_options');
        load(".\optimisation\fmincon_options.mat", 'fmincon_options');

        particleswarm_options.HybridFcn = {strcat(particleswarm_options.HybridFcn, '_options')};
    end
   

    for p = 1:case_iterations
    it_n = num2str(p);
    r_n = strcat('r_',it_n);
    params.results.caseID = string({'particleswarm',r_n});
    
    if wdsweep_n > 0
    wind_direction = params.globcon.wd_range(p);
    end

    if farmsz_n > 0
        turbine_centres = cell2mat(params.globcon.turbine_centres(p));
        if params.globcon.sortlocs == 1
            turbine_centres = sortlocs(turbine_centres, wind_direction);
            params.globcon.turbine_centres{p} = turbine_centres;
        end
    else
        if params.globcon.sortlocs == 1
            turbine_centres = sortlocs(turbine_centres, wind_direction);
            params.farm.turbine_centres = turbine_centres;
        end
    end

    yaw_angles = cell2mat(params.globcon.yaw_init(p));
    ub = cell2mat(params.globcon.ub(p));
    lb = cell2mat(params.globcon.lb(p));
    diameter = cell2mat(params.turb.diameters(p));

    particleswarm_options.MaxIterations = params.globcon.psMaxIter*size(yaw_angles,1); 

    objective = @(yaw_angles)-floris_pwr(wind_speed,density,wind_direction,...
    turbine_centres,yaw_angles,diameter,power_curve,locations);

    
    save("params.mat", "params");
    [yaw_optimal_particleswarm, power_optimal, exit_flag, output] = particleswarm(objective, size(yaw_angles,1), lb, ub, particleswarm_options);
    load("params.mat", "params");

    params.results.particleswarm.(r_n).yaw_optimal = yaw_optimal_particleswarm';
    params.results.particleswarm.(r_n).power_optimal = power_optimal;
    params.results.particleswarm.(r_n).exit_flag = exit_flag;
    params.results.particleswarm.(r_n).output = output;
    params.results.particleswarm.(r_n).pwr_per_turb = floris(wind_speed,density,wind_direction,...
    turbine_centres,yaw_optimal_particleswarm,diameter,power_curve,locations)';

    params.results.particleswarm.(r_n).wind_direction = wind_direction;
    params.results.particleswarm.(r_n).turbine_centres = turbine_centres;
    params.results.particleswarm.(r_n).diameters = diameter;

    save("params.mat", "params");

    msg = strcat(params.results.caseID{1},{' '},'reached stopping constraints for case',{' '},...
                 it_n,'/',string(case_iterations),{' '},'in',{' '},...
                 string(output.iterations),{' '},'iterations');
    tx.Value = msg;
    drawnow()
    end

end


%% OPTIMISE USING GLOBAL SEARCHparticleswarm
if use_globalsearch == 1
    load(".\optimisation\fmincon_options.mat", 'fmincon_options');
    load(".\optimisation\gs_options.mat", 'gs_options');      
    fmincon_options.OutputFcn = @gs_fmincon_outputIterationData;
     

    for p = 1:case_iterations
    it_n = num2str(p);
    r_n = strcat('r_',it_n);
    params.results.caseID = string({'globalsearch',r_n});
    
    params.results.globalsearch.(r_n).yaw_initial = yaw_angles;
    params.results.globalsearch.(r_n).pwr_per_turb_initial = floris(wind_speed,density,wind_direction,...
    turbine_centres,yaw_angles,diameter,power_curve,locations)';
    params.results.globalsearch.(r_n).power_initial = floris_pwr(wind_speed,density,wind_direction,...
    turbine_centres,yaw_angles,diameter,power_curve,locations);
    
    %create options for globalsearch local solver and start points
    gslocal = createOptimProblem('fmincon','x0',yaw_angles,...
                                 'objective',objective,'lb',lb,'ub',ub,...
                                 'options',fmincon_options);
    gs_options.OutputFcn = @gs_outputIterationData;
    gs = GlobalSearch(gs_options);    

    save("params.mat", "params");
    [yaw_optimal_globalsearch, power_optimal] = run(gs, gslocal);
    load("params.mat", "params");

    params.results.globalsearch.(r_n).yaw_optimal = yaw_optimal_globalsearch';
    params.results.globalsearch.(r_n).power_optimal = power_optimal;

    output.iterations = params.results.globalsearch.(r_n).state_vals.iteration(end);
    output.funccount = params.results.globalsearch.(r_n).state_vals.funccount(end);
        params.results.globalsearch.(r_n).output = output;

    params.results.globalsearch.(r_n).pwr_per_turb = floris(wind_speed,density,wind_direction,...
    turbine_centres,yaw_optimal_globalsearch,diameter,power_curve,locations)';
    save("params.mat", "params");

    msg = strcat(params.results.caseID{1},{' '},'reached stopping constraints for case',{' '},...
                 it_n,'/',string(case_iterations),{' '},'in',{' '},...
                 string(output.iterations),{' '},'iterations');
    tx.Value = msg;
    drawnow()
    end

end


%delete the caseID counter

if isfield(params.results, 'caseID') == 1
    params.results = rmfield(params.results, "caseID");
end

if isfield(params.results, 'wind_speed_array') == 1
    params.results = rmfield(params.results, 'wind_speed_array');
end

save("params.mat", "params");

pause(2)
tx.Value = 'Optimisation Complete';
drawnow()

end