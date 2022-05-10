function success = convert_globcon(filename)

%set succesful exit flag to 0, where it will change to 1 upon a succesful
%matrix creation
success = 0;
%filename = "params.mat";

load(filename,'params')

n_turb = params.farm.n_turb;

yaw_init = params.globcon.yaw_init;
mu = params.globcon.mu;
ub = params.globcon.ub;
lb = params.globcon.lb;
eql_n = params.globcon.eql_n;
rnd_n = params.globcon.rnd_n;
wdsweep_n = params.globcon.wdsweep_n;
farmsz_n = params.globcon.farmsz_n;
turbine_centres = params.globcon.turbine_centres;
diameters = params.turb.diameters;


%choose all turbines facing same way if there has been an error in start
%point generation
if eql_n > 0 && rnd_n > 0
    rnd_n = 0;
end



%check if sigma needs to be changed to a vector for the normal distribution
sigma = params.globcon.sigma;
sigma_max = params.globcon.sigma_max;

%make vector of sigma values
if rnd_n > 1

    if sigma_max > sigma
        sigma = linspace(sigma,sigma_max,rnd_n);

    else
        sigma = sigma.*ones(1,rnd_n);

    end

end

%% if farm size is not changing
if farmsz_n < 1
%equal starting yaw values
if eql_n >= 1   
    
    yaw_init = zeros(n_turb,eql_n);
    yaw_init(1,:) = params.globcon.yaw_init(1,:);
    for i = 1:eql_n
        
        %populate each column with the first value, n_turb times
        single_yaw = yaw_init(1,i);
        yaw_init(:,i) = single_yaw.*ones(n_turb,1);
        yaw_init_cell{i} = yaw_init(:,i);
        diameters_cell{i} = diameters*ones(size(yaw_init,1));
        lb_cell{i} = lb*ones(size(yaw_init,1));
        ub_cell{i} = ub*ones(size(yaw_init,1));

    end
yaw_init = yaw_init_cell;
lb = lb_cell;
ub = ub_cell;
diameters = diameters_cell;
success = 1;

%normally distributed starting values
elseif rnd_n >= 1
    yaw_init = zeros(n_turb,rnd_n);
    for i = 1:rnd_n
            yaw_init(:,i) = normrnd(mu,sigma(i),n_turb,1);
            yaw_init_cell{i} = min(max(yaw_init(:,i),lb),ub);
            diameters_cell{i} = diameters*ones(size(yaw_init,1));
            lb_cell{i} = lb*ones(size(yaw_init,1));
            ub_cell{i} = ub*ones(size(yaw_init,1));
    end
yaw_init = yaw_init_cell;
lb = lb_cell;
ub = ub_cell;
diameters = diameters_cell;
success = 1;

elseif (eql_n + rnd_n) < 1
        if wdsweep_n > 0
        case_it = wdsweep_n;
        yaw_init_cell{1} = zeros(n_turb,case_it);
        diameters_cell{i} = diameters*ones(size(yaw_init,1));
        lb_cell{i} = lb*ones(size(yaw_init,1));
        ub_cell{i} = ub*ones(size(yaw_init,1));
        end


yaw_init = yaw_init_cell;
lb = lb_cell;
ub = ub_cell;
diameters = diameters_cell;

success = 1;
end


%% If farm size is changing
elseif farmsz_n >= 1
    case_it = farmsz_n;
    if (eql_n + rnd_n) < 1

        for i = 1:case_it
            fs = size(turbine_centres{i},1);
            fs_yaw_init{i} = zeros(fs,1);
            fs_ub{i} = ub*ones(fs,1);
            fs_lb{i} = lb*ones(fs,1);
            fs_diam{i} = diameters*ones(fs,1);
        end            
        success = 1;
    
    elseif eql_n >= 1 
        yaw_init(1,:) = params.globcon.yaw_init(1,:);
        for i =1:case_it
        fs = size(turbine_centres{i},1);             
        %populate each farm with the initial yaw value
        single_yaw = yaw_init(1,i);
        fs_yaw_init{i} = single_yaw.*ones(fs,1);
        fs_ub{i} = ub*ones(fs,1);
        fs_lb{i} = lb*ones(fs,1);
        fs_diam{i} = diameters*ones(fs,1);
        end
        success = 1;

    elseif rnd_n >= 1
        
        for i =1:case_it
        fs = size(turbine_centres{i},1);             
        fs_yaw_init{i} = normrnd(mu,sigma(i),fs,1);
        fs_yaw_init{i} = min(max(cell2mat(fs_yaw_init(i)),lb),ub);
        fs_ub{i} = ub*ones(fs,1);
        fs_lb{i} = lb*ones(fs,1);
        fs_diam{i} = diameters*ones(fs,1);
        end
        success = 1;


        
    end
    yaw_init = fs_yaw_init;
    diameters = fs_diam;
    n_turb = size(turbine_centres{1},1);
    ub = fs_ub;
    lb = fs_lb;
else
    success = 0;
    return
end


%% IF CONVERSION IS SUCCESFUL, SET UPPER AND LOWER BOUNDS AND SAVE TO PARAMS.MAT

if success == 1
params.globcon.yaw_init = yaw_init;
params.turb.diameters = diameters;
params.globcon.n_turb = n_turb;
params.globcon.ub = ub;
params.globcon.lb = lb;

save(filename, 'params');
else
    return
end
end