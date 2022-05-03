function locations = random_turbines(X,Y,N_TURB,R_MIN)

%Define wind farm size
x_length = X;
y_length = Y;


%Define number of turbines
n_turb = N_TURB;
% Define counter for how many turbines can actually fit in the given area
n_turb_max = 0;
% Define minimum distance between turbines
r_min = R_MIN;



% Construct loop parameters
k_max = 5000 * n_turb;
k = 1;

% Declare arrays to hold the x and y coordinate values for the turbine
% locations
x = nan(1, n_turb_max);
y = nan(1, n_turb_max);

while n_turb_max < n_turb && k < k_max
  % Pick a random location
  x_guess = x_length*rand();
  y_guess = y_length*rand();
  if n_turb_max == 0
    % First point cannot clash with any other points
    n_turb_max = n_turb_max + 1;
    x(n_turb_max) = x_guess;
    y(n_turb_max) = y_guess;
    continue;
  end
  % Find distances between new guess and stored turbine locations
  distances = sqrt((x-x_guess) .^ 2 + (y - y_guess) .^ 2);
  if min(distances) >= r_min
    % If its far anough away from the other points or its the first guess,
    % store it
    n_turb_max = n_turb_max + 1;
    x(n_turb_max) = x_guess;
    y(n_turb_max) = y_guess;
  end
  % Increment the loop counter.
  k = k + 1;
end
% Remove remaining NaN values
x = x(1:n_turb_max);
y = y(1:n_turb_max);
z = zeros(1,n_turb_max);


locations = [x;y;z]';
end