function sortedlocs = sortlocs(locs, wind_direction)

w = deg2rad(wind_direction);

s = size(locs,1);

xmin = min(locs(:,1));
xmax = max(locs(:,1));

ymin = min(locs(:,2));
ymax = max(locs(:,2));

D = max([xmax-xmin, ymax-ymin]);


centre_point_farm = [(xmax-xmin)/2,(ymax-ymin)/2];


%create 2 points defining a line as the wind front approaching the farm
windfront_centre = centre_point_farm + [D*sin(w), D*cos(w)];

windfront1 = windfront_centre + [-D*cos(w), D*sin(w)];
windfront2 = windfront_centre + [D*cos(w), -D*sin(w)];

windfront = [windfront1; windfront2];




%get numerator/denominator for distance from point to line equation
numerators = abs((windfront(2,1)-windfront(1,1)) .* (windfront(1,2)-locs(:,2)) ...
                    - (windfront(1,1)-locs(:,1)) .* (windfront(2,2)-windfront(1,2)));

denominators = ones(s,1).*sqrt((windfront(1,1)-windfront(2,1))^2 +(windfront(1,2)-windfront(2,2))^2);

distances = numerators./denominators;

%sort turbines based on their proximity to the windfront
[~,I] = sort(distances);

sortedlocs = locs(I,:);

end