function [power,speed] = floris(wind_speed,density,wind_direction,...
    turbine_centres,yaw_angles,diameter,power_curve,location)
%FLORIS Calculate turbine powers in a wind farm by FLORIS model
%   This function calculates the power of wind turbines within a
%   wind farm, and also the speed at a user-specified point. 
%   The power and speed calculations follows the FLORIS model
%   from Gebraad, P.M.O., et al. Wind Energy, 19, 95-114, 2016.
%   Note that terrain speed-ups are not included so any difference 
%   in height of turbines is assumed to be due to different tower heights.
%   Note that the boundary layer IS NOT included so turbines are assumed to
%   be in uniform flow.
%
%   To call the function:
%
%   [POWER,SPEED]=floris(WIND_SPEED, DENSITY, WIND_DIRECTION, TURBINE_CENTRES, ...
%       ... YAW_ANGLES, DIAMETERS, POWER_CURVE, LOCATION)
%
%   WIND_SPEED: freestream wind speed at the farm at hub height
%   DENSITY: air density at the hub height
%   WIND_DIRECTION: wind direction in degrees (north is zero and clockwise
%       is positive e.g. 90 would be wind blowing from east)
%   TURBINE CENTRES: Nx3 matrix of coordinates of each turbine (so
%       number of rows is number of turbines) where first column is 
%       x-component of turbine centres, second is y and third is z
%   YAW_ANGLES: Nx1 vector of yaw angle "offset" of each turbine i.e. the
%       angle that each turbine is relative to the freestream wind 
%       direction i.e. if yaw is zero then each turbine is facing directly
%       into the wind
%   DIAMETERS: Nx1 vector of diameter of each turbine
%   POWER_CURVE: 2-column matrix containing list of points on the power
%       curve, where first column is speed (m/s) and second is power (W).
%       Linear interpolation is used to determine power after speed at each
%       turbine is known
%   LOCATION: Mx3 matrix containing coordinates of where to evaluate 
%       M values of wind speed
%   POWER (output): Nx1 vector containing power of each turbine
%   SPEED (output): Mx1 vector of wind-speed at points in LOCATION
%
% Version 1.2 released 22/4/21
%
% Code is copyright of Daniel Poole, University of Bristol, 2021.

%% SET-UP

%Check for errors in inputs
err=checkInputs(wind_speed,wind_direction,turbine_centres,yaw_angles,diameter,power_curve,location,density);
if(err)
    return
end

%Constants
[size1, size2]=size(yaw_angles);
nturbs=max(size1,size2);
yaw_rads=yaw_angles*pi/180;
[nlocs,~]=size(location);

%% PREPROCESSING

%Power curve interpolation at current windspeed (assume it's zero outside
%of the power curve values)
power_curve=sortrows(power_curve);
power=interp1(power_curve(:,1),power_curve(:,2),wind_speed,'linear',0);
cp=power/(0.5*density*pi*(diameter(1)/2)^2*wind_speed^3);


%Betz limit check
if(cp>(16/27-1e-10))
    disp('ERROR: Cp is above Betz limit (to within a tolerance)');
    disp(['Current value is ',num2str(cp)]);
    return
end
if(cp<1e-10)
    disp('ERROR: Cp is negative (to within a tolerance)');
    disp(['Current value is ',num2str(cp)]);
    return
end
    
%Solve for induction factor
[hi_dl, psi] = inductionFactorPrep(wind_speed, power_curve, diameter, density);
afactor = inductionFactor(cp,hi_dl,psi);


%Rotate turbines to get such that wind axis is positive x
rotangle=(90+wind_direction)*pi/180;
xt=(turbine_centres(:,1)*cos(rotangle))-(turbine_centres(:,2)*sin(rotangle));
yt=(turbine_centres(:,1)*sin(rotangle))+(turbine_centres(:,2)*cos(rotangle));
zt=turbine_centres(:,3);
xm=(location(:,1)*cos(rotangle))-(location(:,2)*sin(rotangle));
ym=(location(:,1)*sin(rotangle))+(location(:,2)*cos(rotangle));
zm=location(:,3); 

%Determine all upstream turbines (i.e. those with no occlusion)
uset=ones(nturbs,1);
for nt=1:nturbs
    for jj=1:nturbs
        %Only consider if jj turbine is upstream of nt
        if(xt(jj)<xt(nt))
            %Calculate wake centre at nt
            [ywake,zwake]=wakecentre(xt(nt),afactor,yaw_rads(jj),xt(jj),yt(jj),zt(jj),diameter(jj));
            %Calculate size of wake at nt
            dwake(3)=wakeexpand_q(xt(nt),3,diameter(jj),xt(jj));
            %Determine overlap
            lhs=sqrt((ywake-yt(nt))^2+(zwake-zt(nt))^2);
            rhs=0.5*(dwake(3)+diameter(nt));
            if(lhs < rhs)
                uset(nt)=0;
                break
            end
        end
    end
end

%% POWER CALCULATIONS

%Velocity at upstream turbines is freeestream
uturb=wind_speed*uset;

for nt=1:nturbs
    
    %For turbine nt, if it is not a freestream turbine, find velocity
    %Need to find upstream turbine with most overlap
    if(uset(nt) == 0)
        
        olapmax=-1000000000.0;
        waketot=0.0;
        
        %Loop over turbines upstream of this one
        for jj=1:nturbs
            if(xt(jj)<xt(nt))
                
                %Wake size and location of jj at nt
                [ywake,zwake]=wakecentre(xt(nt),afactor,yaw_rads(jj),xt(jj),yt(jj),zt(jj),diameter(jj));
                dwake(1)=wakeexpand_q(xt(nt),1,diameter(jj),xt(jj));
                dwake(2)=wakeexpand_q(xt(nt),2,diameter(jj),xt(jj));
                dwake(3)=wakeexpand_q(xt(nt),3,diameter(jj),xt(jj));
                
                %Wake decay
                vwake(1)=wakedecay_q(xt(nt),1,diameter(jj),yaw_rads(jj),xt(jj));
                vwake(2)=wakedecay_q(xt(nt),2,diameter(jj),yaw_rads(jj),xt(jj));
                vwake(3)=wakedecay_q(xt(nt),3,diameter(jj),yaw_rads(jj),xt(jj));

                %Wake overlap with each zone
                olap1=overlap(yt(nt),zt(nt),ywake,zwake,0.5*diameter(nt),0.5*dwake(1));
                olap2=overlap(yt(nt),zt(nt),ywake,zwake,0.5*diameter(nt),0.5*dwake(2));
                olap3=overlap(yt(nt),zt(nt),ywake,zwake,0.5*diameter(nt),0.5*dwake(3));
                olap2=olap2-olap1;
                olap3=olap3-olap2-olap1;
                olap=olap1+olap2+olap3;
                
                %Sum wake zones
                turbarea=pi*(diameter(nt)/2.0d0)^2;
                sumwake=(vwake(1)*min(olap1/turbarea,1.0))+(vwake(2)*min(olap2/turbarea,1.0))+(vwake(3)*min(olap3/turbarea,1.0));

                %Wake superposition - root sum squares
                waketot=waketot+(afactor*sumwake)^2;
                
                %Update turbine of greatest overlap
                if( (uset(jj)==1) && (jj ~= nt)) % 
                    if(olap>olapmax)
                        olapmax=olap;
                        iolap=jj;
                    end
                end
                
            end
            
        end
        % INSTEAD OF ASSIGNING WIND VELOCITY AS A FUNCTION OF THE MOST
        % UPSTREAM TURBINE, PLACE THE TURBINE IN A SET, UJJ_SET

        % ORDER THE SET BY X(NT) SMALLEST TO LARGEST
        % EVALUATE




        %Turbine velocity
        uturb(nt)=uturb(iolap)*(1.0-2.0*sqrt(waketot));
    end
	
    %Catch negative uturb (sometimes was giving very small complex number)
    if(abs(imag(uturb(nt)))<1e-6)
        uturb(nt)=real(uturb(nt));
    end
    
    %Power (with yaw modification to match SOWFA results)

    pP = 1.88;
    power(nt)=interp1(power_curve(:,1),power_curve(:,2),uturb(nt),'linear',0);
    yaw_loss=cos(yaw_rads(nt))^pP;
    power(nt)=power(nt)*yaw_loss;
    
    %cp=4.0*afactor*(1.0-afactor)^2*0.768*cos(yaw_rads(nt))^1.88;
    %power(nt)=0.5*1.225*pi*(diameter(nt)/2.0)^2*cp*uturb(nt)^3;
    
end

%% VELOCITY CALCULATIONS

for iloc=1:nlocs

    %Velocity deficit at upstream turbines is freeestream
    udeficit=0.0;

    %Determine all turbines that current point is within the wake of
    ii=0;
    for jj=1:nturbs

        %Can only be within wake of those upstream of current point
        if(xt(jj)<xm(iloc))

            %Wake size and location of jj at point
            [ywake,zwake]=wakecentre(xm(iloc),afactor,yaw_rads(jj),xt(jj),yt(jj),zt(jj),diameter(jj));
            dwake(3)=wakeexpand_q(xm(iloc),3,diameter(jj),xt(jj));

            %If point is within wake, get velocity at point
            if(norm([ym(iloc)-ywake,zm(iloc)-zwake])<=0.5*dwake(3))
                ii=ii+1;

                %Wake centre line
                [ywake,zwake]=wakecentre(xm(iloc),afactor,yaw_rads(jj),xt(jj),yt(jj),zt(jj),diameter(jj));
                dwake(1)=wakeexpand_q(xm(iloc),1,diameter(jj),xt(jj));
                dwake(2)=wakeexpand_q(xm(iloc),2,diameter(jj),xt(jj));
                dwake(3)=wakeexpand_q(xm(iloc),3,diameter(jj),xt(jj));

                %Wake decay
                vwake(1)=wakedecay_q(xm(iloc),1,diameter(jj),yaw_rads(jj),xt(jj));
                vwake(2)=wakedecay_q(xm(iloc),2,diameter(jj),yaw_rads(jj),xt(jj));
                vwake(3)=wakedecay_q(xm(iloc),3,diameter(jj),yaw_rads(jj),xt(jj));

                %Wake zones
                rad=sqrt((ym(iloc)-ywake)^2+(zm(iloc)-zwake)^2);
                if(rad <= dwake(1)/2)
                    wakedecay=vwake(1);
                elseif(rad <=  dwake(2)/2)
                    wakedecay=vwake(2);
                elseif(rad <=  dwake(3)/2)
                    wakedecay=vwake(3);
                else
                    wakedecay=0;
                end

                %Wake velocity
                uwake=uturb(jj)*(1.0-(2.0*afactor*wakedecay));
                udeficit=udeficit+(wind_speed-uwake);
            end
        end
    end

    if(ii==0)
        %Freestream wind
        speed(iloc)=wind_speed;
    else
        %Root sum square for deficits
        speed(iloc)=wind_speed-udeficit/ii; 
    end
    
end

end

%% LOCAL FUNCTIONS

function err=checkInputs(wind_speed,wind_direction,turbine_centres,yaw_angles,diameter,power_curve,location,density)

err=false;

if(wind_speed < 0)
    disp('ERROR: wind speed must be positive');
    disp(['Current value is ',num2str(wind_speed)]);
    err=true;
end

if(density < 0)
    disp('ERROR: density must be positive');
    disp(['Current value is ',num2str(density)]);
    err=true;
end

[nturbs,ndim]=size(turbine_centres);
if(ndim ~= 3)
    disp('ERROR: turbine_centres must be matrix with 3 columns');
    disp(['Currently has ',num2str(ndim),' columns']);
    err=true;
end

[nt1,nt2]=size(yaw_angles);
if(nturbs ~= max(nt1,nt2))
    disp('ERROR: number of yaw angles must be same as number of rows in turbine_centres');
    disp(['yaw_angles currently has ',num2str(max(nt1,nt2)),' entries, while turbine_centres has ',num2str(nturbs),' rows']);
    err=true;
end

[nt1,nt2]=size(diameter);
if(nturbs ~= max(nt1,nt2))
    disp('ERROR: number of diameters must be same as number of rows in turbine_centres');
    disp(['diameter currently has ',num2str(max(nt1,nt2)),' entries, while turbine_centres has ',num2str(nturbs),' rows']);
    err=true;
end

[pc1,pc2]=size(power_curve);
if(pc2 ~= 2)
    disp('ERROR: power_curve must be matrix with 2 columns');
    disp(['Currently has ',num2str(pc2),' columns']);
    err=true;
end

[loc1,loc2]=size(location);
if(loc2 ~= 3)
    disp('ERROR: location must be matrix with 3 columns');
    disp(['Currently has ',num2str(loc2),' columns']);
    err=true;
end

end

% function to choose appropriate value for afactor
function [dl, psi] = inductionFactorPrep(wind_speed, power_curve, diameter, density) 
ws = linspace(min(power_curve(:,1)),max(power_curve(:,1)),1000);
power=interp1(power_curve(:,1),power_curve(:,2),ws,'linear',0);
cp = power./(0.5*density*pi*((diameter(1)/2)^2)*(ws.^3));
maxcp = max(cp);
% find wind speed corresponding to max cp to decide if disk loading is high
% or low
max_ws = ws(cp == maxcp);
if wind_speed <= max_ws
    dl = 1;
else
    dl = 2;
end  
psi = (16/27-1e-10)/maxcp;
end


function afactor=inductionFactor(cp,hi_dl,psi)
sols=roots([4 -8 4 psi*-cp]);
sols=sort(sols);

if hi_dl == 1 || hi_dl == 2

    afactor=sols(hi_dl);
else %if there is an error calculating high or low d_l assume optimal axial induction factor
    afactor=0.3333333;
end

if afactor > 0.5
    afactor = 0.5;
end
end

function [ywake,zwake]=wakecentre(x,a,yawrads,xt,yt,zt,d)
kd= 0.15; %0.15;
ad= -3;
%ad = -0.003;
bd= -0.01;
% % use Buhls equation to find ct
% f = 0.9;
% if (a >= 0) && (a <= 0.4)
%     ct= 4.0*a*f*(1.0-a);
% else
%     ct = 0.8889 +  ((4*f -4.4444)*a) + ((5.5556 - 4*f)*(a^2));    
% end
ct= 4.0*a*(1.0-a);
wakeangle=0.5*(cos(yawrads))^2*sin(yawrads)*ct;
yawoff=(2.0*kd*(x-xt)/d)+1.0;
yawoff=wakeangle*(15.0*yawoff^4+wakeangle^2);
denom=30*kd/d;
denom=denom*((2.0*kd*(x-xt)/d)+1.0)^5;
yawoff=yawoff/denom;
yawoff=yawoff-(wakeangle*d*(15.0+wakeangle^2)/(30.0*kd));
%rotoff=(ad*d)+bd*(x-xt);
%rotoff=0;
rotoff=ad+(bd*(x-xt));
ywake=yt+yawoff+rotoff;
zwake=zt;
end

function dwake=wakeexpand_q(x,q,d,xturb)
ke=0.065;
switch(q)
    case(1)
        me= -0.5;
    case(2)
        me= 0.8; %0.22;
    case(3)
        me= 4;
end
dwake=d+2.0*ke*me*(x-xturb);
dwake=max(dwake,0.0);
end

function vwake=wakedecay_q(x,q,d,gam,xturb)
ke=0.065;
switch(q)
    case(1)
        mu= 0.18;
    case(2)
        mu= 1;
    case(3)
        mu= 5;
end
au= 5.15;
bu= 1.7;%1.66;
m=mu/cos(au+(bu*gam));
vwake=(d/(d+2.0*ke*m*(x-xturb)))^2;
end

function overlap_area=overlap(x0,y0,x1,y1,r0,r1)
d = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
rr0 = r0 * r0;
rr1 = r1 * r1;
%Circles do not overlap
if (d > (r1 + r0))
	overlap_area=0.0;
%Circle1 is completely inside circle0
elseif((d <= abs(r0 - r1)) && (r0 >= r1))
	overlap_area=pi * rr1;
%Circle0 is completely inside circle1
elseif((d <= abs(r0 - r1)) && (r0 < r1))
	overlap_area=pi * rr0;
%Circles partially overlap
else
    phi = (acos((rr0 + (d * d) - rr1) / (2.0 * r0 * d))) * 2.0;
    theta = (acos((rr1 + (d * d) - rr0) / (2.0 * r1 * d))) * 2.0;
    theta=real(theta);
    area1 = 0.5 * theta * rr1 - 0.5 * rr1 * sin(theta);
    area2 = 0.5 * phi * rr0 - 0.5 * rr0 * sin(phi);
    overlap_area= area1 + area2;
end
end
