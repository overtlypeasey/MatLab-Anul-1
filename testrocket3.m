function twoStageGravityTurnLEO_withParabola()
%--------------------------------------------------------------------------
% Two-Stage 2D Gravity Turn Simulation to LEO with Predefined Parabolic Trajectory
%
% - Implements a two-stage rocket ascent with gravity turn
% - Uses 1976 Standard Atmosphere for atmospheric properties
% - Applies constant thrust for each stage
% - Retains reserve fuel in Stage 1 for landing maneuvers
% - Removes Stage 1 mass correctly at separation
% - Prevents the rocket from falling back to the ground by simulating a solid ground
% - Ensures payload remains in stable LEO at 200 km altitude by following a predefined parabola
% - Implements dynamic pitch control to follow the parabolic trajectory
% - Plots Altitude, Velocity, and Heat Flux vs Time along with the predefined parabola
%
% DISCLAIMER:
%   This is a simplified simulation and does not account for many real-world factors
%   such as Earth rotation, 3D flight paths, detailed engine performance, etc.
%--------------------------------------------------------------------------

    %------------------------------
    % EARTH & SIM PARAMETERS
    %------------------------------
    Re          = 6371000;            % [m] Mean Earth radius
    mu          = 3.986004418e14;     % [m^3/s^2] Earth gravitational parameter (G*M)
    tMax        = 2000;               % [s] Maximum simulation time
    dt          = 0.1;                % [s] Time step (reduced for accuracy)
    stage1Sep   = 150;                % [s] Time of Stage 1 separation
    targetAlt   = 200e3;              % [m] Target LEO altitude
    flattenAlt  = 120e3;              % [m] Altitude to begin thrust flattening

    %------------------------------
    % STAGE PARAMETERS
    %------------------------------
    % Stage 1
    thrust1        = 7.6e6;       % [N] Constant thrust for Stage 1 (~7.6 MN)
    mdot1          = 2200;        % [kg/s] Mass flow rate Stage 1
    dryMass1       = 22000;       % [kg] 22-ton dry mass
    propMass1      = 330000;      % [kg] Stage 1 propellant
    stage1Reserve  = 20000;       % [kg] Reserve propellant in Stage 1

    % Stage 2
    thrust2        = 1.2e6;       % [N] Constant thrust for Stage 2 (~1.2 MN)
    mdot2          = 150;         % [kg/s] Mass flow rate Stage 2
    dryMass2       = 5000;        % [kg] 5-ton dry mass
    propMass2      = 90700;       % [kg] Stage 2 propellant (~90.7 t)

    % Payload
    payloadMass    = 10000;       % [kg] Payload mass

    % Combined initial mass (Stage1 + Stage2 + Payload)
    m0 = (dryMass1 + propMass1) + (dryMass2 + propMass2) + payloadMass;

    % Rocket geometry
    rocket_d  = 3.7;                         % [m] Diameter
    crossArea = pi*(rocket_d/2)^2;           % [m^2] Cross-sectional area

    %------------------------------
    % PREDEFINED PARABOLA PARAMETERS
    %------------------------------
    % Define a parabola in Cartesian coordinates: y = a*x^2 + c
    % We choose c = Re to start at Earth's surface
    % At x = x_max, y = Re + targetAlt
    x_max = 200e3;                              % [m] Maximum horizontal distance for parabola
    y_max = Re + targetAlt;                     % [m] Maximum altitude at x_max
    a = (y_max - Re) / x_max^2;                 % [1/m] Parabola coefficient
    c = Re;                                     % [m] y-intercept

    %------------------------------
    % STATE VARIABLES
    %  (r, theta) in polar coordinates
    %  (vr, vth)  Radial and tangential velocities
    %------------------------------
    t   = 0;         % [s] Initial time
    r   = Re;        % [m] Initial radius (Earth surface)
    th  = 0;         % [rad] Initial angle
    vr  = 0;         % [m/s] Initial radial velocity
    vth = 0;         % [m/s] Initial tangential velocity
    m   = m0;        % [kg] Initial mass

    % Track propellant remaining
    stage1PropRem = propMass1;
    stage2PropRem = propMass2;

    % Orbit achieved flag
    orbitAchieved = false;

    % Ground flag
    onGround = true;

    % Data storage: [time, alt, vr, vth, speed, Mach, heatFlux, stageID, x, y, desiredPitch_deg]
    simData = [];

    %------------------------------
    % MAIN SIMULATION LOOP
    %------------------------------
    while (t < tMax) && (~orbitAchieved)
        % Calculate altitude
        alt = r - Re;
        if alt < 0
            fprintf('Rocket has crashed or is below Earth surface!\n');
            break;
        end

        % 1) Check if orbit is achieved
        [inOrbit, ~] = checkOrbitInsertion(r, vr, vth, mu, targetAlt);
        if inOrbit
            orbitAchieved = true;
            fprintf('*** Orbit Achieved at t=%.1f s ***\n', t);
            break;
        end

        % 2) Determine active stage and get current thrust
        if t < stage1Sep
            stageID = 1;
            if stage1PropRem > 0
                thrustCurrent = thrust1;
                dmCurrent = mdot1;
            else
                thrustCurrent = 0;
                dmCurrent = 0;
            end
        else
            stageID = 2;
            if stage2PropRem > 0
                thrustCurrent = thrust2;
                dmCurrent = mdot2;
            else
                thrustCurrent = 0;
                dmCurrent = 0;
            end
        end

        % 3) Calculate local gravity
        gLocal = mu / (r^2);

        % 4) Get atmospheric properties from 1976 Standard Atmosphere
        [T, aLoc, ~, rho] = standardAtmos1976(alt);

        % 5) Calculate speed and Mach number
        speed = sqrt(vr^2 + vth^2);
        if aLoc > 1e-6
            Mach = speed / aLoc;
        else
            Mach = 0;
        end

        % 6) Get drag coefficient based on Mach number
        Cd = getDragCoefficient(Mach);

        % 7) Calculate drag force
        F_drag = 0.5 * Cd * crossArea * rho * speed^2;

        % 8) Convert current polar coordinates to Cartesian for parabola
        x = r * cos(th);
        y = r * sin(th);

        % 9) Determine desired pitch angle based on predefined parabola
        % Calculate dy/dx at current x for the parabola y = a*x^2 + c
        dy_dx = 2 * a * x;
        desiredPitch_deg = atan(dy_dx) * (180/pi);  % Convert to degrees
        desiredPitch_rad = deg2rad(desiredPitch_deg); % Convert to radians

        % Ensure desiredPitch does not exceed 90 degrees
        if desiredPitch_deg > 90
            desiredPitch_deg = 90;
            desiredPitch_rad = deg2rad(90);
        end

        % 10) Decompose thrust into radial and tangential components based on desired pitch angle
        thrustRad = thrustCurrent * cos(desiredPitch_rad);
        thrustTan = thrustCurrent * sin(desiredPitch_rad);

        % 11) Decompose drag into radial and tangential components
        if speed < 1e-3
            dragRad = 0;
            dragTan = 0;
        else
            dragRad = -F_drag * (vr / speed);
            dragTan = -F_drag * (vth / speed);
        end

        % 12) Calculate net radial and tangential forces
        FR = thrustRad + dragRad - m * gLocal;    % Radial force
        FT = thrustTan + dragTan;                % Tangential force

        % 13) Calculate accelerations in polar coordinates
        ar = (FR / m) - (vth^2 / r);             % Radial acceleration
        at = (FT / m) + ((vr * vth) / r);        % Tangential acceleration

        % 14) Calculate heat flux (Sutton-Graves approximation)
        if (speed > 0) && (rho > 1e-10) && (~onGround)
            kSG    = 1.83;        % [W/(m^2*(kg/m^3)^(1/2)*(m/s)^3.05)]
            R_nose = 1.0;         % [m] Nose radius
            heatFlux = kSG * sqrt(rho / R_nose) * (speed^3.05);
        else
            heatFlux = 0;
        end

        % 15) Update velocities and positions using Euler integration
        vr_new  = vr  + ar * dt;
        vth_new = vth + at * dt;
        r_new   = r   + vr_new * dt;
        th_new  = th  + (vth_new / r_new) * dt;

        % 16) Update mass based on propellant consumption
        if stageID == 1
            % Do not burn below Stage 1 reserve
            usableProp = max(0, stage1PropRem - stage1Reserve);
            if usableProp > (dmCurrent * dt)
                stage1PropRem = stage1PropRem - dmCurrent * dt;
                m_new = m - dmCurrent * dt;
            else
                % Burn only the usable propellant
                dmUsed = usableProp;
                stage1PropRem = stage1PropRem - dmUsed;
                m_new = m - dmUsed;
            end
        else
            % Stage 2 burns normally until propellant is exhausted
            if stage2PropRem > (dmCurrent * dt)
                stage2PropRem = stage2PropRem - dmCurrent * dt;
                m_new = m - dmCurrent * dt;
            else
                % Burn only remaining propellant
                dmUsed = stage2PropRem;
                stage2PropRem = 0;
                m_new = m - dmUsed;
            end
        end

        % 17) Stage separation at stage1Sep (t = 150 s)
        if (t < stage1Sep) && (t + dt >= stage1Sep)
            % Remove Stage 1 dry mass and reserve propellant
            leftoverStage1 = dryMass1 + stage1Reserve;  % all Stage1 mass left
            m_new = m_new - leftoverStage1;
            fprintf('Stage 1 separated at t=%.1f s, Altitude=%.1f km\n', t, alt/1e3);
        end

        % 18) Update state variables for next iteration
        vr  = vr_new;
        vth = vth_new;
        r   = r_new;
        th  = th_new;
        m   = max(m_new, 0);  % Prevent negative mass
        t   = t + dt;

        % 19) Store simulation data
        simData(end+1,:) = [ t, alt, vr, vth, speed, Mach, heatFlux, stageID, x, y, desiredPitch_deg ];
    end

    % Final data point
    simData(end+1,:) = [ t, (r_new - Re), vr_new, vth_new, sqrt(vr_new^2 + vth_new^2), ...
                         Mach, heatFlux, stageID, x, y, desiredPitch_deg ];

    % If orbit was not achieved within the loop
    if ~orbitAchieved
        [inOrbit, ~] = checkOrbitInsertion(r_new, vr_new, vth_new, mu, targetAlt);
        if inOrbit
            fprintf('*** Orbit Achieved at final step, t=%.1f s ***\n', t);
        else
            fprintf('End of simulation reached without achieving orbit.\n');
        end
    end

    %------------------------------
    % OUTPUT & PLOTS
    %------------------------------
    % Extract data for plotting
    timeVec        = simData(:,1);
    altVec         = simData(:,2);
    vrVec          = simData(:,3);
    vthVec         = simData(:,4);
    speedVec       = simData(:,5);
    MachVec        = simData(:,6);
    heatFluxVec    = simData(:,7);
    stageIDVec     = simData(:,8);
    xVec           = simData(:,9);
    yVec           = simData(:,10);
    desiredPitchVec = simData(:,11);

    % Define predefined parabola for plotting
    xDesired = linspace(0, max(xVec)*1.1, 1000);    % Extend a bit beyond max x
    yDesired = a * xDesired.^2 + c;

    % Display a summary table (sampled every 10%)
    fprintf('\n   t(s)  Alt(km)   vr(m/s)  vth(m/s)  Speed(m/s)  Mach  HeatFlux(W/m^2)  Stage  x(m)    y(m)      Pitch(deg)\n');
    fprintf('-------------------------------------------------------------------------------------------------------------\n');
    step = floor(length(timeVec)/10);
    for i = 1 : step : length(timeVec)
        fprintf('%6.1f  %7.2f   %8.1f  %8.1f  %10.1f   %4.2f      %8.1f      %d  %7.1f  %7.1f    %10.1f\n', ...
                timeVec(i), altVec(i)/1e3, vrVec(i), vthVec(i), speedVec(i), MachVec(i), heatFluxVec(i), stageIDVec(i), xVec(i), yVec(i), desiredPitchVec(i));
    end

    % Plotting
    figure('Name','Payload Ascent to LEO with Predefined Parabolic Trajectory','Color','w','Position',[100 100 1200 800]);

    % Altitude vs Time
    subplot(3,1,1);
    plot(timeVec, altVec/1e3, 'b-', 'LineWidth', 2); hold on;
    xlabel('Time (s)');
    ylabel('Altitude (km)');
    title('Altitude vs. Time');
    grid on;
    legend('Actual Trajectory');

    % Velocity vs Time
    subplot(3,1,2);
    plot(timeVec, speedVec, 'r-', 'LineWidth', 2); hold on;
    xlabel('Time (s)');
    ylabel('Speed (m/s)');
    title('Speed vs. Time');
    grid on;
    legend('Actual Speed');

    % Heat Flux vs Time
    subplot(3,1,3);
    plot(timeVec, heatFluxVec, 'm-', 'LineWidth', 2); hold on;
    xlabel('Time (s)');
    ylabel('Heat Flux (W/m^2)');
    title('Heat Flux vs. Time (Sutton-Graves Approx)');
    grid on;
    legend('Heat Flux');

    % Plot the predefined parabola
    figure('Name','Rocket Trajectory vs. Predefined Parabola','Color','w','Position',[150 150 800 600]);
    plot(xVec, yVec, 'b-', 'LineWidth', 2); hold on;
    plot(xDesired, yDesired, 'r--', 'LineWidth', 2);
    xlabel('Horizontal Distance (m)');
    ylabel('Vertical Distance (m)');
    title('Rocket Trajectory vs. Predefined Parabola');
    legend('Actual Trajectory', 'Predefined Parabola');
    grid on;

end

%% HELPER FUNCTIONS

%--------------------------------------------------------------------------
% Dynamic Pitch Control Function
% Determines the desired pitch angle based on current position to follow the predefined parabola
%--------------------------------------------------------------------------
% (This function is now integrated within the main loop as 'controlPitchAngle')

%--------------------------------------------------------------------------
% Check if orbit insertion is achieved
%--------------------------------------------------------------------------
function [inOrbit, info] = checkOrbitInsertion(r, vr, vth, mu, targetAlt)
    % Define tolerances
    Re     = 6371000;             % [m] Earth radius
    alt    = r - Re;
    vCirc  = sqrt(mu / r);        % Circular orbit velocity
    tolAlt = 10e3;                % [m] Tolerance for altitude
    tolVr  = 50;                  % [m/s] Tolerance for radial velocity
    tolVth = 100;                 % [m/s] Tolerance for tangential velocity

    if (abs(alt - targetAlt) < tolAlt) && (abs(vr) < tolVr) && (abs(vth - vCirc) < tolVth)
        inOrbit = true;
        info    = 'Orbit Achieved';
    else
        inOrbit = false;
        info    = 'Not in Orbit';
    end
end

%--------------------------------------------------------------------------
% Get Drag Coefficient based on Mach number
%--------------------------------------------------------------------------
function Cd_out = getDragCoefficient(Mach)
    if Mach < 0.8
        Cd_out = 0.20;
    elseif Mach < 1.2
        Cd_out = 0.30;  % Transonic bump
    elseif Mach < 5.0
        Cd_out = 0.25;  % Supersonic
    else
        Cd_out = 0.20;  % Hypersonic
    end
end

%--------------------------------------------------------------------------
% 1976 Standard Atmosphere (simplified up to ~47 km, clamp beyond)
%--------------------------------------------------------------------------
function [T, a, p, rho] = standardAtmos1976(h)
    % Returns Temperature (K), Speed of Sound (m/s), Pressure (Pa), Density (kg/m^3)
    % for geometric altitude h (m)

    g0    = 9.80665;        % [m/s^2] Sea-level gravity
    R_air = 287.0531;       % [J/(kg·K)] Specific gas constant for air
    gamma = 1.4;            % Ratio of specific heats

    % Define atmospheric layers
    hBase = [    0,  11000, 20000, 32000, 47000 ];            % [m]
    TBase = [288.15, 216.65, 216.65, 228.65, 270.65];        % [K]
    pBase = [101325, 22632.06, 5474.889, 868.0187, 110.9063];% [Pa]
    L     = [-0.0065, 0.0, 0.0010, 0.0028];                  % [K/m]

    if h < hBase(2)
        % Layer 1: 0-11 km, lapse rate = -6.5 K/km
        Tb = TBase(1);
        pb = pBase(1);
        Lb = L(1);
        hb = hBase(1);
        T  = Tb + Lb * (h - hb);
        p  = pb * (T / Tb)^(-g0 / (Lb * R_air));
    elseif h < hBase(3)
        % Layer 2: 11-20 km, isothermal
        Tb = TBase(2);
        pb = pBase(2);
        Lb = L(2);
        hb = hBase(2);
        T  = Tb;
        p  = pb * exp(-g0 * (h - hb) / (R_air * Tb));
    elseif h < hBase(4)
        % Layer 3: 20-32 km, lapse rate = +1.0 K/km
        Tb = TBase(3);
        pb = pBase(3);
        Lb = L(3);
        hb = hBase(3);
        T  = Tb + Lb * (h - hb);
        p  = pb * (T / Tb)^(-g0 / (Lb * R_air));
    elseif h < hBase(5)
        % Layer 4: 32-47 km, lapse rate = +2.8 K/km
        Tb = TBase(4);
        pb = pBase(4);
        Lb = L(4);
        hb = hBase(4);
        T  = Tb + Lb * (h - hb);
        p  = pb * (T / Tb)^(-g0 / (Lb * R_air));
    else
        % Above 47 km, clamp temperature and pressure
        Tb = TBase(5);
        pb = pBase(5);
        T  = Tb;
        p  = pb;
    end

    % Calculate density
    rho = p / (R_air * T);

    % Calculate speed of sound
    a = sqrt(gamma * R_air * T);
end
