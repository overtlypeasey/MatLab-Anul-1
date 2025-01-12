function twoStageGravityTurnLEO_corrected()
%--------------------------------------------------------------------------
% Two-Stage 2D Gravity Turn Simulation to LEO with Realistic Thrust Curves
%
% - Implements a two-stage rocket ascent with gravity turn
% - Uses 1976 Standard Atmosphere for atmospheric properties
% - Includes realistic thrust curves for each stage
% - Retains reserve fuel in Stage 1 for landing maneuvers
% - Removes Stage 1 mass correctly at separation
% - Prevents the rocket from falling back to the ground by simulating a solid ground
% - Ends simulation when payload is in stable LEO
% - Plots Altitude, Velocity, and Heat Flux vs Time
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
    dt          = 1.0;                % [s] Time step
    pitchKick   = 10;                 % [s] Time to begin gravity turn
    stage1Sep   = 150;                % [s] Time of Stage 1 separation
    targetAlt   = 200e3;              % [m] Target LEO altitude
    flattenAlt  = 120e3;              % [m] Altitude to begin thrust flattening

    %------------------------------
    % STAGE PARAMETERS
    %------------------------------
    % Stage 1
    thrust1_func   = @thrustCurveStage1;  % Function handle for Stage 1 thrust
    burnTime1      = 150;                 % [s] Total burn time for Stage 1
    mdot1          = 2200;                % [kg/s] Mass flow rate Stage 1
    dryMass1       = 22000;               % [kg] 22-ton dry mass
    propMass1      = 330000;              % [kg] Stage 1 propellant
    stage1Reserve  = 20000;               % [kg] Reserve propellant in Stage 1

    % Stage 2
    thrust2_func   = @thrustCurveStage2;  % Function handle for Stage 2 thrust
    burnTime2      = 900;                 % [s] Total burn time for Stage 2
    mdot2          = 150;                 % [kg/s] Mass flow rate Stage 2 (Adjusted below)
    dryMass2       = 5000;                % [kg] 5-ton dry mass (realistic value)
    propMass2      = 90700;               % [kg] Stage 2 propellant (~90.7 t)

    % Payload
    payloadMass    = 10000;               % [kg] Payload mass

    % Combined initial mass (Stage1 + Stage2 + Payload)
    m0 = (dryMass1 + propMass1) + (dryMass2 + propMass2) + payloadMass;

    % Rocket geometry
    rocket_d  = 3.7;                         % [m] Diameter
    crossArea = pi*(rocket_d/2)^2;           % [m^2] Cross-sectional area

    %------------------------------
    % STATE VARIABLES
    %  (r, theta) in spherical coordinates
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

    % Data storage: [time, alt, vr, vth, speed, Mach, heatFlux, stageID]
    simData = [];

    %--------------------------------------------------------------------------
    % MAIN SIMULATION LOOP
    %--------------------------------------------------------------------------
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
            if t <= burnTime1
                thrustCurrent = thrust1_func(t);
                dmCurrent = mdot1;
            else
                thrustCurrent = 0;
                dmCurrent = 0;
            end
        else
            stageID = 2;
            timeSinceStage2Start = t - stage1Sep;
            if (timeSinceStage2Start <= burnTime2) && (stage2PropRem > 0)
                thrustCurrent = thrust2_func(timeSinceStage2Start);
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

        % 8) Determine thrust direction (gravity turn logic)
        if onGround
            % While on ground, check if thrust exceeds weight
            if thrustCurrent > m * gLocal
                onGround = false;
                % Calculate initial acceleration
                netThrust = thrustCurrent - m * gLocal;
                ar = netThrust / m;
                at = 0;  % Initial tangential acceleration
            else
                % Thrust insufficient to lift off
                ar = 0;
                at = 0;
                vr = 0;
                vth = 0;
            end
        else
            % Gravity turn logic
            if t < pitchKick
                % Purely radial thrust (vertical)
                thrustRad = thrustCurrent;
                thrustTan = 0;
            else
                if speed > 1e-3
                    velAngle = atan2(vth, vr);
                else
                    velAngle = 0;
                end
                % Zero-lift turn: align thrust with velocity vector
                thrustRad = thrustCurrent * cos(velAngle);
                thrustTan = thrustCurrent * sin(velAngle);

                % After flattenAlt altitude and if Stage 2 is active, make thrust more tangential
                if (stageID == 2) && (alt > flattenAlt)
                    thrustRad = 0.1 * thrustCurrent * cos(velAngle);
                    thrustTan = 1.0 * thrustCurrent * sin(velAngle);
                end
            end

            % Decompose drag into radial and tangential components
            if speed < 1e-3
                dragRad = 0;
                dragTan = 0;
            else
                dragRad = -F_drag * (vr / speed);
                dragTan = -F_drag * (vth / speed);
            end

            % Net radial and tangential forces
            FR = thrustRad + dragRad - m * gLocal;    % Radial force
            FT = thrustTan + dragTan;                % Tangential force

            % Accelerations in polar coordinates
            ar = (FR / m) - (vth^2 / r);
            at = (FT / m) + ((vr * vth) / r);
        end

        % 9) Calculate heat flux (Sutton-Graves approximation)
        if (speed > 0) && (rho > 1e-10) && (~onGround)
            kSG    = 1.83;        % [W/(m^2*(kg/m^3)^(1/2)*(m/s)^3.05)]
            R_nose = 1.0;         % [m] Nose radius
            heatFlux = kSG * sqrt(rho / R_nose) * (speed^3.05);
        else
            heatFlux = 0;
        end

        % 10) Calculate net forces only if not on ground
        if ~onGround
            % Update velocities and positions
            vr_new  = vr  + ar * dt;
            vth_new = vth + at * dt;
            r_new   = r   + vr_new * dt;
            th_new  = th  + (vth_new / r_new) * dt;
        else
            % On ground, no movement
            vr_new  = 0;
            vth_new = 0;
            r_new   = r;
            th_new  = th;
            heatFlux = 0;  % No heating on ground
        end

        % 11) Update mass based on propellant consumption
        if ~onGround
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
        else
            % On ground, do not burn propellant
            m_new = m;
        end

        % 12) Stage separation at stage1Sep (t = 150 s)
        if (t < stage1Sep) && (t + dt >= stage1Sep)
            % Remove Stage 1 dry mass and reserve propellant
            leftoverStage1 = dryMass1 + stage1Reserve;  % all Stage1 mass left
            m_new = m_new - leftoverStage1;
            fprintf('Stage 1 separated at t=%.1f s, Altitude=%.1f km\n', t, alt/1e3);
        end

        % 13) Update state variables for next iteration
        if ~onGround
            vr  = vr_new;
            vth = vth_new;
            r   = r_new;
            th  = th_new;
        end
        m   = max(m_new, 0);  % Prevent negative mass
        t   = t + dt;

        % 14) Store simulation data
        simData(end+1,:) = [ t, alt, vr, vth, speed, Mach, heatFlux, stageID ];
    end

    % Final data point
    simData(end+1,:) = [ t, (r_new - Re), vr_new, vth_new, sqrt(vr_new^2 + vth_new^2), ...
                         Mach, heatFlux, stageID ];

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
    timeVec   = simData(:,1);
    altVec    = simData(:,2);
    vrVec     = simData(:,3);
    vthVec    = simData(:,4);
    speedVec  = simData(:,5);
    MachVec   = simData(:,6);
    heatFlux  = simData(:,7);
    stageID   = simData(:,8);

    % Display a summary table
    fprintf('\n   t(s)  Alt(km)   vr(m/s)  vth(m/s)  Speed(m/s)  Mach  HeatFlux(W/m^2)  Stage\n');
    fprintf('-------------------------------------------------------------------------------\n');
    for i = 1 : floor(length(timeVec)/10) : length(timeVec)
        fprintf('%6.1f  %7.2f   %8.1f  %8.1f  %10.1f   %4.2f      %8.1f      %d\n', ...
                timeVec(i), altVec(i)/1e3, vrVec(i), vthVec(i), speedVec(i), MachVec(i), heatFlux(i), stageID(i));
    end

    % Plotting
    figure('Name','Payload Ascent to LEO','Color','w','Position',[100 100 900 800]);

    % Altitude vs Time
    subplot(3,1,1);
    plot(timeVec, altVec/1e3, 'b-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Altitude (km)');
    title('Altitude vs. Time');
    grid on;

    % Velocity vs Time
    subplot(3,1,2);
    plot(timeVec, speedVec, 'r-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Speed (m/s)');
    title('Speed vs. Time');
    grid on;

    % Heat Flux vs Time
    subplot(3,1,3);
    plot(timeVec, heatFlux, 'm-', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Heat Flux (W/m^2)');
    title('Heat Flux vs. Time (Sutton-Graves Approx)');
    grid on;

end

%% HELPER FUNCTIONS

%--------------------------------------------------------------------------
% Thrust Curve for Stage 1
% Implements a realistic thrust curve with ramp-up and ramp-down
%--------------------------------------------------------------------------
function thrust = thrustCurveStage1(t)
    % Define time segments for Stage 1 thrust
    % Ramp-up: 0-10 s
    % Constant thrust: 10-140 s
    % Ramp-down: 140-150 s
    if t < 10
        thrust = 7.6e6 * (t / 10);  % Linear ramp-up
    elseif t < 140
        thrust = 7.6e6;              % Constant thrust
    elseif t <= 150
        thrust = 7.6e6 * (1 - (t - 140) / 10);  % Linear ramp-down
    else
        thrust = 0;                  % No thrust after burnout
    end
end

%--------------------------------------------------------------------------
% Thrust Curve for Stage 2
% Implements a realistic thrust curve with ramp-up and ramp-down
% Adjusted to reflect realistic Isp
%--------------------------------------------------------------------------
function thrust = thrustCurveStage2(t)
    % Adjusted thrust curve for Stage 2 to match realistic Isp of 300 s
    % Thrust = 441300 N (calculated based on mdot2=150 kg/s and Isp=300 s)
    
    % Define time segments for Stage 2 thrust
    % Ramp-up: 0-20 s
    % Constant thrust: 20-222 s (Adjusted burn time based on propMass2 and mdot2)
    % Ramp-down: 222-242 s
    if t < 20
        thrust = 441300 * (t / 20);   % Linear ramp-up
    elseif t < 222
        thrust = 441300;              % Constant thrust
    elseif t <= 242
        thrust = 441300 * (1 - (t - 222) / 20);  % Linear ramp-down
    else
        thrust = 0;                  % No thrust after burnout
    end
end

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
