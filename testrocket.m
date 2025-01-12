clear
clc
close all

% function twoStageGravityTurnLEO_withStage1Removal()
%--------------------------------------------------------------------------
% Two-Stage 2D Gravity Turn around a Spherical Earth (No Rotation),
% using 1976 Standard Atmosphere, aiming for LEO.
%
% CHANGES/KEY POINTS:
%   1) The first-stage (Stage 1) has a 22-ton dry mass (22,000 kg).
%   2) We remove Stage 1 (its dry mass + leftover prop) at t=150 s.
%   3) The simulation continues until the payload is stable in LEO
%      or we hit a time limit.
%   4) We plot altitude, speed, and heat flux for the final payload.
%
% STILL SIMPLIFIED:
%   - No Earth rotation, no inclination, no real throttle schedules.
%   - Heat flux uses a Sutton-Graves type formula (often for reentry),
%     here just to show how to track heating over time in an ascent sim.
%--------------------------------------------------------------------------

    %------------------------------
    % EARTH & SIM PARAMETERS
    %------------------------------
    Re          = 6371000;            % [m] mean Earth radius
    mu          = 3.986004418e14;     % [m^3/s^2] Earth GM
    tMax        = 2000;               % [s] max simulation time
    dt          = 1.0;                % [s] time step
    pitchKick   = 10;                 % [s] time for initial pitch maneuver
    stage1Sep   = 150;                % [s] time of stage-1 separation
    targetAlt   = 200e3;              % [m] target LEO altitude
    flattenAlt  = 120e3;              % [m] altitude above which we pitch nearly horizontal

    %------------------------------
    % STAGE PARAMETERS
    %------------------------------
    % Stage 1
    thrust1        = 7.6e6;       % [N] ~7.6 MN
    mdot1          = 2200;        % [kg/s] mass flow
    dryMass1       = 22000;       % [kg] 22-ton dry mass
    propMass1      = 330e3;       % [kg] Stage1 prop
    stage1Reserve  = 20e3;        % [kg] leftover unburned prop (for booster landing, etc.)

    % Stage 2
    thrust2        = 0.9e6;       % [N] ~0.9 MN
    mdot2          = 150;         % [kg/s]
    dryMass2       = 20000;       % [kg]
    propMass2      = 100e3;       % [kg]

    % Combined initial mass (Stage1 + Stage2 + payload if any)
    m0 = (dryMass1 + propMass1) + (dryMass2 + propMass2);

    % Rocket geometry
    rocket_d  = 3.7;                         % [m] diameter
    crossArea = pi*(rocket_d/2)^2;           % [m^2] cross-sectional area

    %------------------------------
    % STATE VARIABLES
    %  (r, theta) in spherical coords
    %  (vr, vth)  radial, tangential velocity
    %------------------------------
    t   = 0;
    r   = Re;   % start at Earth surface
    th  = 0;
    vr  = 0;
    vth = 0;
    m   = m0;

    % Track each stage's prop
    stage1PropRem = propMass1;
    stage2PropRem = propMass2;

    % For an orbit insertion check
    orbitAchieved = false;

    % Data storage: [time, r, alt, vr, vth, Mach, heatFlux, stageID]
    simData = [];

    %--------------------------------------------------------------------------
    % MAIN LOOP
    %--------------------------------------------------------------------------
    while (t < tMax) && (~orbitAchieved)
        % altitude
        alt = r - Re;
        if alt < 0
            fprintf('Rocket below ground!\n');
            break;
        end

        % 1) Check if we achieved orbit
        [inOrbit, ~] = checkOrbitInsertion(r, vr, vth, mu, targetAlt);
        if inOrbit
            orbitAchieved = true;
            fprintf('*** Orbit Achieved at t=%.1f s ***\n', t);
            break;
        end

        % 2) Determine active stage
        if t < stage1Sep
            stageID       = 1;
            thrustCurrent = thrust1;
            mdotCurrent   = mdot1;
        else
            stageID = 2;
            if stage2PropRem > 0
                thrustCurrent = thrust2;
                mdotCurrent   = mdot2;
            else
                thrustCurrent = 0;  % coasting
                mdotCurrent   = 0;
            end
        end

        % 3) local gravity
        gLocal = mu / (r^2);

        % 4) 1976 Standard Atmosphere
        [T, aLoc, ~, rho] = standardAtmos1976(alt);

        % 5) speed & Mach
        speed = sqrt(vr^2 + vth^2);
        Mach  = (aLoc>1e-6)*(speed / max(aLoc,1e-6));

        % 6) Mach-based drag coefficient
        Cd = getDragCoefficient(Mach);

        % 7) Drag force
        F_drag = 0.5 * Cd * crossArea * rho * speed^2;

        % 8) Thrust direction logic (gravity turn)
        if t < pitchKick
            % purely radial at liftoff
            thrustRad = thrustCurrent;
            thrustTan = 0;
        else
            % zero-lift turn => align thrust w/ velocity
            if speed > 1e-3
                velAngle  = atan2(vth, vr);
            else
                velAngle  = 0;
            end
            thrustRad = thrustCurrent * cos(velAngle);
            thrustTan = thrustCurrent * sin(velAngle);

            % If stage=2 & alt>flattenAlt => push mostly horizontal
            if (stageID == 2) && (alt > flattenAlt)
                thrustRad = 0.1 * thrustCurrent * cos(velAngle);
                thrustTan = 1.0 * thrustCurrent * sin(velAngle);
            end
        end

        % 9) Decompose drag in polar coords
        if speed < 1e-3
            dragRad = 0;
            dragTan = 0;
        else
            dragRad = -F_drag * (vr/speed);
            dragTan = -F_drag * (vth/speed);
        end

        % 10) Net radial/tangential force
        FR = thrustRad + dragRad - m*gLocal;
        FT = thrustTan + dragTan;

        % 11) Accelerations in polar coords
        ar = (FR/m) - (vth^2 / r);
        at = (FT/m) + (vr*vth / r);

        % 12) Euler integration
        vr_new  = vr  + ar*dt;
        vth_new = vth + at*dt;
        r_new   = r   + vr_new*dt;
        th_new  = th  + (vth_new / r_new)*dt;

        % 13) Propellant usage
        dm = mdotCurrent * dt;
        if stageID == 1
            % Do not burn below stage1Reserve
            usableProp = max(0, stage1PropRem - stage1Reserve);
            if usableProp > dm
                stage1PropRem = stage1PropRem - dm;
                m_new = m - dm;
            else
                dm_used = usableProp;  % might be 0 if at reserve
                stage1PropRem = stage1PropRem - dm_used;
                m_new = m - dm_used;
            end
        else
            % Stage 2
            if stage2PropRem > dm
                stage2PropRem = stage2PropRem - dm;
                m_new = m - dm;
            else
                dm_used = stage2PropRem;
                stage2PropRem = 0;
                m_new = m - dm_used;
            end
        end

        % If crossing stage1Sep, remove stage1 dry mass + leftover prop
        if (t < stage1Sep) && (t + dt >= stage1Sep)
            leftoverStage1 = dryMass1 + stage1PropRem;  % all stage1 mass left
            m_new = m_new - leftoverStage1;
        end

        % 14) Heat flux (Sutton-Graves approximation)
        heatFlux = 0;
        if speed > 0 && rho > 1e-10
            kSG   = 1.83;    % empirical factor
            R_n   = 1.0;     % [m] nose radius
            heatFlux = kSG * sqrt(rho / R_n) * (speed^3.05);
        end

        % 15) Advance time
        t_next = t + dt;

        % 16) Store data
        simData(end+1,:) = [ t, r, alt, vr, vth, Mach, heatFlux, stageID ];

        % 17) Update states
        vr  = vr_new;
        vth = vth_new;
        r   = r_new;
        th  = th_new;
        m   = max(m_new, 0);
        t   = t_next;

        if r < Re
            fprintf('Below Earth surface, break.\n');
            break;
        end
    end

    % Final data point
    simData(end+1,:) = [ t, r, (r - Re), vr, vth, Mach, heatFlux, stageID ];

    % Check orbit if not flagged
    if ~orbitAchieved
        [inOrbit, ~] = checkOrbitInsertion(r, vr, vth, mu, targetAlt);
        if inOrbit
            fprintf('*** Orbit Achieved at final step (t=%.1f s) ***\n', t);
        else
            fprintf('End of sim, orbit not achieved.\n');
        end
    end

    %--------------------------------------------------------------------------
    % PLOTTING: Altitude, Speed, Heat Flux vs Time
    %--------------------------------------------------------------------------
    timeVec   = simData(:,1);
    altVec    = simData(:,3);
    vrVec     = simData(:,4);
    vthVec    = simData(:,5);
    heatFlux  = simData(:,7);

    speedVec  = sqrt(vrVec.^2 + vthVec.^2);

    figure('Name','Payload Ascent to LEO','Color','w','Position',[100 100 900 600]);

    % Altitude vs Time
    subplot(3,1,1);
    plot(timeVec, altVec*1e-3,'b-','LineWidth',2);
    xlabel('Time (s)'); ylabel('Altitude (km)');
    title('Altitude vs. Time'); grid on;

    % Velocity vs Time
    subplot(3,1,2);
    plot(timeVec, speedVec,'r-','LineWidth',2);
    xlabel('Time (s)'); ylabel('Speed (m/s)');
    title('Speed vs. Time'); grid on;

    % Heat Flux vs Time
    subplot(3,1,3);
    plot(timeVec, heatFlux,'m-','LineWidth',2);
    xlabel('Time (s)'); ylabel('Heat Flux (W/m^2)');
    title('Heat Flux vs. Time (Sutton-Graves Approx)'); grid on;
% end


%% HELPER FUNCTIONS

%--------------------------------------------------------------------------
% Check if orbit insertion is achieved near target altitude
%--------------------------------------------------------------------------
function [inOrbit, info] = checkOrbitInsertion(r, vr, vth, mu, targetAlt)
    % We define "orbit achieved" if:
    %   1) altitude ~ targetAlt ± 10 km
    %   2) radial velocity ~ 0 ± 50 m/s
    %   3) tangential speed ~ sqrt(mu/r) ± 100 m/s
    Re     = 6371000;
    alt    = r - Re;
    vCirc  = sqrt(mu/r);
    tolAlt = 10e3;
    tolVr  = 50;
    tolVth = 100;

    if abs(alt - targetAlt) < tolAlt && abs(vr) < tolVr && abs(vth - vCirc) < tolVth
        inOrbit = true;
        info    = "Orbit Achieved";
    else
        inOrbit = false;
        info    = "Not in orbit";
    end
end

%--------------------------------------------------------------------------
% Mach-based drag coefficient
%--------------------------------------------------------------------------
function Cd_out = getDragCoefficient(Mach)
    if Mach < 0.8
        Cd_out = 0.20;
    elseif Mach < 1.2
        Cd_out = 0.30;  % transonic bump
    elseif Mach < 5.0
        Cd_out = 0.25;  % supersonic
    else
        Cd_out = 0.20;  % hypersonic
    end
end

%--------------------------------------------------------------------------
% 1976 Standard Atmosphere (simplified up to ~47 km, clamp above)
%--------------------------------------------------------------------------
function [T, a, p, rho] = standardAtmos1976(h)
    g0    = 9.80665;    
    R_air = 287.0531;  
    gamma = 1.4;

    hBase = [    0,  11000, 20000, 32000, 47000 ];
    TBase = [288.15, 216.65, 216.65, 228.65, 270.65];
    pBase = [101325, 22632.06, 5474.889, 868.0187, 110.9063];
    L     = [-0.0065, 0.0, 0.0010, 0.0028];

    if h < hBase(2)
        Tb = TBase(1); pb = pBase(1); Lb = L(1); hb = hBase(1);
        T  = Tb + Lb*(h - hb);
        p  = pb*(T/Tb)^(-g0/(Lb*R_air));
    elseif h < hBase(3)
        Tb = TBase(2); pb = pBase(2); Lb = L(2); hb = hBase(2);
        T  = Tb;
        p  = pb*exp(-g0*(h-hb)/(R_air*Tb));
    elseif h < hBase(4)
        Tb = TBase(3); pb = pBase(3); Lb = L(3); hb = hBase(3);
        T  = Tb + Lb*(h - hb);
        p  = pb*(T/Tb)^(-g0/(Lb*R_air));
    elseif h < hBase(5)
        Tb = TBase(4); pb = pBase(4); Lb = L(4); hb = hBase(4);
        T  = Tb + Lb*(h - hb);
        p  = pb*(T/Tb)^(-g0/(Lb*R_air));
    else
        % clamp above 47 km
        Tb = TBase(5);
        pb = pBase(5);
        T  = Tb;
        p  = pb;
    end

    rho = p/(R_air*T);
    a   = sqrt(gamma*R_air*T);
end
