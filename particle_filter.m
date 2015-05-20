function [] = particle_filter()
    % Define the related parameters
    % Information for the whole 
    % sampling interval Ts
    Ts = 2;
    % Sampling time
    Times = 300;
    % Sampling numbers
    K = Times / Ts;
    T = 2; % targets numbers
    M = 4; % sensors numbers

    % Information for sensors, in this example there are 4 sensors.
    x_sen0 = 0;
    y_sen0 = 0;
    cost_sen0 = 0.5;
    sigma_r_sen0 = 500;
    sigma_theta_sen0 = 2;
    sigma_rdot_sen0 = 6;

    x_sen1 = -30000;
    y_sen1 = 0;
    cost_sen1 = 0.1;
    sigma_r_sen1 = 500;
    sigma_theta_sen1 = 2;
    sigma_rdot_sen1 = 6;

    x_sen2 = 15000;
    y_sen2 = 26000;
    cost_sen2 = 0.1;
    sigma_r_sen2 = 500;
    sigma_theta_sen2 = 2;
    sigma_rdot_sen2 = 6;

    x_sen3 = 15000;
    y_sen3 = -26000;
    cost_sen3 = 0.1;
    sigma_r_sen3 = 250;
    sigma_theta_sen3 = 1;
    sigma_rdot_sen3 = 3;


    x_sen = [x_sen0, x_sen1, x_sen2, x_sen3];
    y_sen = [y_sen0, y_sen1, y_sen2, y_sen3];
    c_usage = [cost_sen0, cost_sen1, cost_sen2, cost_sen3];
    a = zeros(K, M);
    u = zeros(K, M);
    % sensors sensible area
    R_min = 25000;
    R_max = 100000;
    % sensors start/stop cost
    cost_start = 0.1;
    % the detection probability
    Pd = 0.85;

    % false alarm probability
    Pf = 0.47;
    % false alarm rate
    beta = 1.0 * 10 ^ -8;
    alpha = 0.5;
    N = 2000;
    % numbers of particles
    
    H = zeros(K, M, T, 3);


    % targets Information
    x_target0 = 2500;
    y_target0 = 25000;
    vx_target0 = 0;
    vy_target0 = -222;
    x_target1 = -5000;
    y_target1 = -20000;
    vx_target1 = 120;
    vy_target1 = 0;
    x_target = zeros(Times / Ts, 2);
    vx_target = zeros(Times / Ts, 2);
    y_target = zeros(Times / Ts, 2);
    vy_target = zeros(Times / Ts, 2);

    x_target_hat = zeros(Times / Ts, 2);
    vx_target_hat = zeros(Times / Ts, 2);
    y_target_hat = zeros(Times / Ts, 2);
    vy_target_hat = zeros(Times / Ts, 2);
    x_target(1, :) = [x_target0, x_target1];
    y_target(1, :) = [y_target0, y_target1];
    vx_target(1, :) = [vx_target0, vx_target1];
    vy_target(1, :) = [vy_target0, vy_target1];


    function [r, theta, rdot] = observe(i, m, k)
        % observations from target i and sensor m at time k.
        if sqrt((x_target(k,i) - x_sen(m))^2 + (y_target(k,i) - y_sen(m))^2) > R_min  && sqrt((x_target(k,i) - x_sen(m))^2 + (y_target(k,i) - y_sen(m))^2) < R_max;
            r = sqrt((x_target(k,i) - x_sen(m))^2 + (y_target(k,i) - y_sen(m))^2) + randn;
            theta = atan((y_target(k,i) - y_sen(m)) / (x_target(k,i) - y_sen(m))) + randn;
            % the result is radians
            r_dot = ((x_target(k,i) - x_sen(m)) * vx_target(k, i) + (y_target(k, i) - y_sen(m)) * vy_target(k, i)) / sqrt((x_target(k, i) - x_sen(m))^2 + (y_target(k, i) - y_sen(m))^2) + randn;
        end
    end

    function FalseAlarm = false_alarm()
        % notice the condition that false_num = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate the false alarm
        false_num = poissrnd(1); %poisson distrubition with 1 mean
        FalseAlarm = zeros(false_num, 3);
        for i = 1:false_num
            FalseAlarm(i, 1) = R_min + (R_max - R_min) * rand();
            FalseAlarm(i, 2) = -pi + 2 * pi * rand;
            FalseAlarm(i, 3) = 120 + (222 - 120) * rand;
            % the r_dot false is between target 0 and 1;
            % there may be some problems for the choices.
        end
    end

    function OneStepCost = one_step_cost(i, m, k)
        error_cost = 0;
        sensor_cost = 0;
        for i = 1:T
            error_cost = error_cost + (x_target_hat(k,i) - x_target(k, i))^2 + (y_target_hat(k, i) - y_target(k, i))^2;
        end
        for m = 1:M
            sensor_cost = sensor_cost + c_usage(m) * a(k, m) + cost_start * (a(k, m) - u(k, m));
        end
        OneStepCost = alpha * error_cost + sensor_cost;
    end

    function detection = det_num()
        % Calculate the number of detection in the detected probability
        u = rand;
        if rand > Pd
            detection = 0;
        else
            detection = 1;
        end
    end

    function Hy = detected(m, k)
        FalseAlarm = false_alarm();
        [false_num, no_use] = size(FalseAlarm);
        Hy = zeros(T + false_num, 3);
        for i = 1:T
            detect_result(i) = det_num()
            if detect_result(i) = 1
                Hy(i, :) = observe(i, m, k);
            else
                Hy(i, :) = zeros(1, 3);
            end
        end
        for i = 1:false_num
            Hy(T + i, :) = FalseAlarm(i);
        end
    end
    


    function weight_pmk = weight_cal(p, m, k)
        %
        % calculate weight using 8, 9, 10
        weight_pmk = 0;
        Hy = detected();
        [N_hm, no_use] = size(Hy)
        dm = ones(1, T);
        em = ones(1, T);
        nf = N_hm - T;
        nd = 0;
        for i = 1:T
            if Hy(i, :) == [0, 0, 0]
                ;
            else
                nd = nd + 1;
            end
        end
        % hypothese index for all kinds of conditions, here is a small problem, it is only ok when there are two targets
        t = 0;
        for i = 1:N_hm
            for j = 1:N_hm
                if i ~= j
                    Hl(1, t) = i;
                    Hl(2, t) = j;
                    t = t + 1;
                end
            end
        end

        [no_use, N_allindex] = size(Hl)
        for l = 1:N_allindex
            weight_l = beta^(nf) * (1 - Pd)^(T - nd) * Pd^nd; 
            for i = 1:T
                r_im_4 = ((x_sen(m) - x_target_hat(k))^2 + (y_sen(m) - y_target_hat(k))^2)^4
                sigma_r_i = Qr * r_im_4;
                sigma_rdot_i = Qrdot * r_im_4;
                sigma_theta_i = Qtheta * r_im_4;
                dm(i) = ((r_hat(m, i) - Hy(Hl(i, l), 1)) / sigma_r_i)^2 + ((theta_hat(m, i) - Hy(Hl(i, l), 2)) / sigma_theta_i)^2 + ((rdot_hat(m, i) - Hy(Hl(i, l), 3)) / sigma_rdot_i)^2;
                em(i) = (exp(-dm(i) / 2)) / (sqrt((2*pi)^3) * sigma_r_i * sigma_theta_i * sigma_rdot_i);
                if Hy(Hl(i, l), :) ~= [0 0 0]
                    weight_l = weight_l * em(i);
                end
            end
            weight_pmk = weight_pmk + weight_l;
        end
       % in this part the hypo is now consedered wrongly, it needs to change again. 





    function pf()
        % Initialization
        for i = 1:N
            S(0, i, :) = x_target_hat(k, :);
            w(0, i, :) = 1 / N * ones(1, M);
        end
        k = 1;
        for m = 1:M
            detected(m, k);
        end


        % Prediction
    end




end



