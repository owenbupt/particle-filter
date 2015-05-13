%%
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

    sigma_r_sen0 = 500;
    sigma_theta_sen0 = 2;
    sigma_rdot_sen0 = 6;

    x_sen1 = -30000;
    y_sen1 = 0;

    sigma_r_sen1 = 500;
    sigma_theta_sen1 = 2;
    sigma_rdot_sen1 = 6;

    x_sen2 = 15000;
    y_sen2 = 26000;

    sigma_r_sen2 = 500;
    sigma_theta_sen2 = 2;
    sigma_rdot_sen2 = 6;

    x_sen3 = 15000;
    y_sen3 = -26000;

    sigma_r_sen3 = 250;
    sigma_theta_sen3 = 1;
    sigma_rdot_sen3 = 3;


    x_sen = [x_sen0, x_sen1, x_sen2, x_sen3];
    y_sen = [y_sen0, y_sen1, y_sen2, y_sen3];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sensors sensible area
    R_min = 25000;
    R_max = 100000;

    % the detection probability
    Pd = 0.85;

    % false alarm probability
    Pf = 0.47;
    % false alarm rate
    beta = 1.0 * 10 ^ -8;
    N = 2000;
    % numbers of particles


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
    x_target(1, :) = [x_target0, x_target1];
    y_target(1, :) = [y_target0, y_target1];
    vx_target(1, :) = [vx_target0, vx_target1];
    vy_target(1, :) = [vy_target0, vy_target1];
    x_target_hat = x_target;
    y_target_hat = y_target;
    vx_target_hat = vx_target;
    vy_target_hat = vy_target;
    
    x_kpi = zeros(K, N, T, 4);
    vx = [randn, randn];
    vy = [randn, randn];
    Qr = 1;
    Qrdot = 1;
    Qtheta = 1;

    % Initialization
    weight_pk = zeros(K, N);
    for p = 1:N
        for i = 1:T
            x_kpi(1, p, i, :)= [x_target_hat(1, i); vx_target_hat(1, i); y_target_hat(1, i); vy_target_hat(1, i)];
        end
        weight_pk(1, p) = 1 / N;
    end
    for k = 2:K
    % Prediction
        temp = zeros(4, T);
        for i = 1:T
            temp(:, i) = x_kpi(k - 1, p, i, :);
        end
        temp(:,1)
        for i = 1:T
            x_kpi(k, p, i, :) = [1 Ts 0 0; 0 1 0 0; 0 0 1 Ts; 0 0 0 1]*temp(:, i) + [Ts^2 / 2, 0; Ts, 0; 0, Ts^2 / 2; 0 Ts] * [vx(i); vy(i)];
        end

    % Weight update
        for p = 1:N
            for m = 1:M
                weight_pmk = weight_cal(p, m, k);
                weight_pk(k, p) = weight_pk(k, p) * weight_pmk;
            end
        end
        %%
        % normalization
        for p = 1:N
            weight_pk(k, p) = weight_pk(k, p) / sum(weight_pk(k, :));
        end

        % resampling
        temp_wei = zeros(1, N);
        temp_S = zeros(1, N);
        for i = 1:N
            temp_u = rand;
            sum_wei = 0;
            for j = 1:N
                sum_wei = sum_wei + weight_pk(k, j);
                if sum_wei > temp_u
                    temp_wei(i) = weight_pk(k, j);
                    temp_S = S(k, p, :);
                    break;
                end
            end
        end
        weight_pk(k,:) = temp_wei;
        S(k, p, :) = temp_S;
    end



