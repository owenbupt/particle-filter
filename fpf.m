function [x_target_hat, vx_target_hat, y_target_hat, vy_target_hat, S_kpi, weight_kp] = fpf(k, Achose, S_kpi, weight_kp, x_target, vx_target, y_target, vy_target, x_target_hat, vx_target_hat, y_target_hat, vy_target_hat)
    % if nargin == 0
        %%
        % Define the related parameters
        % Information for the whole
        % sampling interval Ts 
        Ts = 2;
        % Sampling time
        % Sampling numbers
        T = 2; % targets numbers
        M = 4; % sensors numbers

        % Information for sensors, in this example there are 4 sensors.
        x_sen0 = 0;
        y_sen0 = 0;

        x_sen1 = -30000;
        y_sen1 = 0;

        x_sen2 = 15000;
        y_sen2 = 26000;


        x_sen3 = 15000;
        y_sen3 = -26000;



        x_sen = [x_sen0, x_sen1, x_sen2, x_sen3];
        y_sen = [y_sen0, y_sen1, y_sen2, y_sen3];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sensors sensible area
        R_min = 0;
        R_max = 100000000;

        % the detection probability
        Pd = 0.8;
        %%%%%%%%%%%%%change it to 1 for debug

        % false alarm probability
    
        % false alarm rate
        false_rate_beta = 1.0 * 10 ^ -8;
        N = 200;
        % numbers of particles


        % % targets Information
        % x_target0 = 2500;
        % y_target0 = 25000;
        % vx_target0 = 0;
        % vy_target0 = -222;
        % x_target1 = -5000;
        % y_target1 = -20000;
        % vx_target1 = 120;
        % vy_target1 = 0;
        % x_target = zeros(K, T);
        % vx_target = zeros(K, T);
        % y_target = zeros(K, T);
        % vy_target = zeros(K, T);
        % x_target(1, :) = [x_target0, x_target1];
        % y_target(1, :) = [y_target0, y_target1];
        % vx_target(1, :) = [vx_target0, vx_target1];
        % vy_target(1, :) = [vy_target0, vy_target1];
        % x_target_hat = x_target;
        % y_target_hat = y_target;
        % vx_target_hat = vx_target;
        % vy_target_hat = vy_target;

        % S_kpi = zeros(K, N, T, 4);
        % vx = [randn, randn];
        % vy = [randn, randn];
        Qr = 1;
        Qrdot = 1;
        Qtheta = 1;
        false_mean = 1;
        % weight_kp = ones(K, N);
    % end

    % plot(1:K, vx_target(:,1));
    % plot(x_target(:, 1), y_target(:, 1));

    %%
    % Initialization
    % there is some problems in this part, the particles are not genearated by
    % random,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ao jiao de fen ge xian%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Prediction
    for p = 1:N
        temp = zeros(4, T);
        for i = 1:T
            temp(:, i) = S_kpi(k - 1, p, i, :);
        end
        for i = 1:T
            S_kpi(k, p, i, :) = [1 Ts 0 0; 0 1 0 0; 0 0 1 Ts; 0 0 0 1]*temp(:, i) + [Ts^2 / 2, 0; Ts, 0; 0, Ts^2 / 2; 0 Ts] * [randn; randn];
        end
    end
    % Weight update
    % Problems: Only for the condition there are 4 sensors, try to find the matlab function to use a variables with a number
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% add the action status for the m chosen%%%%
    for m = 1:M
        if Achose(m) == 1
            FalseArr = false_alarm(R_min, R_max, false_mean);

            Hy = detected(x_target, y_target, vx_target, vy_target, x_sen, y_sen, k, m, T, Pd, R_min, R_max, FalseArr);
            if m == 1
                Hy1 = Hy;
            elseif m == 2
                Hy2 = Hy;
            elseif m == 3
                Hy3 = Hy;
            else
                Hy4 = Hy;
            end
        end
    end

    for p = 1:N
        for m = 1:M
            if Achose(m) == 1
                if m == 1
                    Hy = Hy1;
                    weight_pmk = weight_cal(S_kpi, x_sen, y_sen, p, m, k, Qr, Qtheta, Qrdot, T, Hy, false_rate_beta, Pd);
                    weight_kp(k, p) = weight_kp(k, p) * weight_pmk;
                elseif m == 2
                    Hy = Hy2;
                    weight_pmk = weight_cal(S_kpi, x_sen, y_sen, p, m, k, Qr, Qtheta, Qrdot, T, Hy, false_rate_beta, Pd);
                    weight_kp(k, p) = weight_kp(k, p) * weight_pmk;
                elseif m == 3
                    Hy = Hy3;
                    weight_pmk = weight_cal(S_kpi, x_sen, y_sen, p, m, k, Qr, Qtheta, Qrdot, T, Hy, false_rate_beta, Pd);
                    weight_kp(k, p) = weight_kp(k, p) * weight_pmk;
                elseif m == 4
                    Hy = Hy4;
                    weight_pmk = weight_cal(S_kpi, x_sen, y_sen, p, m, k, Qr, Qtheta, Qrdot, T, Hy, false_rate_beta, Pd);
                    weight_kp(k, p) = weight_kp(k, p) * weight_pmk;
                end
            end
        end
    end
    %%
    % normalization
    weight_kp(k, :) = weight_kp(k, :) ./ sum(weight_kp(k, :));
    % estimalation
    for i = 1:T
        x_target_hat(k, i) = sum(S_kpi(k, :, i, 1) .* weight_kp(k, :));
        vx_target_hat(k, i) = sum(S_kpi(k, :, i, 2) .* weight_kp(k, :));
        y_target_hat(k, i) = sum(S_kpi(k, :, i, 3) .* weight_kp(k, :));
        vy_target_hat(k, i) = sum(S_kpi(k, :, i, 4) .* weight_kp(k, :));
    end
    % resampling
    temp_wei = zeros(1, N);
    temp_S = zeros(1, N, T, 4);
    for i = 1:N
        temp_u = rand;
        for j = 1:N
            sum_wei = sum(weight_kp(k, 1:j));
            if sum_wei > temp_u
                temp_wei(1, i) = weight_kp(k, j);
                temp_S(1, i, :, :) = S_kpi(k, j, :, :);
                break;
            end
        end
    end
    weight_kp(k,:) = temp_wei(1, :) ./ sum(temp_wei(1,:));
    S_kpi(k, :, :, :) = temp_S(1, :, :, :);
end
