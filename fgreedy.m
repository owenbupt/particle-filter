%% fgreedy: the function form of the greedy algorthim
function [greedy_cost] = fgreedy(k_start, Achose, S_kpi, weight_kp, x_target, vx_target, y_target, vy_target, x_target_hat, vx_target_hat, y_target_hat, vy_target_hat)
     %%
    % Define the related parameters
    % Information for the whole
    % sampling interval Ts
    Ts = 2;
    % Sampling time
    Times = 30;
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
    Pd = 0.8;
    %%%%%%%%%%%%%change it to 1 for debug

    % false alarm probability
    Pf = 0.47;
    % false alarm rate
    false_rate_beta = 1.0 * 10 ^ -8;
    N = 200;
    % numbers of particles

    % vx = [randn, randn];
    % vy = [randn, randn];
    Qr = 1;
    Qrdot = 1;
    Qtheta = 1;
    false_mean = 1;



    % A = zeros(10, 4);
    Achose = zeros(1, 4);
    A_km = zeros(K, M);
    % for i = 1:10
    %     a = dec2bin(i, 4);
    %     A(i,:) = a;
    % end
    greedy_cost = 0;
    for k = k_start:K
        for t = 1:T
            dis_min = exp(100000); 
            for m = 1:M
                dis_tm = r(x_sen(m), y_sen(m), x_target(t), y_target(t));
                if dis_tm < dis_min
                    chosn_m(t) = m;
                    dis_min = dis_tm;
                end
            end
            Achose(chosn_m(t)) = 1;
            A(k, :) = Achose(1, :);
        end
        [x_target_hat, vx_target_hat, y_target_hat, vy_target_hat S_kpi weight_kp] = fpf(k, Achose, S_kpi, weight_kp, x_target, vx_target, y_target, vy_target, x_target_hat, vx_target_hat, y_target_hat, vy_target_hat);
        greedy_cost = greedy_cost + one_step_cost(k, A_km, x_target_hat, y_target_hat, x_target, y_target);
    end

end