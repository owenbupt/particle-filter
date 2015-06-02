%% fgreedy: the function form of the greedy algorthim
function [greedy_cost] = fgreedy(k_start, Achose, S_kpi, weight_kp, x_target, vx_target, y_target, vy_target, x_target_hat, vx_target_hat, y_target_hat, vy_target_hat)
     %%
    % Define the related parameters
    % Information for the whole
    % sampling interval Ts
    Ts = 2;
    % Sampling time
    Times = 100;
    % Sampling numbers
    K = Times / Ts;
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





    % A = zeros(10, 4);
%     Achose = zeros(1, 4);
    A_km = zeros(K, M);
    % for i = 1:10
    %     a = dec2bin(i, 4);
    %     A(i,:) = a;
    % end
    greedy_cost = 0;
    chosn_m = zeros(1, T);
    for k = k_start:K
        for t = 1:T
            dis_min = exp(100000); 
            for m = 1:M
                dis_tm = r(x_sen(m), y_sen(m), x_target(k, t), y_target(k, t));
                if dis_tm < dis_min
                    chosn_m(t) = m;
                    dis_min = dis_tm;
                end
            end
            Achose(chosn_m(t)) = 1;
            A_km(k, :) = Achose(1, :);
%             A(k, :) = Achose(1, :);
        end
        [x_target_hat, vx_target_hat, y_target_hat, vy_target_hat, S_kpi, weight_kp] = fpf(k, Achose, S_kpi, weight_kp, x_target, vx_target, y_target, vy_target, x_target_hat, vx_target_hat, y_target_hat, vy_target_hat);
        greedy_cost = greedy_cost + one_step_cost(k, A_km, x_target_hat, y_target_hat, x_target, y_target);
    end

end