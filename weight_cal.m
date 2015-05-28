function weight_pmk = weight_cal(S_kpi, x_sen, y_sen, p, m, k, Qr, Qtheta, Qrdot, T, Hy,  false_rate_beta, Pd)
    %%
    % Problems: p is not used, so the function is the same result for
    % all P, I guess there are some problems for this condition
    % calculate weight using 8, 9, 10

    weight_pmk = 0;

    [N_hm, no_use] = size(Hy);
    dm = ones(1, T);
    em = ones(1, T);
    nf = N_hm - T;
    nd = 0;
    r_hat = zeros(1, T);
    theta_hat = zeros(1, T);
    rdot_hat = zeros(1, T);
    x_hat_kpi = zeros(1, T);
    y_hat_kpi = zeros(1, T);
    vx_hat_kpi = zeros(1, T);
    vy_hat_kpi = zeros(1, T);

    % x_sen = x_sen / 10^3;
    % y_sen = y_sen / 10^3;
    % for i = 1:T
    %     x_hat_kpi(i) = S_kpi(k, p, i, 1) / 10^3;
    %     vx_hat_kpi(i) = S_kpi(k, p, i, 2);
    %     y_hat_kpi(i) = S_kpi(k, p, i, 3) / 10^3;
    %     vy_hat_kpi(i) = S_kpi(k, p, i, 4);
    % end
    for i = 1:T
        x_hat_kpi(i) = S_kpi(k, p, i, 1);
        vx_hat_kpi(i) = S_kpi(k, p, i, 2);
        y_hat_kpi(i) = S_kpi(k, p, i, 3);
        vy_hat_kpi(i) = S_kpi(k, p, i, 4);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:T
        r_hat(i) = sqrt((x_hat_kpi(i) - x_sen(m))^2 + (y_hat_kpi(i) - y_sen(m))^2);
        theta_hat(i) = atan((y_hat_kpi(i) - y_sen(m)) / (x_hat_kpi(i) - y_sen(m)));
        % the result is radians
        rdot_hat(i) = ((x_hat_kpi(i) - x_sen(m)) * vx_hat_kpi(i) + (y_hat_kpi(i) - y_sen(m)) * vy_hat_kpi(i)) / sqrt((x_hat_kpi(i) - x_sen(m))^2 + (y_hat_kpi(i) - y_sen(m))^2);
        if ~any(Hy(i, :)) == 0
            nd = nd + 1;    
        end
    end
    % hypothese index for all kinds of conditions, here is a small problem, it is only ok when there are two targets
    t = 1;
    Hl = zeros(2, N_hm * (N_hm - 1));
    for i = 1:N_hm
        for j = 1:N_hm
            if i ~= j
                Hl(1, t) = i;
                Hl(2, t) = j;
                t = t + 1;
            end
        end
    end

    [~, N_allindex] = size(Hl);
    for l = 1:N_allindex
        weight_l = false_rate_beta^(nf) * (1 - Pd)^(T - nd) * Pd^nd; 
        for i = 1:T
            r_im_4 = (sqrt((x_sen(m) - x_hat_kpi(i))^2 + (y_sen(m) - y_hat_kpi(i))^2))^4;
            sigma_r_i = Qr * sqrt(r_im_4);
            sigma_rdot_i = Qrdot * sqrt(r_im_4);
            sigma_theta_i = Qtheta * sqrt(r_im_4);
            sigma_r_i = 500;
            sigma_rdot_i = 6;
            sigma_theta_i = 2 / 180 * pi;
            dm(i) = ((r_hat(i) - Hy(Hl(i, l), 1)) / sigma_r_i)^2 + ((theta_hat(i) - Hy(Hl(i, l), 2)) / sigma_theta_i)^2 + ((rdot_hat(i) - Hy(Hl(i, l), 3)) / sigma_rdot_i)^2;
            em(i) = (exp(-dm(i) / 2)) / (sqrt((2*pi)^3) * sigma_r_i * sigma_theta_i * sigma_rdot_i)*10^25;
            %% there is differences for different number of Hl
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~any(Hy(Hl(i, l), :)) == 0
                weight_l = weight_l * em(i);
            end
        end
        weight_pmk = weight_pmk + weight_l;
    end
end
