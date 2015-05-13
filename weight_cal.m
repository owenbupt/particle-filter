    function weight_pmk = weight_cal(p, m, k)
        %
        % calculate weight using 8, 9, 10
       
        weight_pmk = 0;
        Hy = detected(m, k);
        [N_hm, no_use] = size(Hy);
        dm = ones(1, T);
        em = ones(1, T);
        nf = N_hm - T;
        nd = 0;
        for i = 1:T
            for m = 1:M
                r_hat(k, m, i) = sqrt((x_target_hat(k,i) - x_sen(m))^2 + (y_target_hat(k,i) - y_sen(m))^2) + randn;
                theta_hat(k, m, i) = atan((y_target_hat(k,i) - y_sen(m)) / (x_target_hat(k,i) - y_sen(m))) + randn;
                % the result is radians
                rdot_hat(k, m, i) = ((x_target_hat(k,i) - x_sen(m)) * vx_target_hat(k, i) + (y_target_hat(k, i) - y_sen(m)) * vy_target_hat(k, i)) / sqrt((x_target_hat(k, i) - x_sen(m))^2 + (y_target_hat(k, i) - y_sen(m))^2) + randn;
            end
            if Hy(i, :) ~= [0, 0, 0]
                nd = nd + 1;
            end
        end
        % hypothese index for all kinds of conditions, here is a small problem, it is only ok when there are two targets
        t = 1;
        for i = 1:N_hm
            for j = 1:N_hm
                if i ~= j
                    Hl(1, t) = i;
                    Hl(2, t) = j;
                    t = t + 1;
                end
            end
        end

        [no_use, N_allindex] = size(Hl);
        for l = 1:N_allindex
            weight_l = beta^(nf) * (1 - Pd)^(T - nd) * Pd^nd; 
            for i = 1:T
                r_im_4 = ((x_sen(m) - x_target_hat(k))^2 + (y_sen(m) - y_target_hat(k))^2)^4;
                sigma_r_i = Qr * r_im_4;
                sigma_rdot_i = Qrdot * r_im_4;
                sigma_theta_i = Qtheta * r_im_4;
                dm(i) = ((r_hat(k, m, i) - Hy(Hl(i, l), 1)) / sigma_r_i)^2 + ((theta_hat(k, m, i) - Hy(Hl(i, l), 2)) / sigma_theta_i)^2 + ((rdot_hat(k, m, i) - Hy(Hl(i, l), 3)) / sigma_rdot_i)^2;
                em(i) = (exp(-dm(i) / 2)) / (sqrt((2*pi)^3) * sigma_r_i * sigma_theta_i * sigma_rdot_i);
                if Hy(Hl(i, l), :) ~= [0, 0, 0]
                    weight_l = weight_l * em(i);
                end
            end
            weight_pmk = weight_pmk + weight_l;
        end
       % in this part the hypo is now consedered wrongly, it needs to change again. 
    end