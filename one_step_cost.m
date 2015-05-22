function OneStepCost = one_step_cost(k, A_km, x_target_hat, y_target_hat, x_target, y_target)
    T = 2;
    M = 4;
    %%%%%%%%%%%%%%%%%to fix%%%%%%%%%%%%%%%%%%%%
    cost_sen0 = 0.5;
    cost_sen1 = 0.1;
    cost_sen2 = 0.1;
    cost_sen3 = 0.1;
    alpha = 0.5;
    % sensors start/stop cost
    cost_start = 0.1;
    c_usage = [cost_sen0, cost_sen1, cost_sen2, cost_sen3];
    % it is only used in the sensor management function
    error_cost = 0;
    sensor_cost = 0;
    for i = 1:T
        error_cost = error_cost + (x_target_hat(k,i) - x_target(k, i))^2 + (y_target_hat(k, i) - y_target(k, i))^2;
    end
    for m = 1:M
        sensor_cost = sensor_cost + c_usage(m) * A_km(k, m) + cost_start * abs(A_km(k, m) - A_km(k - 1, m));
        % the absolute value of difference between A(k) and A(k - 1) means the 
        % action changed
    end
    OneStepCost = alpha * error_cost / 10^6 + sensor_cost;
end