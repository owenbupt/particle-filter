function O = observe(x_target, y_target, vx_target, vy_target, x_sen, y_sen, i, m, k, R_min, R_max)
    % observations from target i and sensor m at time k.
    if sqrt((x_target(k,i) - x_sen(m))^2 + (y_target(k,i) - y_sen(m))^2) > R_min  && sqrt((x_target(k,i) - x_sen(m))^2 + (y_target(k,i) - y_sen(m))^2) < R_max
        r = sqrt((x_target(k,i) - x_sen(m))^2 + (y_target(k,i) - y_sen(m))^2) + randn;
        theta = atan((y_target(k,i) - y_sen(m)) / (x_target(k,i) - y_sen(m))) + randn;
        % the result is radians
        r_dot = ((x_target(k,i) - x_sen(m)) * vx_target(k, i) + (y_target(k, i) - y_sen(m)) * vy_target(k, i)) / sqrt((x_target(k, i) - x_sen(m))^2 + (y_target(k, i) - y_sen(m))^2) + randn;
        O = [r theta r_dot];
        disp('observations:');
        disp(O);
    else
        O = [0 0 0];
        disp(sqrt((x_target(k,i) - x_sen(m))^2 + (y_target(k,i) - y_sen(m))^2));
        disp('No observations');        
    end
end