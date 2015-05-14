function FalseAlarm = false_alarm()
    % notice the condition that false_num = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change the false number to 0 to help debug
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate the false alarm
        % sensors sensible area
    R_min = 25000;
    R_max = 100000;
    false_num = poissrnd(0); %poisson distrubition with 1 mean
    FalseAlarm = zeros(false_num, 3);
    for i = 1:false_num
        FalseAlarm(i, 1) = R_min + (R_max - R_min) * rand();
        FalseAlarm(i, 2) = -pi + 2 * pi * rand;
        FalseAlarm(i, 3) = 120 + (222 - 120) * rand;
        % the r_dot false is between target 0 and 1;
        % there may be some problems for the choices.
    end
end