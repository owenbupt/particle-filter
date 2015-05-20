function FalseArr = false_alarm(R_min, R_max, false_mean)
    % notice the condition that false_num = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % change the false number to 0 to help debug
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate the false alarm
        % sensors sensible area
    false_num = poissrnd(false_mean); %poisson distrubition with 1 mean
    %false_num
    FalseArr = zeros(false_num, 3);
    for i = 1:false_num
        % consider the uniform function to change the 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        FalseArr(i, 1) = R_min + (R_max - R_min) * rand;
        FalseArr(i, 2) = -pi + 2 * pi * rand;
        FalseArr(i, 3) = 120 + (222 - 120) * rand;
        % the r_dot false is between target 0 and 1;
        % there may be some problems for the choices.
    end
end
