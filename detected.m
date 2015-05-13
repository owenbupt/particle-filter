    function Hy = detected(m, k)
        %%
         % Get the all heypothese
        FalseAlarm = false_alarm();
        [false_num, no_use] = size(FalseAlarm);
        Hy = zeros(T + false_num, 3);
        for i = 1:T
            detect_result = det_num();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if detect_result == 1
                Hy(i, :) = observe(x_target, y_target, vx_target, vy_target, i, m, k);
            else
                Hy(i, :) = zeros(1, 3);
            end
        end
        for i = 1:false_num
            Hy(T + i, :) = FalseAlarm(i);
        end
            function detection = det_num()
        % Calculate the number of detection in the detected probability
        u = rand;
        if u > Pd
            detection = 0;
        else
            detection = 1;
        end
    end
    end