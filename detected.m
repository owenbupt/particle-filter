    function Hy = detected(x_target, y_target, vx_target, vy_target, x_sen, y_sen, k, m, T, Pd, R_min, R_max, FalseArr)
        %%
         % Get the all heypothese
        
        [false_num, no_use] = size(FalseArr);
        Hy = zeros(T + false_num, 3);
        for i = 1:T
            % Calculate the number of detection in the detected probability
            u = rand;
            if u > Pd
                detect_result = 0;
            else
                detect_result = 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if detect_result == 1
                Hy(i, :) = observe(x_target, y_target, vx_target, vy_target, x_sen, y_sen, i, m, k, R_min, R_max);
            else
                Hy(i, :) = zeros(1, 3);
                % disp('Pd no observe');
            end
        end

        for i = 1:false_num
            Hy(T + i, :) = FalseArr(i);
        end           
    
    end
