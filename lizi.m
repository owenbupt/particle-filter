x = 0.1; % initial state                      初始状�?

Q = 1; % process noise covariance                         过程噪声方差

R = 1; % measurement noise covariance                 测量噪声方差

tf = 50; % simulation length                                    模拟长度

N = 100; % number of particles in the particle filter                   颗粒过滤器中的粒子数


xhat = x;     %xhat=x=0.1

P = 2;

xhatPart = x;    %xhatPart=x=0.1


% Initialize the particle filter. 初始化粒子滤波，xpart值用来在不同时刻生成粒子

for i = 1 : N

    xpart(i) = x + sqrt(P) * randn;   %  randn产生标准正�?分布的随机数或矩阵的函数�?

end      %初始化xpart(i)为生成的100个随机粒�?

xArr = [x];    %xArr=x=0.1

xhatPartArr = [xhatPart];   %xhatPartArr = [xhatPart]=0.1


close all;


for k = 1 : tf %tf为时间长度，k可以理解为时间轴上的k时刻

% System simulation系统仿真

% x数据为时刻k的真实状态�?

    x = 0.5 * x + 25 * x / (1 + x^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn; %状�?方程(1)

    y = x^2 / 20 + sqrt(R) * randn;%观测方程(2)。观测方程是在观测�?和待估参数之间建立的函数关系式�?


    % Particle filter 生成100个粒子并根据预测和观测�?差�?计算各个粒子的权�?

    for i = 1 : N

        xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i)^2) + 8 * cos(1.2*(k-1)) + sqrt(Q) * randn;
        

        ypart = xpartminus(i)^2 / 20;

        vhat = y - ypart; %观测和预测的�?

        q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R); %根据差�?给出100个粒子对应的权重

    end


    % Normalize the likelihood of each a priori estimate.       每一个先验估计正常化的可能�?

    qsum = sum(q);
    
    qsum
    for i = 1 : N

        q(i) = q(i) / qsum;%归一化权重，较大的权重除以qsum不为零，大部分较小的权重除以qsum后为�?

    end




    % Resample.         重新取样

    for i = 1 : N

        u = rand; % uniform random number between 0 and 1      0�?之间的均�?��机数

        qtempsum = 0;

        for j = 1 : N

            qtempsum = qtempsum + q(j);

            if qtempsum >= u

                %重采样对低权重进行剔除，同时保留高权重，防止�?��的办�?

                xpart(i) = xpartminus(j);

                break;

            end

        end

    end


% The particle filter estimate is the mean of the particles.      粒子滤波的估计是颗粒的平均�?

xhatPart = mean(xpart); %经过粒子滤波处理后的均�?


% Plot the estimated pdf's at a specific time.      绘制在特定的时间估计的概率密度函�?






% Save data in arrays for later plotting

xArr = [xArr x];

xhatPartArr = [xhatPartArr xhatPart];

end


t = 0 : tf;

figure;

plot(t, xArr, 'b.', t, xhatPartArr, 'g'); 
%此图2对应xArr为真值，xhatPartArr为粒子滤波�?

xlabel('time step'); ylabel('state');

legend('True state', 'Particle filter estimate');