% =========================================================================
% 机器人动力学参数辨识 - 有限傅里叶激励轨迹优化 (严格包含缸体物理限位)
% =========================================================================
clear; clc;

%% 1. 傅里叶级数基础设置
num_joints = 3;       % 3个活动关节: theta1, theta4, theta5
N = 5;                % 傅里叶级数谐波项数 (通常选 5)
f_base = 0.1;         % 轨迹基频 (Hz)，周期 T = 1/f_base = 10秒
wf = 2 * pi * f_base; % 基频角速度

% 优化变量 x 的构成: 每个关节 (q0, a1~a5, b1~b5) 共 11 个参数，总计 33 个变量
total_vars = num_joints * (2 * N + 1);

% 初始化猜测值 (给一个较小的随机扰动，防止陷入初始死区)
x0 = randn(total_vars, 1) * 0.5; 
% 为了让初值更合理，把基准位置 q0 设在工作空间中间附近
x0(1)  = pi/4;  % theta1 初始位置
x0(12) = pi/3;  % theta4 初始位置
x0(23) = pi/3;  % theta5 初始位置

%% 2. 物理限位约束定义 (从 PLC 代码推导)
% 速度限制 (rad/s) 和 加速度限制 (rad/s^2)
dq_max = [1.0; 1.0; 1.0];  
d2q_max = [3.0; 3.0; 3.0]; 

%% 3. 设置优化器 fmincon
disp('==================================================');
disp('开始进行傅里叶激励轨迹优化...');
disp('正在寻找最小条件数，同时严格保证推杆不越界...');
disp('由于需要搜索33维空间，可能需要十分钟左右，请耐心等待...');
disp('==================================================');

options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...          % SQP 在带强非线性约束时表现优秀
    'MaxFunctionEvaluations', 100000, ...
    'MaxIterations', 1000, ...
    'StepTolerance', 1e-5);

% 定义目标函数和非线性约束
objective_func = @(x) cost_function(x, num_joints, N, wf);
nonlcon_func   = @(x) trajectory_constraints(x, num_joints, N, wf, dq_max, d2q_max);

% 运行优化
[x_opt, fval] = fmincon(objective_func, x0, [], [], [], [], [], [], nonlcon_func, options);

fprintf('\n✅ 优化完成！观测矩阵的最小条件数为: %.4f\n', fval);

% 提取优化后的系数并保存
save('Optimized_Fourier_Coeffs.mat', 'x_opt', 'N', 'wf', 'num_joints');
disp('优化系数已保存至 Optimized_Fourier_Coeffs.mat');

%% =========================================================================
% 【目标函数】：计算整条轨迹下 W_base 矩阵的归一化条件数
function cost = cost_function(x, num_joints, N, wf)
    num_samples = 100; % 在一个周期内采样 100 个点
    t_array = linspace(0, 2*pi/wf, num_samples);
    
    % 预分配观测矩阵 (3 个方程 * 样本数, 28 个 Base Parameters)
    W_traj = zeros(num_joints * num_samples, 28); 
    
    for i = 1:num_samples
        t = t_array(i);
        [q, dq, d2q] = calc_fourier_kinematics(x, t, num_joints, N, wf);
        
        % 调用上一步自动生成的快速数值计算函数
        W_t = calc_W_base(q(1), q(2), q(3), dq(1), dq(2), dq(3), d2q(1), d2q(2), d2q(3));
        W_traj((i-1)*3+1 : i*3, :) = W_t;
    end
    
    % --- 核心修改：按列归一化，消除质量与惯量单位量级的巨大差异 ---
    norm_W = vecnorm(W_traj);
    norm_W(norm_W == 0) = 1; % 防止除以 0
    W_traj_norm = W_traj ./ norm_W;
    
    % 计算归一化后的条件数
    cost = cond(W_traj_norm);
    
    % 惩罚项：如果条件数算出 NaN 或 Inf (通常因为矩阵严重降秩)，给予极大惩罚
    if isnan(cost) || isinf(cost)
        cost = 1e8;
    end
end

%% =========================================================================
% 【非线性约束】：严格限制关节角度、速度、加速度，以及 底层推杆的长度
function [c, ceq] = trajectory_constraints(x, num_joints, N, wf, dq_max, d2q_max)
    num_samples = 80; % 密集检查约束
    t_array = linspace(0, 2*pi/wf, num_samples);
    
    % 预分配空间(虽然不是必须，但用固定模式追加更稳妥)
    c = [];
    ceq = [];
    
    % 载入缸体机构参数 (用于计算 t2, t3)
    Le = 414; K2 = 270; R1 = 200; 
    M1 = 312.5; M2 = 140.5; R2 = 318;
    
    for i = 1:num_samples
        t = t_array(i);
        [q, dq, d2q] = calc_fourier_kinematics(x, t, num_joints, N, wf);
        
        theta1 = q(1); theta4 = q(2); theta5 = q(3);
        
        % 1. 关节 1 的 PLC 角度软限位 [-pi/12, 7*pi/12]
        c = [c; theta1 - (pi/2 + pi/12)]; % theta1 < 上限 (c <= 0)
        c = [c; -pi/12 - theta1];         % theta1 > 下限
        
        % 2. 运动学映射：计算推杆 2 和 3 的长度
        theta2 = 5*pi/6 - theta4;
        theta3 = -pi/6 - theta5;
        
        % ======= 核心修改区 =======
        % 推杆 2 长度 t2 
        t2_inner = Le^2 - (K2 - R1*cos(theta2))^2;
        
        % 约束 1: 强制要求根号内部必须大于等于0 (即 -t2_inner <= 0)
        c = [c; -t2_inner]; 
        
        % 使用 max() 保证即使搜索过程违规，也不会产生复数导致优化器崩溃
        t2 = -R1*sin(theta2) + sqrt(max(t2_inner, 0));
        
        % 约束 2 & 3: 推杆 2 的 PLC 行程软限位 [110, 395]
        c = [c; t2 - 395]; 
        c = [c; 110 - t2];
        % =========================
        
        % 推杆 3 长度 t3 (内部为平方和，必定大于0，无需 max 保护)
        t3 = sqrt((M1 - R2*cos(theta3))^2 + (M2 - R2*sin(theta3))^2);
        % 推杆 3 的 PLC 行程软限位 [208, 492]
        c = [c; t3 - 492];
        c = [c; 208 - t3];
        
        % 3. 速度和加速度约束
        c = [c; abs(dq) - dq_max];
        c = [c; abs(d2q) - d2q_max];
    end
end

%% =========================================================================
% 【辅助函数】：根据傅里叶系数 x 和时间 t 解析计算位置、速度、加速度
function [q, dq, d2q] = calc_fourier_kinematics(x, t, num_joints, N, wf)
    q = zeros(num_joints, 1);
    dq = zeros(num_joints, 1);
    d2q = zeros(num_joints, 1);
    
    idx = 1;
    for i = 1:num_joints
        q0 = x(idx); idx = idx + 1;
        a  = x(idx : idx+N-1); idx = idx + N;
        b  = x(idx : idx+N-1); idx = idx + N;
        
        q_i = q0; dq_i = 0; d2q_i = 0;
        
        for k = 1:N
            q_i   = q_i   + (a(k)/(wf*k)) * sin(wf*k*t) - (b(k)/(wf*k)) * cos(wf*k*t);
            dq_i  = dq_i  + a(k) * cos(wf*k*t) + b(k) * sin(wf*k*t);
            d2q_i = d2q_i - a(k) * wf * k * sin(wf*k*t) + b(k) * wf * k * cos(wf*k*t);
        end
        
        q(i) = q_i; dq(i) = dq_i; d2q(i) = d2q_i;
    end
end