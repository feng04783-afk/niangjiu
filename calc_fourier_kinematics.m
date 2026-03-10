function [q, dq, d2q] = calc_fourier_kinematics(x, t, num_joints, N, wf)
% 根据傅里叶系数 x 和时间 t 解析计算位置、速度、加速度
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
        
        q(i) = q_i; 
        dq(i) = dq_i; 
        d2q(i) = d2q_i;
    end
end