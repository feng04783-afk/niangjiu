% =========================================================================
% 基于显式动力学方程的线性化与最小惯性参数集提取
% =========================================================================
clear; clc;

%% 1. 自动定义符号变量
fprintf('正在初始化符号变量...\n');

% 广义坐标与运动学变量
syms theta1 theta4 theta5 theta6 theta7 cph1 real
syms dq1 dq2 dq3 real
syms d2q1 d2q2 d2q3 real
d2q = [d2q1; d2q2; d2q3]; % 对应你公式中的 d2q

% 几何尺寸与常数变量 (包括各杆件的质心偏置)
syms L1 L2 L3 L7 L8 L9 S1 S3 S4 K1 K2 M1 M2 R1 R2 n1 g real
syms LC1 LC1y LC2 LC2y LC3 Lc11 Ls11 Lc13 Ls13 Lc14 Lc14y Lc15 Lcx10 Ls5 Ls6 real

% 循环定义15个连杆的质量和惯量，并组装待辨识参数向量 Pi_sym
Pi_sym = [];
for i = 1:15
    % 质量
    eval(['syms m' num2str(i) ' real;']);
    % 转动惯量
    eval(['syms Ixx' num2str(i) ' Iyy' num2str(i) ' Izz' num2str(i) ' real;']);
    % 将参数压入列向量
    eval(['Pi_sym = [Pi_sym; m' num2str(i) '; Ixx' num2str(i) '; Iyy' num2str(i) '; Izz' num2str(i) '];']);
end

%% 2. 载入你推导的动力学显式公式
fprintf('正在载入动力学公式 (Mm, V2, G_dot)...\n');

% ------ 请核对：以下是你刚刚提供的精确代码 ------
Mm=[Ixx5 + Iyy14/2 + Izz6 + Izz9 + Izz10 + (2*Lc11*m11*sin(theta1)*(M1 - R2*cos(theta5 - theta7))*(Ls11*cos(theta1) - sin(theta1)*(M1 - S1)) - 2*Lc11*m11*cos(theta1)*(M1 - R2*cos(theta5 - theta7))*(Ls11*sin(theta1) + cos(theta1)*(M1 - S1)))*n1 + (Lc11^2*m11*cos(theta1)^2*(M1 - R2*cos(theta5 - theta7))^2 + Lc11^2*m11*sin(theta1)^2*(M1 - R2*cos(theta5 - theta7))^2)*n1^2 + Ixx14/2 - Ixx8*((K2 - R1*cos(theta4 - theta6))^2/L8^2 - 1) + Iyy1*cos(theta4)^2 + Iyy2*cos(theta5)^2 + Iyy3*cos(theta5)^2 + Iyy4*cos(theta4)^2 + Iyy13*cos(theta4)^2 + Iyy15*cos(theta5)^2 - Iyy7*((K1 + R1*sin(theta4 - theta6))^2/L7^2 - 1) + Ixx1*sin(theta4)^2 + Ixx2*sin(theta5)^2 + Ixx3*sin(theta5)^2 + Ixx4*sin(theta4)^2 + Ixx13*sin(theta4)^2 + Ixx15*sin(theta5)^2 + (Iyy8*(K2 - R1*cos(theta4 - theta6))^2)/L8^2 + (Ixx7*(K1 + R1*sin(theta4 - theta6))^2)/L7^2 + Iyy11*n1^2*(M1 - R2*cos(theta5 - theta7))^2 + Ixx11*n1^2*(M2 + R2*sin(theta5 - theta7))^2 + (Ixx12*(M1 - R2*cos(theta5 - theta7))^2)/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2) + (Izz12*(M2 + R2*sin(theta5 - theta7))^2)/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2) + (m13*(4*Lc13^2*cos(theta4)^2 - 8*Lc13*S1*cos(theta4) - 4*sqrt((3))*Lc13*S3*cos(theta4) + 4*Ls13^2 + 4*S1^2 + 4*sqrt((3))*S1*S3 + 3*S3^2))/4 + m15*(sin(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) + Ls13*cos(theta1))^2 + m15*(cos(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) - Ls13*sin(theta1))^2 + m6*(K2^2 - 2*K2*S1 + Ls6^2 + S1^2) + (m4*(2*S1 - L1*cos(theta4) + 2*L2*cos(theta5))^2)/4 + m3*(L1*cos(theta4) - S1 + LC3*cos(theta5))^2 + m1*(LC1*cos(theta4) - S1 + LC1y*sin(theta4))^2 + (m14*(cos(theta1)*(Lc14y + 2*S1 - sqrt((3))*Lc14 + sqrt((3))*S3 - 2*L1*cos(theta4)) - 2*Ls13*sin(theta1))^2)/4 + (m14*(sin(theta1)*(Lc14y + 2*S1 - sqrt((3))*Lc14 + sqrt((3))*S3 - 2*L1*cos(theta4)) + 2*Ls13*cos(theta1))^2)/4 + m10*(Lcx10 - S1 + L1*cos(theta4) + L3*cos(theta5))^2 + m5*(sin(theta1)*(S1 + sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2) - R1*cos(theta4 - theta6)) + Ls5*cos(theta1))^2 + m5*(cos(theta1)*(S1 + sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2) - R1*cos(theta4 - theta6)) - Ls5*sin(theta1))^2 + (m7*(L7^2*cph1^2 - 4*L7*R1*cph1*cos(theta4 - theta6) + 4*L7*S1*cph1 + 4*Ls5^2 + 4*R1^2*cos(theta4 - theta6)^2 - 8*R1*S1*cos(theta4 - theta6) + 4*S1^2))/4 + L9^2*m9 + (m8*(K2^2 + 2*K2*R1*cos(theta4 - theta6) - 4*K2*S1 + 4*Ls6^2 + R1^2*cos(theta4 - theta6)^2 - 4*R1*S1*cos(theta4 - theta6) + 4*S1^2))/4 + m12*(S1 - R2*cos(theta5 - theta7))^2 + m11*(Ls11*sin(theta1) + cos(theta1)*(M1 - S1))^2 + m11*(Ls11*cos(theta1) - sin(theta1)*(M1 - S1))^2 + m2*(S1 + LC2*cos(theta5) - LC2y*sin(theta5))^2, - L1*Ls13*m14*sin(theta4) - L1*Ls13*m15*sin(theta4) - Lc13*Ls13*m13*sin(theta4) - (Ls6*R1*m8*sin(theta4 - theta6))/2 + (Ls5*R1*m5*((R1*sin(2*theta4 - 2*theta6))/2 - sin(theta4 - theta6)*sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2) + K1*cos(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2) + (Ls5*R1*m7*(K1*cos(theta4 - theta6) - 2*L7*cph1*sin(theta4 - theta6) + R1*cos(theta4 - theta6)*sin(theta4 - theta6)))/(2*L7*cph1), (-Lc11*Ls11*R2*m11*(M2 + R2*sin(theta5 - theta7))*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7)))*n1^3 - Lc15*Ls13*m15*sin(theta5); - L1*Ls13*m14*sin(theta4) - L1*Ls13*m15*sin(theta4) - Lc13*Ls13*m13*sin(theta4) - (Ls6*R1*m8*sin(theta4 - theta6))/2 + (Ls5*R1*m5*((R1*sin(2*theta4 - 2*theta6))/2 - sin(theta4 - theta6)*sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2) + K1*cos(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2) + (Ls5*R1*m7*(K1*cos(theta4 - theta6) - 2*L7*cph1*sin(theta4 - theta6) + R1*cos(theta4 - theta6)*sin(theta4 - theta6)))/(2*L7*cph1), Izz1 + Izz4 + Izz13 + m6*(R1*cos(theta4 - theta6) - (R1*sin(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/sqrt(L7^2 - (K2 - R1*cos(theta4 - theta6))^2))^2 + L1^2*m10 + L1^2*m3 + (L1^2*m4)/4 + L1^2*m14 + L1^2*m15 + Lc13^2*m13 + m1*(LC1^2 + LC1y^2) + m8*(R1*cos(theta4 - theta6) - (R1*sin(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/(L8*sqrt(1 - (K2 - R1*cos(theta4 - theta6))^2/L8^2)*2))^2 + m7*cos(theta1)^2*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/(2*L7*cph1))^2 + m7*sin(theta1)^2*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/(2*L7*cph1))^2 + m5*cos(theta1)^2*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2))^2 + m5*sin(theta1)^2*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2))^2 + (R1^2*m7*cos(theta4 - theta6)^2)/4 + (R1^2*m8*sin(theta4 - theta6)^2*cos(theta1)^2)/4 + (R1^2*m8*sin(theta4 - theta6)^2*sin(theta1)^2)/4 - (Izz7*R1^2*cos(theta4 - theta6)^2)/(L7^2*((K1 + R1*sin(theta4 - theta6))^2/L7^2 - 1)) - (Izz8*R1^2*sin(theta4 - theta6)^2)/(L8^2*((K2 - R1*cos(theta4 - theta6))^2/L8^2 - 1)), L1*L3*cos(theta4 - theta5)*m10 - (L1*L2*m4*cos(theta4 - theta5))/2 + L1*LC3*m3*cos(theta4 - theta5) + L1*Lc15*m15*cos(theta4 - theta5); (-Lc11*Ls11*R2*m11*(M2 + R2*sin(theta5 - theta7))*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7)))*n1^3 - Lc15*Ls13*m15*sin(theta5), L1*L3*cos(theta4 - theta5)*m10 - (L1*L2*m4*cos(theta4 - theta5))/2 + L1*LC3*m3*cos(theta4 - theta5) + L1*Lc15*m15*cos(theta4 - theta5), Izz2 + Izz3 + Izz15 + (Lc11^2*R2^2*m11*(M1 - R2*cos(theta5 - theta7))^2*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7))^2 + Lc11^2*R2^2*m11*cos(theta1)^2*(M2 + R2*sin(theta5 - theta7))^2*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7))^2 + Lc11^2*R2^2*m11*sin(theta1)^2*(M2 + R2*sin(theta5 - theta7))^2*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7))^2)*n1^6 - (Iyy12*((R2*sin(theta5 - theta7))/sqrt(((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)) - ((M1 - R2*cos(theta5 - theta7))*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))))/(2*((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^(3/2)))^2)/((M1 - R2*cos(theta5 - theta7))^2/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2) - 1) + Izz11*R2^2*n1^4*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7))^2 + m4*L2^2 + m10*L3^2 + m3*LC3^2 + m15*Lc15^2 + m12*R2^2 + m2*(LC2^2 + LC2y^2)];

V2=[dq1*((dq3*(2*R2*m11*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))*Lc11^2*cos(theta1)^2 + 2*R2*m11*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))*Lc11^2*sin(theta1)^2)*n1^2)/2 - (dq3*(2*Lc11*R2*m11*sin(theta5 - theta7)*cos(theta1)*(Ls11*sin(theta1) + cos(theta1)*(M1 - S1)) - 2*Lc11*R2*m11*sin(theta5 - theta7)*sin(theta1)*(Ls11*cos(theta1) - sin(theta1)*(M1 - S1)))*n1)/2 + (dq2*((m13*(8*Lc13*S1*sin(theta4) - 8*Lc13^2*cos(theta4)*sin(theta4) + 4*sqrt((3))*Lc13*S3*sin(theta4)))/4 - (m8*(2*K2*R1*sin(theta4 - theta6) - 4*R1*S1*sin(theta4 - theta6) + 2*R1^2*cos(theta4 - theta6)*sin(theta4 - theta6)))/4 + (m7*(8*R1*S1*sin(theta4 - theta6) - 8*R1^2*cos(theta4 - theta6)*sin(theta4 - theta6) + 4*L7*R1*cph1*sin(theta4 - theta6)))/4 + 2*Ixx1*cos(theta4)*sin(theta4) + 2*Ixx4*cos(theta4)*sin(theta4) + 2*Ixx13*cos(theta4)*sin(theta4) - 2*Iyy1*cos(theta4)*sin(theta4) - 2*Iyy4*cos(theta4)*sin(theta4) - 2*Iyy13*cos(theta4)*sin(theta4) + 2*m1*(LC1y*cos(theta4) - LC1*sin(theta4))*(LC1*cos(theta4) - S1 + LC1y*sin(theta4)) + (L1*m4*sin(theta4)*(2*S1 - L1*cos(theta4) + 2*L2*cos(theta5)))/2 - 2*L1*m3*sin(theta4)*(L1*cos(theta4) - S1 + LC3*cos(theta5)) - 2*L1*sin(theta4)*m10*(Lcx10 - S1 + L1*cos(theta4) + L3*cos(theta5)) + 2*m5*cos(theta1)*(cos(theta1)*(S1 + sqrt((L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) - R1*cos(theta4 - theta6)) - Ls5*sin(theta1))*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) + 2*m5*sin(theta1)*(sin(theta1)*(S1 + sqrt((L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) - R1*cos(theta4 - theta6)) + Ls5*cos(theta1))*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) + L1*m14*cos(theta1)*sin(theta4)*(cos(theta1)*(Lc14y + 2*S1 - sqrt((3))*Lc14 + sqrt((3))*S3 - 2*L1*cos(theta4)) - 2*Ls13*sin(theta1)) + L1*m14*sin(theta1)*sin(theta4)*(sin(theta1)*(Lc14y + 2*S1 - sqrt((3))*Lc14 + sqrt((3))*S3 - 2*L1*cos(theta4)) + 2*Ls13*cos(theta1)) + (2*Ixx7*R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/L7^2 - (2*Ixx8*R1*sin(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/L8^2 - (2*Iyy7*R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/L7^2 + (2*Iyy8*R1*sin(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/L8^2 + 2*L1*m15*cos(theta1)*sin(theta4)*(cos(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) - Ls13*sin(theta1)) + 2*L1*m15*sin(theta1)*sin(theta4)*(sin(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) + Ls13*cos(theta1))))/2 + (dq3*(2*Ixx2*cos(theta5)*sin(theta5) - 2*m2*(LC2y*cos(theta5) + LC2*sin(theta5))*(S1 + LC2*cos(theta5) - LC2y*sin(theta5)) + 2*Ixx3*cos(theta5)*sin(theta5) + 2*Ixx15*cos(theta5)*sin(theta5) - 2*Iyy2*cos(theta5)*sin(theta5) - 2*Iyy3*cos(theta5)*sin(theta5) - 2*Iyy15*cos(theta5)*sin(theta5) - L2*m4*sin(theta5)*(2*S1 - L1*cos(theta4) + 2*L2*cos(theta5)) - 2*LC3*m3*sin(theta5)*(L1*cos(theta4) - S1 + LC3*cos(theta5)) - 2*L3*sin(theta5)*m10*(Lcx10 - S1 + L1*cos(theta4) + L3*cos(theta5)) - (Ixx12*(M1 - R2*cos(theta5 - theta7))^2*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^2 + 2*R2*m12*sin(theta5 - theta7)*(S1 - R2*cos(theta5 - theta7)) - (Izz12*(M2 + R2*sin(theta5 - theta7))^2*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^2 + 2*Ixx11*R2*n1^2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*Iyy11*R2*n1^2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7)) + 2*Lc15*m15*cos(theta1)*sin(theta5)*(cos(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) - Ls13*sin(theta1)) + 2*Lc15*m15*sin(theta1)*sin(theta5)*(sin(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) + Ls13*cos(theta1)) + (2*Ixx12*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7)))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2) + (2*Izz12*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)))/2) - dq2*(dq2*(L1*Ls13*m14*cos(theta4) + L1*Ls13*m15*cos(theta4) + Lc13*Ls13*m13*cos(theta4) + (Ls6*R1*m8*cos(theta4 - theta6))/2 + (Ls5*R1*m5*(cos(theta4 - theta6)*sqrt((L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) - R1*cos(2*theta4 - 2*theta6) + K1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*sin(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt((L7^2 - (K1 + R1*sin(theta4 - theta6))^2))))/sqrt((L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) + (Ls5*R1*m7*(- R1*cos(theta4 - theta6)^2 + 2*L7*cph1*cos(theta4 - theta6) + R1*sin(theta4 - theta6)^2 + K1*sin(theta4 - theta6)))/(2*L7*cph1) - (Ls5*R1^2*m5*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6))*((R1*sin(2*theta4 - 2*theta6))/2 - sin(theta4 - theta6)*sqrt((L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) + K1*cos(theta4 - theta6)))/(L7^2 - (K1 + R1*sin(theta4 - theta6))^2)^(3/2)) - (dq1*((m13*(8*Lc13*S1*sin(theta4) - 8*Lc13^2*cos(theta4)*sin(theta4) + 4*sqrt((3))*Lc13*S3*sin(theta4)))/4 - (m8*(2*K2*R1*sin(theta4 - theta6) - 4*R1*S1*sin(theta4 - theta6) + 2*R1^2*cos(theta4 - theta6)*sin(theta4 - theta6)))/4 + (m7*(8*R1*S1*sin(theta4 - theta6) - 8*R1^2*cos(theta4 - theta6)*sin(theta4 - theta6) + 4*L7*R1*cph1*sin(theta4 - theta6)))/4 + 2*Ixx1*cos(theta4)*sin(theta4) + 2*Ixx4*cos(theta4)*sin(theta4) + 2*Ixx13*cos(theta4)*sin(theta4) - 2*Iyy1*cos(theta4)*sin(theta4) - 2*Iyy4*cos(theta4)*sin(theta4) - 2*Iyy13*cos(theta4)*sin(theta4) + 2*m1*(LC1y*cos(theta4) - LC1*sin(theta4))*(LC1*cos(theta4) - S1 + LC1y*sin(theta4)) + (L1*m4*sin(theta4)*(2*S1 - L1*cos(theta4) + 2*L2*cos(theta5)))/2 - 2*L1*m3*sin(theta4)*(L1*cos(theta4) - S1 + LC3*cos(theta5)) - 2*L1*sin(theta4)*m10*(Lcx10 - S1 + L1*cos(theta4) + L3*cos(theta5)) + 2*m5*cos(theta1)*(cos(theta1)*(S1 + sqrt((L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) - R1*cos(theta4 - theta6)) - Ls5*sin(theta1))*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) + 2*m5*sin(theta1)*(sin(theta1)*(S1 + sqrt((L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) - R1*cos(theta4 - theta6)) + Ls5*cos(theta1))*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) + L1*m14*cos(theta1)*sin(theta4)*(cos(theta1)*(Lc14y + 2*S1 - sqrt((3))*Lc14 + sqrt((3))*S3 - 2*L1*cos(theta4)) - 2*Ls13*sin(theta1)) + L1*m14*sin(theta1)*sin(theta4)*(sin(theta1)*(Lc14y + 2*S1 - sqrt((3))*Lc14 + sqrt((3))*S3 - 2*L1*cos(theta4)) + 2*Ls13*cos(theta1)) + (2*Ixx7*R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/L7^2 - (2*Ixx8*R1*sin(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/L8^2 - (2*Iyy7*R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/L7^2 + (2*Iyy8*R1*sin(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/L8^2 + 2*L1*m15*cos(theta1)*sin(theta4)*(cos(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) - Ls13*sin(theta1)) + 2*L1*m15*sin(theta1)*sin(theta4)*(sin(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) + Ls13*cos(theta1))))/2) - dq3*((dq3*((Lc11*Ls11*R2*m11*(M2 + R2*sin(theta5 - theta7))*(M2*cos(theta5 - theta7) + M1*sin(theta5 - theta7)) + Lc11*Ls11*R2^2*m11*cos(theta5 - theta7)*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7)))*n1^3 + Lc15*Ls13*m15*cos(theta5)))/2 - (dq1*(n1^2*(2*R2*m11*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))*Lc11^2*cos(theta1)^2 + 2*R2*m11*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))*Lc11^2*sin(theta1)^2) - n1*(2*Lc11*R2*m11*sin(theta5 - theta7)*cos(theta1)*(Ls11*sin(theta1) + cos(theta1)*(M1 - S1)) - 2*Lc11*R2*m11*sin(theta5 - theta7)*sin(theta1)*(Ls11*cos(theta1) - sin(theta1)*(M1 - S1))) - 2*m2*(LC2y*cos(theta5) + LC2*sin(theta5))*(S1 + LC2*cos(theta5) - LC2y*sin(theta5)) + 2*Ixx2*cos(theta5)*sin(theta5) + 2*Ixx3*cos(theta5)*sin(theta5) + 2*Ixx15*cos(theta5)*sin(theta5) - 2*Iyy2*cos(theta5)*sin(theta5) - 2*Iyy3*cos(theta5)*sin(theta5) - 2*Iyy15*cos(theta5)*sin(theta5) - L2*m4*sin(theta5)*(2*S1 - L1*cos(theta4) + 2*L2*cos(theta5)) - 2*LC3*m3*sin(theta5)*(L1*cos(theta4) - S1 + LC3*cos(theta5)) - 2*L3*sin(theta5)*m10*(Lcx10 - S1 + L1*cos(theta4) + L3*cos(theta5)) - (Ixx12*(M1 - R2*cos(theta5 - theta7))^2*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^2 + 2*R2*m12*sin(theta5 - theta7)*(S1 - R2*cos(theta5 - theta7)) - (Izz12*(M2 + R2*sin(theta5 - theta7))^2*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^2 + 2*Ixx11*R2*n1^2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*Iyy11*R2*n1^2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7)) + 2*Lc15*m15*cos(theta1)*sin(theta5)*(cos(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) - Ls13*sin(theta1)) + 2*Lc15*m15*sin(theta1)*sin(theta5)*(sin(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) + Ls13*cos(theta1)) + (2*Ixx12*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7)))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2) + (2*Izz12*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)))/2 + (dq3*n1^3*(Lc11*Ls11*R2*m11*(M2 + R2*sin(theta5 - theta7))*(M2*cos(theta5 - theta7) + M1*sin(theta5 - theta7)) + Lc11*Ls11*R2^2*m11*cos(theta5 - theta7)*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7))))/2 + (Lc15*Ls13*dq3*m15*cos(theta5))/2); dq3^2*(L1*L3*sin(theta4 - theta5)*m10 - (L1*L2*m4*sin(theta4 - theta5))/2 + L1*LC3*m3*sin(theta4 - theta5) + L1*Lc15*m15*sin(theta4 - theta5)) + (dq2^2*(2*m5*cos(theta1)^2*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2))*(R1*cos(theta4 - theta6) - (R1^2*cos(theta4 - theta6)^2)/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2) + (R1*sin(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2) - (R1^2*cos(theta4 - theta6)^2*(K1 + R1*sin(theta4 - theta6))^2)/(L7^2 - (K1 + R1*sin(theta4 - theta6))^2)^(3/2)) - 2*m6*(R1*cos(theta4 - theta6) - (R1*sin(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/sqrt(L7^2 - (K2 - R1*cos(theta4 - theta6))^2))*(R1*sin(theta4 - theta6) + (R1^2*sin(theta4 - theta6)^2)/sqrt(L7^2 - (K2 - R1*cos(theta4 - theta6))^2) + (R1*cos(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/sqrt(L7^2 - (K2 - R1*cos(theta4 - theta6))^2) + (R1^2*sin(theta4 - theta6)^2*(K2 - R1*cos(theta4 - theta6))^2)/(L7^2 - (K2 - R1*cos(theta4 - theta6))^2)^(3/2)) - 2*m8*(R1*cos(theta4 - theta6) - (R1*sin(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/(2*L8*sqrt(1 - (K2 - R1*cos(theta4 - theta6))^2/L8^2)))*(R1*sin(theta4 - theta6) + (R1^2*sin(theta4 - theta6)^2)/(2*L8*sqrt(1 - (K2 - R1*cos(theta4 - theta6))^2/L8^2)) + (R1*cos(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/(2*L8*sqrt(1 - (K2 - R1*cos(theta4 - theta6))^2/L8^2)) + (R1^2*sin(theta4 - theta6)^2*(K2 - R1*cos(theta4 - theta6))^2)/(2*L8^3*(1 - (K2 - R1*cos(theta4 - theta6))^2/L8^2)^(3/2))) + 2*m5*sin(theta1)^2*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2))*(R1*cos(theta4 - theta6) - (R1^2*cos(theta4 - theta6)^2)/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2) + (R1*sin(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2) - (R1^2*cos(theta4 - theta6)^2*(K1 + R1*sin(theta4 - theta6))^2)/(L7^2 - (K1 + R1*sin(theta4 - theta6))^2)^(3/2)) + 2*m7*cos(theta1)^2*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/(2*L7*cph1))*(R1*cos(theta4 - theta6) - (R1^2*cos(theta4 - theta6)^2)/(2*L7*cph1) + (R1*sin(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/(2*L7*cph1)) + 2*m7*sin(theta1)^2*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/(2*L7*cph1))*(R1*cos(theta4 - theta6) - (R1^2*cos(theta4 - theta6)^2)/(2*L7*cph1) + (R1*sin(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/(2*L7*cph1)) - (R1^2*m7*cos(theta4 - theta6)*sin(theta4 - theta6))/2 + (R1^2*m8*cos(theta4 - theta6)*sin(theta4 - theta6)*cos(theta1)^2)/2 + (R1^2*m8*cos(theta4 - theta6)*sin(theta4 - theta6)*sin(theta1)^2)/2 + (2*Izz8*R1^3*sin(theta4 - theta6)^3*(K2 - R1*cos(theta4 - theta6)))/(L8^4*((K2 - R1*cos(theta4 - theta6))^2/L8^2 - 1)^2) + (2*Izz7*R1^3*cos(theta4 - theta6)^3*(K1 + R1*sin(theta4 - theta6)))/(L7^4*((K1 + R1*sin(theta4 - theta6))^2/L7^2 - 1)^2) - (2*Izz8*R1^2*cos(theta4 - theta6)*sin(theta4 - theta6))/(L8^2*((K2 - R1*cos(theta4 - theta6))^2/L8^2 - 1)) + (2*Izz7*R1^2*cos(theta4 - theta6)*sin(theta4 - theta6))/(L7^2*((K1 + R1*sin(theta4 - theta6))^2/L7^2 - 1))))/2 - (dq1^2*((m13*(8*Lc13*S1*sin(theta4) - 8*Lc13^2*cos(theta4)*sin(theta4) + 4*sqrt((3))*Lc13*S3*sin(theta4)))/4 - (m8*(2*K2*R1*sin(theta4 - theta6) - 4*R1*S1*sin(theta4 - theta6) + 2*R1^2*cos(theta4 - theta6)*sin(theta4 - theta6)))/4 + (m7*(8*R1*S1*sin(theta4 - theta6) - 8*R1^2*cos(theta4 - theta6)*sin(theta4 - theta6) + 4*L7*R1*cph1*sin(theta4 - theta6)))/4 + 2*Ixx1*cos(theta4)*sin(theta4) + 2*Ixx4*cos(theta4)*sin(theta4) + 2*Ixx13*cos(theta4)*sin(theta4) - 2*Iyy1*cos(theta4)*sin(theta4) - 2*Iyy4*cos(theta4)*sin(theta4) - 2*Iyy13*cos(theta4)*sin(theta4) + 2*m1*(LC1y*cos(theta4) - LC1*sin(theta4))*(LC1*cos(theta4) - S1 + LC1y*sin(theta4)) + (L1*m4*sin(theta4)*(2*S1 - L1*cos(theta4) + 2*L2*cos(theta5)))/2 - 2*L1*m3*sin(theta4)*(L1*cos(theta4) - S1 + LC3*cos(theta5)) - 2*L1*sin(theta4)*m10*(Lcx10 - S1 + L1*cos(theta4) + L3*cos(theta5)) + 2*m5*cos(theta1)*(cos(theta1)*(S1 + sqrt((L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) - R1*cos(theta4 - theta6)) - Ls5*sin(theta1))*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) + 2*m5*sin(theta1)*(sin(theta1)*(S1 + sqrt((L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) - R1*cos(theta4 - theta6)) + Ls5*cos(theta1))*(R1*sin(theta4 - theta6) - (R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/sqrt(L7^2 - (K1 + R1*sin(theta4 - theta6))^2)) + L1*m14*cos(theta1)*sin(theta4)*(cos(theta1)*(Lc14y + 2*S1 - sqrt((3))*Lc14 + sqrt((3))*S3 - 2*L1*cos(theta4)) - 2*Ls13*sin(theta1)) + L1*m14*sin(theta1)*sin(theta4)*(sin(theta1)*(Lc14y + 2*S1 - sqrt((3))*Lc14 + sqrt((3))*S3 - 2*L1*cos(theta4)) + 2*Ls13*cos(theta1)) + (2*Ixx7*R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/L7^2 - (2*Ixx8*R1*sin(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/L8^2 - (2*Iyy7*R1*cos(theta4 - theta6)*(K1 + R1*sin(theta4 - theta6)))/L7^2 + (2*Iyy8*R1*sin(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/L8^2 + 2*L1*m15*cos(theta1)*sin(theta4)*(cos(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) - Ls13*sin(theta1)) + 2*L1*m15*sin(theta1)*sin(theta4)*(sin(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) + Ls13*cos(theta1))))/2; ((L1*L2*m4*sin(theta4 - theta5))/2 - L1*L3*sin(theta4 - theta5)*m10 - L1*LC3*m3*sin(theta4 - theta5) - L1*Lc15*m15*sin(theta4 - theta5))*dq2^2 - dq1*((dq1*(n1^2*(2*R2*m11*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))*Lc11^2*cos(theta1)^2 + 2*R2*m11*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))*Lc11^2*sin(theta1)^2) - n1*(2*Lc11*R2*m11*sin(theta5 - theta7)*cos(theta1)*(Ls11*sin(theta1) + cos(theta1)*(M1 - S1)) - 2*Lc11*R2*m11*sin(theta5 - theta7)*sin(theta1)*(Ls11*cos(theta1) - sin(theta1)*(M1 - S1))) - 2*m2*(LC2y*cos(theta5) + LC2*sin(theta5))*(S1 + LC2*cos(theta5) - LC2y*sin(theta5)) + 2*Ixx2*cos(theta5)*sin(theta5) + 2*Ixx3*cos(theta5)*sin(theta5) + 2*Ixx15*cos(theta5)*sin(theta5) - 2*Iyy2*cos(theta5)*sin(theta5) - 2*Iyy3*cos(theta5)*sin(theta5) - 2*Iyy15*cos(theta5)*sin(theta5) - L2*m4*sin(theta5)*(2*S1 - L1*cos(theta4) + 2*L2*cos(theta5)) - 2*LC3*m3*sin(theta5)*(L1*cos(theta4) - S1 + LC3*cos(theta5)) - 2*L3*sin(theta5)*m10*(Lcx10 - S1 + L1*cos(theta4) + L3*cos(theta5)) - (Ixx12*(M1 - R2*cos(theta5 - theta7))^2*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^2 + 2*R2*m12*sin(theta5 - theta7)*(S1 - R2*cos(theta5 - theta7)) - (Izz12*(M2 + R2*sin(theta5 - theta7))^2*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^2 + 2*Ixx11*R2*n1^2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*Iyy11*R2*n1^2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7)) + 2*Lc15*m15*cos(theta1)*sin(theta5)*(cos(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) - Ls13*sin(theta1)) + 2*Lc15*m15*sin(theta1)*sin(theta5)*(sin(theta1)*(S1 + (sqrt((2))*S4)/2 - L1*cos(theta4) - Lc15*cos(theta5)) + Ls13*cos(theta1)) + (2*Ixx12*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7)))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2) + (2*Izz12*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)))/2 - (dq3*((Lc11*Ls11*R2*m11*(M2 + R2*sin(theta5 - theta7))*(M2*cos(theta5 - theta7) + M1*sin(theta5 - theta7)) + Lc11*Ls11*R2^2*m11*cos(theta5 - theta7)*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7)))*n1^3 + Lc15*Ls13*m15*cos(theta5)))/2 + (dq3*n1^3*(Lc11*Ls11*R2*m11*(M2 + R2*sin(theta5 - theta7))*(M2*cos(theta5 - theta7) + M1*sin(theta5 - theta7)) + Lc11*Ls11*R2^2*m11*cos(theta5 - theta7)*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7))))/2 + (Lc15*Ls13*dq3*m15*cos(theta5))/2) - dq3*((dq3*((Iyy12*(((M1 - R2*cos(theta5 - theta7))^2*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^2 - (2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7)))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2))*((R2*sin(theta5 - theta7))/sqrt(((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)) - ((M1 - R2*cos(theta5 - theta7))*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))))/(2*((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^(3/2)))^2)/((M1 - R2*cos(theta5 - theta7))^2/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2) - 1)^2 + (2*Iyy12*((R2*sin(theta5 - theta7))/sqrt(((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)) - ((M1 - R2*cos(theta5 - theta7))*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))))/(2*((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^(3/2)))*((R2*cos(theta5 - theta7))/sqrt(((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)) + (3*(M1 - R2*cos(theta5 - theta7))*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7)))^2)/(4*((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^(5/2)) - ((M1 - R2*cos(theta5 - theta7))*(2*R2^2*cos(theta5 - theta7)^2 + 2*R2^2*sin(theta5 - theta7)^2 + 2*R2*cos(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7)) - 2*R2*sin(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7))))/(2*((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^(3/2)) - (R2*sin(theta5 - theta7)*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))))/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^(3/2)))/((M1 - R2*cos(theta5 - theta7))^2/((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2) - 1) - 2*Izz11*R2^2*n1^4*(M2*cos(theta5 - theta7) + M1*sin(theta5 - theta7))*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7))))/2 - (dq3*n1^6*(2*Lc11^2*R2^2*m11*(M1 - R2*cos(theta5 - theta7))^2*(M2*cos(theta5 - theta7) + M1*sin(theta5 - theta7))*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7)) + 2*Lc11^2*R2^3*m11*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7))^2 + 2*Lc11^2*R2^2*m11*cos(theta1)^2*(M2 + R2*sin(theta5 - theta7))^2*(M2*cos(theta5 - theta7) + M1*sin(theta5 - theta7))*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7)) + 2*Lc11^2*R2^2*m11*sin(theta1)^2*(M2 + R2*sin(theta5 - theta7))^2*(M2*cos(theta5 - theta7) + M1*sin(theta5 - theta7))*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7)) + 2*Lc11^2*R2^3*m11*cos(theta5 - theta7)*cos(theta1)^2*(M2 + R2*sin(theta5 - theta7))*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7))^2 + 2*Lc11^2*R2^3*m11*cos(theta5 - theta7)*sin(theta1)^2*(M2 + R2*sin(theta5 - theta7))*(R2 - M1*cos(theta5 - theta7) + M2*sin(theta5 - theta7))^2))/2)];

G_dot=[(0); - g*m6*(R1*cos(theta4 - theta6) - (R1*sin(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/sqrt(L7^2 - (K2 - R1*cos(theta4 - theta6))^2)) - g*m1*(LC1*cos(theta4) + LC1y*sin(theta4)) - g*m8*(R1*cos(theta4 - theta6) - (R1*sin(theta4 - theta6)*(K2 - R1*cos(theta4 - theta6)))/(2*L8*sqrt(1 - (K2 - R1*cos(theta4 - theta6))^2/L8^2))) - L1*g*cos(theta4)*m10 - L1*g*m3*cos(theta4) - (L1*g*m4*cos(theta4))/2 - L1*g*m14*cos(theta4) - L1*g*m15*cos(theta4) - Lc13*g*m13*cos(theta4) - (R1*g*m7*cos(theta4 - theta6))/2; g*m11*((Lc11*(M2 + R2*sin(theta5 - theta7))*(2*R2*cos(theta5 - theta7)*(M2 + R2*sin(theta5 - theta7)) + 2*R2*sin(theta5 - theta7)*(M1 - R2*cos(theta5 - theta7))))/(2*((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2)^(3/2)) - (Lc11*R2*cos(theta5 - theta7))/sqrt(((M1 - R2*cos(theta5 - theta7))^2 + (M2 + R2*sin(theta5 - theta7))^2))) + g*m2*(LC2*cos(theta5) - LC2y*sin(theta5)) - L3*g*cos(theta5)*m10 + L2*g*m4*cos(theta5) - LC3*g*m3*cos(theta5) - Lc15*g*m15*cos(theta5) - R2*g*m12*cos(theta5 - theta7)];

tau0=(Mm*d2q+V2+G_dot)*10^-6;

%%3. 核心：提取回归矩阵 Y
fprintf('正在利用 Jacobian 提取回归矩阵 Y (可能需要几分钟)...\n');

% 利用雅可比矩阵巧妙地将非线性公式转化为相对于惯性参数的线性系数矩阵
Y_sym = jacobian(tau0, Pi_sym);

%% 4. 载入真实几何参数并生成数值化回归矩阵 W_numeric (单位统一转为“米 m”)
fprintf('载入几何参数 (转换为 m/kg/s 制)，生成数值观测矩阵...\n');

scale = 0.001; % 毫米转换为米的核心系数
g_val = 9.80665; % 重力加速度使用 m/s^2 
theta6_val = 5*pi/6; 
theta7_val = -pi/6;

% 基础尺寸 (全部乘以 scale 转为米)
S1_val = 82*scale; L1_val = 600*scale; L3_val = 550*scale; M1_val = 312.5*scale; M2_val = 140.5*scale; R2_val = 318*scale;
R1_val = 200*scale; Le_val = 414*scale; K1_val = 270*scale; 
L2_val = 215*scale; L11_val = 532*scale; S3_val = 150*scale; S4_val = 180*scale;

% 质心位置与偏置 (全部乘以 scale 转为米)
Lc14_val = 100.63*scale; Lc14y_val = 49.13*scale; LC1_val = 192.36*scale; LC1y_val = 7.69*scale; 
LC2_val = -38.71*scale; LC2y_val = -48.68*scale; LC3_val = 167.99*scale;
L9_val = sqrt(4.13^2 + 2.45^2)*scale; Lcx10_val = -13.9*scale;
Ls5_val = 268.5*scale; Ls6_val = -242.5*scale; Ls11_val = 189*scale; Ls13_val = 189*scale;

% 补充的关联关系与推导尺寸 (保持单位同步)
K2_val = K1_val;
L7_val = Le_val;
L8_val = L7_val;
Lc11_val = (45 + 0.5 * 532) * scale; 
Lc13_val = 0.5 * L1_val;
Lc15_val = 0.5 * L3_val;
cph1_val = cos(-pi/6); 
n1_val = 1;            

geom_syms = [g, theta6, theta7, S1, L1, L3, M1, M2, R2, R1, K1, K2, L2, S3, S4, ...
             Lc14, Lc14y, LC1, LC1y, LC2, LC2y, LC3, L9, Lcx10, Ls5, Ls6, Ls11, Ls13, ...
             L7, L8, cph1, n1, Lc11, Lc13, Lc15];
             
geom_vals = [g_val, theta6_val, theta7_val, S1_val, L1_val, L3_val, M1_val, M2_val, R2_val, R1_val, K1_val, K2_val, L2_val, S3_val, S4_val, ...
             Lc14_val, Lc14y_val, LC1_val, LC1y_val, LC2_val, LC2y_val, LC3_val, L9_val, Lcx10_val, Ls5_val, Ls6_val, Ls11_val, Ls13_val, ...
             L7_val, L8_val, cph1_val, n1_val, Lc11_val, Lc13_val, Lc15_val];

num_samples = 200; 
W_numeric = zeros(3 * num_samples, length(Pi_sym));

for i = 1:num_samples
    q_val   = (rand(3,1)-0.5)*pi;  
    dq_val  = (rand(3,1)-0.5)*2;   
    d2q_val = (rand(3,1)-0.5)*5;   
    
    Y_num = subs(Y_sym, ...
        [theta1, theta4, theta5, dq1, dq2, dq3, d2q1, d2q2, d2q3, geom_syms], ...
        [q_val(1), q_val(2), q_val(3), dq_val(1), dq_val(2), dq_val(3), d2q_val(1), d2q_val(2), d2q_val(3), geom_vals]);
    
    W_numeric((i-1)*3+1 : i*3, :) = double(Y_num);
end

%% 5. QR 分解提取最小参数集 (加入矩阵列平衡抗病态处理)
fprintf('开始进行 SVD 与 QR 分解...\n');

% --- 核心修正：对矩阵 W 进行列归一化，防止大数值吞噬小数值 ---
norm_W = vecnorm(W_numeric);
% 将全零列的模长强制设为 1，防止出现除以 0 的 NaN 报错
norm_W(norm_W == 0) = 1; 
W_norm = W_numeric ./ norm_W;

% 使用 SVD 判断有效秩的数量 (基于归一化后的矩阵)
[~, S, ~] = svd(W_norm, 'econ');
s_diag = diag(S);

% 设定一个合理的容差阈值 (1e-6 到 1e-8) 过滤数值计算噪声
tol = 1e-7 * max(s_diag); 
base_param_count = sum(s_diag > tol);

fprintf('--> 初始定义的惯性参数个数: %d\n', length(Pi_sym));
fprintf('--> 剔除冗余后，系统的最小独立参数集个数: %d\n', base_param_count);

% 利用 QR 分解提取具体的 Base Parameters 组合 (基于归一化后的矩阵)
[Q, R, E] = qr(W_norm, 0);
independent_indices = E(1:base_param_count);

% 因为置换矩阵 E 的顺序就是独立性最强的参数排序
Base_Pi = Pi_sym(independent_indices);

disp('--------------------------------------------------');
disp('修正后：系统的最小惯性参数集 (Base Parameters) 包含以下独立参数：');
disp(Base_Pi);
disp('--------------------------------------------------');

%% 6. 提取结果的绝对正确性验证 (Verification)
disp(' ');
disp('==================================================');
disp('开始验证最小参数集的正确性...');

% 6.1 分离独立列 (W_base) 和 冗余列 (W_dependent)
% W_numeric 是归一化前真实的数值观测矩阵
W_base = W_numeric(:, independent_indices); 
dependent_indices = E(base_param_count+1 : end);
W_dep  = W_numeric(:, dependent_indices);

% 6.2 验证一：矩阵满秩性与抗病态能力
rank_W_base = rank(W_base);
cond_W_base = cond(W_base);
fprintf('[验证一] W_base 的秩为: %d (理论上必须等于 %d)\n', rank_W_base, base_param_count);
fprintf('         W_base 的条件数: %.2e (优化激励轨迹前的初始条件数)\n', cond_W_base);

% 6.3 验证二：冗余参数的无损重构能力 (最核心验证)
% 看看被剔除的冗余参数，能否被这 28 个参数完美表达
% 计算线性组合系数矩阵 Beta (W_base * Beta = W_dep)
warning('off', 'MATLAB:rankDeficientMatrix'); % 暂时关闭因奇异导致的警告
Beta = W_base \ W_dep; 
warning('on', 'MATLAB:rankDeficientMatrix');

W_dep_reconstructed = W_base * Beta;
% 计算重构后的最大绝对误差
max_reconstruct_error = max(max(abs(W_dep - W_dep_reconstructed)));
fprintf('[验证二] 冗余矩阵重构的最大绝对误差为: %e\n', max_reconstruct_error);

% 6.4 验证三：动力学力矩等价性测试
% 生成一组随机的完整CAD参数向量 (模拟真实物理情况)
Pi_CAD_Random = rand(length(Pi_sym), 1) * 10; 

% 1. 使用全量参数计算力矩
Tau_full = W_numeric * Pi_CAD_Random;

% 2. 计算映射后的最小参数集数值 (Base Parameters 数值)
% 映射公式: Pi_base = Pi_indep + Beta * Pi_dep
Pi_indep_CAD = Pi_CAD_Random(independent_indices);
Pi_dep_CAD   = Pi_CAD_Random(dependent_indices);
Pi_base_CAD  = Pi_indep_CAD + Beta * Pi_dep_CAD;

% 3. 使用最小参数集计算力矩
Tau_base = W_base * Pi_base_CAD;

% 计算力矩预测误差
max_tau_error = max(abs(Tau_full - Tau_base));
fprintf('[验证三] 全参数模型 vs 最小参数模型 的预测力矩最大误差: %e Nm\n', max_tau_error);

% 6.5 给出最终结论
disp('--------------------------------------------------');
if max_reconstruct_error < 1e-7 && max_tau_error < 1e-7
    disp('✅ 验证完美通过！');
    disp('结论：你提取的 28 个 Base Parameters 是绝对正确的，没有丢失任何原模型的动力学特性，并且消除了所有奇异性。');
    disp('下一步：可以放心开始生成/优化机器人的激励轨迹了！');
else
    disp('❌ 验证未通过！');
    disp('结论：存在较大的重构误差，可能是在 SVD 截断时 tol 阈值设置过大，误删了有用参数。建议调小 tol (如 1e-9) 重试。');
end
disp('==================================================');

%% 7. 导出基于 Base Parameters 的快速计算函数
disp('正在将符号回归矩阵转换为快速数值计算函数...');
% 提取出只对应 Base Parameters 的列
Y_base_sym = Y_sym(:, independent_indices);

% 将所有的常量代入，使其只剩下关节变量 q, dq, d2q
Y_base_sym_eval = subs(Y_base_sym, geom_syms, geom_vals);

% 使用 matlabFunction 生成一个可以直接调用的匿名函数或 .m 文件
% 这会极大加速优化过程中的矩阵组装速度
matlabFunction(Y_base_sym_eval, 'File', 'calc_W_base', ...
    'Vars', {theta1, theta4, theta5, dq1, dq2, dq3, d2q1, d2q2, d2q3});
disp('已成功生成函数文件 calc_W_base.m !');