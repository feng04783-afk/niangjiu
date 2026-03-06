function [theta1,t1,t2,t3,theta21]=kine_p_2(dpt,delta_q2_all)%dpt：末端执行器的位姿坐标（x,y,z）
% delta_q2_all：参数修正量数组（用于补偿机械臂加工 / 装配误差）
%theta1,t1,t2,t3,theta21（关节角度、关键长度 / 位移参数）
% syms S1 S2 L1 L3;
X_offest = -25.4;
Y_offest = 50.520;
Z_offest = 0.042;
X_offest = 0;
Y_offest = 0;
Z_offest = 0;

delta_q2 = sum(delta_q2_all, 2);
% delta_q2 = - delta_q2; 
K1 = 270; 
K2 = 270 + delta_q2(17); 
R1 = 200 + delta_q2(18); 
R2 = 318 + delta_q2(21); 
L7 = 414 + delta_q2(19); 
L8 = 414 + delta_q2(19); 
M1 = 312.5 + delta_q2(22); 
M2 = 140.5;
%delta_q2 是误差修正量，给每个基础尺寸（如连杆长度、半径）叠加补偿值，提升运动学求解精度。


l_AB = 215 + delta_q2(11);
l_AD = 600 + delta_q2(5);
l_CD = 215;
l_BC = 600 + delta_q2(10);
l_AE = 150 + delta_q2(13);
l_DF = 150;
l_EF = 600 + delta_q2(12);
l_DK = 550 + delta_q2(7);
l_DG = 180 + delta_q2(15);
l_HK = 180;
l_GH = 550 + delta_q2(14);
%变量命名（l_AB/l_AD）符合机构学惯例，代表连杆 A-B、A-D 的长度。
thetaA = pi ;
thetaB = 5 * pi / 12 ;
thetaC = 5 * pi / 6 ;
thetaD = 5 * pi / 6 ;

a1=80.5037 + delta_q2(3);
d1=312.5177 + delta_q2(2);
L2=600 + delta_q2(5);
L3=550 + delta_q2(7);
a4=203;

x=dpt(1);y=dpt(2);z=dpt(3);

if x<0&&y<0
    theta1=atan(y/x);
elseif x<0&&y>0 
    theta1=atan(y/x)+2*pi;
elseif x>0
    theta1=atan(y/x)+pi;
elseif x==0&&y<0 
    theta1=0.5*pi;
elseif x==0&&y>0 
    theta1=1.5*pi;
elseif y==0&&x>0 
    theta1=pi;
elseif y==0&&x<0 
    theta1=0;
else
    theta1=0;
end
%  global theta1;
%根据末端坐标 (x,y) 计算基座的水平旋转角 theta1
if x~=0||(x==0&&y==0)
    %大臂转角
    theta2=acos(((z-d1-a4-X_offest)^2+((x-Z_offest*sin(theta1))/cos(theta1)-a1+Y_offest)^2+L2^2-L3^2)/(2*L2*((z-d1-a4-X_offest)^2+((x-Z_offest*sin(theta1))/cos(theta1)-a1+Y_offest)^2)^0.5))...
    -atan((z-d1-a4-X_offest)/((x-Z_offest*sin(theta1))/cos(theta1)-a1+Y_offest));
    %小臂转角
    theta5=acos(((z-d1-a4-X_offest)^2+((x-Z_offest*sin(theta1))/cos(theta1)-a1+Y_offest)^2+L3^2-L2^2)/(2*L3*((z-d1-a4-X_offest)^2+((x-Z_offest*sin(theta1))/cos(theta1)-a1+Y_offest)^2)^0.5))...
    +atan((z-d1-a4-X_offest)/((x-Z_offest*sin(theta1))/cos(theta1)-a1+Y_offest));
else
    %大臂转角
    theta2=acos(((z-d1-a4-X_offest)^2+((y+Z_offest*cos(theta1))/sin(theta1)-a1+Y_offest)^2+L2^2-L3^2)/(2*L2*((z-d1-a4-X_offest)^2+((y+Z_offest*cos(theta1))/sin(theta1)-a1+Y_offest)^2)^0.5))...
    -atan((z-d1-a4-X_offest)/((x-Z_offest*sin(theta1))/cos(theta1)-a1+Y_offest));
    %小臂转角
    theta5=acos(((z-d1-a4-X_offest)^2+((y+Z_offest*cos(theta1))/sin(theta1)-a1+Y_offest)^2+L3^2-L2^2)/(2*L3*((z-d1-a4-X_offest)^2+((y+Z_offest*cos(theta1))/sin(theta1)-a1+Y_offest)^2)^0.5))...
    +atan((z-d1-a4-X_offest)/((x-Z_offest*sin(theta1))/cos(theta1)-a1+Y_offest));
end
Q1=acos(((z-d1-a4-X_offest)^2+((x-Z_offest*sin(theta1))/cos(theta1)-a1+Y_offest)^2+L3^2-L2^2)/(2*L3*((z-d1-a4-X_offest)^2+((x-Z_offest*sin(theta1))/cos(theta1)-a1+Y_offest)^2)^0.5));
% Q2=atan(((x-Z_offest*sin(theta1))/cos(theta1)-a1+Y_offest)/(z-d1-a4-X_offest));
Q3=atan((z-d1-a4-X_offest)/((x-Z_offest*sin(theta1))/cos(theta1)-a1+Y_offest));
% q=[theta1,theta4,theta5].';
theta2 = pi - theta2;
theta5 = pi + theta5;
theta3=theta5-theta2;
theta4=2.5*pi-theta5;
ADC = theta3;
l_AC = sqrt(l_AD^2 + l_CD^2 - 2*l_AD*l_CD*cos(ADC));
DAC = acos((l_AD^2 + l_AC^2 - l_CD^2)/(2*l_AD*l_AC));
BAC = acos((l_AB^2 + l_AC^2 - l_BC^2)/(2*l_AB*l_AC));
BAD = DAC + BAC;
theta2_ = 2*pi + theta2 - BAD;
theta6 = theta2_ - 2*pi + thetaD;
theta5 = 2*pi - thetaC - theta2;
t1=R1*cos(theta5)+sqrt(L7^2-(K1-R1*sin(theta5))^2);
t2=-R1*sin(theta5)+sqrt(L8^2-(K2+R1*cos(theta5))^2);
t3=sqrt((M1+R2*cos(theta6))^2+(M2+R2*sin(theta6))^2);

DKH = theta4 - 5*pi/4;
l_DH = sqrt(l_DK^2 + l_HK^2 - 2*l_DK*l_HK*cos(DKH));
KDH = acos((l_DK^2 + l_DH^2 - l_HK^2)/(2*l_DK*l_DH));
GDH = acos((l_DG^2 + l_DH^2 - l_GH^2)/(2*l_DG*l_DH));
KDG = KDH + GDH;
FDC = KDG + thetaB - thetaA;
ADF = theta3 - FDC;
l_AF = sqrt(l_AD^2 + l_DF^2 - 2*l_AD*l_DF*cos(ADF));
DAF = acos((l_AD^2 + l_AF^2 - l_DF^2)/(2*l_AD*l_AF));
EAF = acos((l_AE^2 + l_AF^2 - l_EF^2)/(2*l_AE*l_AF));
EAD = DAF + EAF;
theta21 = 2*pi + theta2 - EAD;

theta1 = theta1 - delta_q2(1);
t2 = t2 - delta_q2(16);
t3 = t3 - delta_q2(20);
end

 