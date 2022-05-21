function ida_landing
% Inverse Dynamic Analysis of Airplane Landing Gear

m1 = 2; m2 = 1; m3 = 1; m4 = 1; m5 = 1; g = 9.81;
l1 = 13.6; l2 = 5; l3 = 3; l4=3.1; l5=5; l12=3.151; l13=6.324; zhi12 = 1.064; zhi13= 0.08707;
J1 = m1*l1^2/12;J2 = m2*l2^2/12;J3 = m3*l3^2/12;J4 = m4*l4^2/12 ; J5 = m5*l5^2/12 ;
M = diag([m1 m1 J1 m2 m2 J2 m3 m3 J3 m4 m4 J4 m5 m5 J5]);
h = [0 -m1*g 0 0 -m2*g 0 0 -m3*g 0 0 -m4*g 0 0 -m5*g 0]';

load nlgdata
N = size(t,2);
torques = zeros(1,N);
for i = 1:N
    [Phi,D] = constraints(t(i),pcoordsall(:,i));
    D_driver = [zeros(1,14),1];
    D_all = [D' D_driver'];
    rhs = M*acoordsall(:,i)-h;
    lambda_all = D_all\rhs;
    h_driver = D_driver'*lambda_all(15);
    torques(i) = h_driver(15);
end
plot(t(1:940),torques(1,1:940))
end

function [Phi,D] = constraints(t,q)
l1 = 13.6;l2 = 5;l3 = 3;l4=3.1;l5=5;l12=3.151;l13=6.324;zhi12 = 1.064;zhi13= 0.08707;
 
x1 = q(1);y1 = q(2);theta1 = q(3);
x2 = q(4);y2 = q(5);theta2 = q(6);
x3 = q(7);y3 = q(8);theta3 = q(9);
x4 = q(10);y4=q(11);theta4 = q(12);
x5 = q(13);y5=q(14);theta5 = q(15);
Fx= -7.589 ; Fy = 3.975 ;


Phi = [x5 - 0.5*l5*cos(theta5) - Fx;
       y5 - 0.5*l5*sin(theta5) - Fy;
       x5 - x2 + 0.5*l5*cos(theta5)+ 0.5*l2*cos(theta2);
       y5 - y2 + 0.5*l5*sin(theta5)+ 0.5*l2*sin(theta2);
       x2 - x1 + 0.5*l2*cos(theta2)- l12*cos(theta1 + zhi12);
       y2 - y1 + 0.5*l2*sin(theta2)- l12*sin(theta1 + zhi12);
       x1 + 0.5*l1*cos(theta1);
       y1 + 0.5*l1*sin(theta1);
       x5 - x4 + 0.5*l5*cos(theta5)+ 0.5*l4*cos(theta4);
       y5 - y4 + 0.5*l5*sin(theta5)+ 0.5*l4*sin(theta4);
       x4 - x3 + 0.5*l4*cos(theta4)+ 0.5*l3*cos(theta3);
       y4 - y3 + 0.5*l4*sin(theta4)+ 0.5*l3*sin(theta3);
       x3 - x1 + 0.5*l3*cos(theta3)- l13*cos(theta1 - zhi13);
       y3 - y1 + 0.5*l3*sin(theta3)- l13*sin(theta1 - zhi13)];

D=  [ 0,  0,                        0,  0,  0,                   0,  0,  0,                   0,  0,  0,                   0, 1, 0,  (l5*sin(theta5))/2;
  0,  0,                        0,  0,  0,                   0,  0,  0,                   0,  0,  0,                   0, 0, 1, -(l5*cos(theta5))/2;
 0,  0,                        0, -1,  0, -(l2*sin(theta2))/2,  0,  0,                   0,  0,  0,                   0, 1, 0, -(l5*sin(theta5))/2;
  0,  0,                        0,  0, -1,  (l2*cos(theta2))/2,  0,  0,                   0,  0,  0,                   0, 0, 1,  (l5*cos(theta5))/2;
 -1,  0,  l12*sin(theta1 + zhi12),  1,  0, -(l2*sin(theta2))/2,  0,  0,                   0,  0,  0,                   0, 0, 0,                   0;
  0, -1, -l12*cos(theta1 + zhi12),  0,  1,  (l2*cos(theta2))/2,  0,  0,                   0,  0,  0,                   0, 0, 0,                   0;
  1,  0,      -(l1*sin(theta1))/2,  0,  0,                   0,  0,  0,                   0,  0,  0,                   0, 0, 0,                   0;
  0,  1,       (l1*cos(theta1))/2,  0,  0,                   0,  0,  0,                   0,  0,  0,                   0, 0, 0,                   0;
  0,  0,                        0,  0,  0,                   0,  0,  0,                   0, -1,  0, -(l4*sin(theta4))/2, 1, 0, -(l5*sin(theta5))/2;
  0,  0,                        0,  0,  0,                   0,  0,  0,                   0,  0, -1,  (l4*cos(theta4))/2, 0, 1,  (l5*cos(theta5))/2;
  0,  0,                        0,  0,  0,                   0, -1,  0, -(l3*sin(theta3))/2,  1,  0, -(l4*sin(theta4))/2, 0, 0,                   0;
  0,  0,                        0,  0,  0,                   0,  0, -1,  (l3*cos(theta3))/2,  0,  1,  (l4*cos(theta4))/2, 0, 0,                   0;
 -1,  0,  l13*sin(theta1 - zhi13),  0,  0,                   0,  1,  0, -(l3*sin(theta3))/2,  0,  0,                   0, 0, 0,                   0;
  0, -1, -l13*cos(theta1 - zhi13),  0,  0,                   0,  0,  1,  (l3*cos(theta3))/2,  0,  0,                   0, 0, 0,                   0];

   end