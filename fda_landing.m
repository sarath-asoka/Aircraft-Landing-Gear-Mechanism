function fda_landing()
% Forward dynamic analysis of Airplane Landing Gear
global Minv h

m1 = 2; m2 = 1; m3 = 1; m4 = 1; m5 = 1; g = 9.81;
l1 = 13.6; l2 = 5; l3 = 3; l4=3.1; l5=5; l12=3.151; l13=6.324; zhi12 = 1.064; zhi13= 0.08707;
J1 = m1*l1^2/12;J2 = m2*l2^2/12;J3 = m3*l3^2/12;J4 = m4*l4^2/12 ; J5 = m5*l5^2/12 ;
M = diag([m1 m1 J1 m2 m2 J2 m3 m3 J3 m4 m4 J4 m5 m5 J5]);
Minv = inv(M);
h = [0 -m1*g 0 0 -m2*g 0 0 -m3*g 0 0 -m4*g 0 0 -m5*g 0]';
Y0 = [-0.593 -6.774 1.396 -4.301 -2.766 5.1667 -0.947 -0.925 0.255 -3.898 -.912 6.027 -6.493 1.728 5.166 zeros(1,15)]';
options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[T,Y] = ode45(@eom,0:0.001:0.94,Y0,options);
err = zeros(size(T,1),1);
for i = 1:size(T,1)
    [Phi,D] = constraints(T(i),Y(i,1:15)');
    err(i)=norm(Phi);
end
plot(T,err);
end

function dy = eom(t,y)
global Minv h
[Phi,D] = constraints(t,y(1:15,1));
alpha = 10;beta = 10;
lambda = (D*Minv*D')\(gamma(y(1:15,1),y(16:30,1))-D*Minv*h-2*alpha*D*y(16:30,1)-beta^2*Phi);
dy(1:15,1) = y(16:30,1);
dy(16:30,1) = Minv*(h+D'*lambda);
norm(Phi)
end

function output = gamma(q,qdot)
l1 = 13.6;l2 = 5;l3 = 3;l4=3.1;l5=5;l12=3.151;l13=6.324;zhi12 = 1.064;zhi13= 0.08707;

theta1 = q(3);theta2 = q(6);theta3 = q(9);theta4 = q(12);theta5 = q(15);
theta1dot = qdot(3);theta2dot = qdot(6);theta3dot = qdot(9);theta4dot = qdot(12);theta5dot = qdot(15);

output = [ (-l5*cos(theta5)*theta5dot^2)/2;
    (-l5*sin(theta5)*theta5dot^2)/2 ;
    (l2*cos(theta2)*theta2dot^2)/2 + (l5*cos(theta5)*theta5dot^2)/2 ;
    (l2*sin(theta2)*theta2dot^2)/2 + (l5*sin(theta5)*theta5dot^2)/2 ;
    (-l12*cos(theta1 + zhi12)*theta1dot^2) + (l2*cos(theta2)*theta2dot^2)/2
    (-l12*sin(theta1 + zhi12)*theta1dot^2) + (l2*sin(theta2)*theta2dot^2)/2
    (l1*cos(theta1)*theta1dot^2)/2;
    (l1*sin(theta1)*theta1dot^2)/2;
    (l4*cos(theta4)*theta4dot^2)/2 + (l5*cos(theta5)*theta5dot^2)/2 ;
    (l4*sin(theta4)*theta4dot^2)/2 + (l5*sin(theta5)*theta5dot^2)/2 ;
    (l3*cos(theta3)*theta3dot^2)/2 + (l4*cos(theta4)*theta4dot^2)/2 ;
    (l3*sin(theta3)*theta3dot^2)/2 + (l4*sin(theta4)*theta4dot^2)/2 ;
    (-l13*cos(theta1 - zhi13)*theta1dot^2) + (l3*cos(theta3)*theta3dot^2)/2 ;
    (-l13*sin(theta1 - zhi13)*theta1dot^2) + (l3*sin(theta3)*theta3dot^2)/2 ;
    ] ; 
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
   
