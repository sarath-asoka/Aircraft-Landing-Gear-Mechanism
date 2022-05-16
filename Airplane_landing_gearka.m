function Airplane_landing_gearka
%--------------------------------------------------------------------------
% Kinematic analysis of Airplane Nose landing gear mechanism in body-coordinates
% Position, velocity, and acceleration analysis from t = 0 to t = 5s
% Coordinate-partitioning method
%-------------------------------------------------------------------------- 
clc
%position analysis
q=zeros(15,1);t=0:0.001:0.94;N=size(t,2);
pcoordsall = zeros(15,N);
%% Initial guess solution
q(3) = 1.396;q(1) = -.593;q(2) = -6.774;q(4) = -4.301;q(5) = -2.766;q(6) = 5.1667;q(7) = -.947;q(8) = -.925;
q(9) = 0.255;q(10) = -3.898;q(11) = -.912;q(12) = 6.027;q(13) = -6.493;q(14) = 1.728;
%% Kinematic analysis
for i=1:N
    [phi1,phi1dot,phi1ddot] = driver(t(i));
    q(15)  = phi1;
    pcoordsall(:,i) = Newton_Raphson(q,15,1e-4,100,@constraints,t(i));
    fprintf('\n got q at t(i) %d\n',i); 
end
plot(t,pcoordsall(1,:))
hold on
plot(t,pcoordsall(13,:))

% Velocity analysis
qdot = zeros(15,1);vcoordsall = zeros(15,N);
for i = 1:N
    [phi1,phi1dot,phi1ddot] = driver(t(i));
    qdot(15) = phi1dot;
    [Phi,D] = constraints(t(i),pcoordsall(:,i));
    Dnew = D(:,1:14);
    rhs = -D(:,15)*qdot(15);
    qdotnew = Dnew\rhs;
    vcoordsall(:,i) = [qdotnew(1:14,1)' qdot(15)]';
end
plot(t,vcoordsall(7,:))
legend('vel')

% Acceleration Analysis
qddot = zeros(15,1);acoordsall = zeros(15,N);
for i = 1:N
    [phi1,phi1dot,phi1ddot] = driver(t(i));
    qddot(15) = phi1ddot;
    [Phi,D] = constraints(t(i),pcoordsall(:,i));
    Dnew = D(:,1:14);
    rhsa = -D(:,15)*qddot(3)+gamma(pcoordsall(:,i),vcoordsall(:,i));
    qddotnew = Dnew\rhsa;
    acoordsall(:,i) = [qddotnew(1:14,1)' qddot(15)]';
end
plot(t,acoordsall(7,:))
legend('ac')
save nlgdata.mat t pcoordsall vcoordsall acoordsall

end


function [phi1,phi1dot,phi1ddot] = driver(t)
if t<=2
    phi1 = 5.1667+t^3-t^4/2;
    phi1dot = 3*t^2-2*t^3;
    phi1ddot = 6*t-6*t^2;

end
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

function sol = Newton_Raphson(q,n,tol,iter_max,fun,t)
coord = zeros(iter_max,n+1);
flag = 0;sol = [];
for i = 1:iter_max
   [Phi,D] = fun(t,q);
   D = D(:,1:14);
%    fprintf('Determinant no. %d is %f',i,det(D)); 
   err = sqrt(Phi'*Phi);
   fprintf('\nerror no. %d is %4.5f',i,err); 
   coord(i,:) = [err q'];
   if err<tol
       flag = 1;
       sol = coord(i,2:n+1)';
       break;
   end
   delta_q = -D\Phi;
   delta_q = [delta_q(1:14,1)' 0]';
   q = q + delta_q;
end
if flag == 0
    disp('Convergence failed in Newton-Raphson');
    return;
end
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

