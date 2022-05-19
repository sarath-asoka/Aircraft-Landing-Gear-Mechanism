clear all;clc
myVideo = VideoWriter('myVideoFile'); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)
load nlgdata
l1 = 13.6;l2 = 5;l3 = 3;l4=3.1;l5=5;l12=3.151;l13=6.324;zhi12 = 1.064;zhi13= 0.08707;
figure
axis square
axis([-20 20 -20 20]);
line([0 0],[2 -15],'LineWidth',3,'color','k');


N = size(t,2);
drag = [];shock_strut =  [];link2 = [];link3 = [];link4 = [];Fpos =[];Apos = [];Bpos = [];Cpos = [];Dpos = [];Epos= [];timedisplay = [];

phi1 = pcoordsall(3,1);phi2 = pcoordsall(6,1);phi3 = pcoordsall(9,1);
phi4 = pcoordsall(12,1);phi5 = pcoordsall(15,1);
yF =-13.546 ;xF =-1.185 ;
xA = -3.205;yA= -5.013 ;
xA1 =-.46 ;yA1 = -5.253;
xB = 0.504 ;yB = -.546;
xB1 =-.044 ;yB1 = -.498;
xC = -2.398;yC = -1.305;
xD = -5.347;yD = -0.519;
xE = -7.589 ;yE = 3.975;
drag = line([xE xD],[yE yD],'linewidth',3);
shock_strut = line([0 xF],[0 yF],'linewidth',3);
link2 = line([xA xD],[yA yD],'linewidth',3);
link3 = line([xB xC],[yB yC],'linewidth',3);
link4 = line([xC xD],[yC yD],'linewidth',3);
link5 = line([xA1 xA],[yA1 yA],'linewidth',3);
link6 = line([xB1 xB],[yB1 yB],'linewidth',3);


Apos = rectangle('Position',[xA-0.1,yA-0.1,1,1],'Curvature',[1,1],'FaceColor','r');
Bpos = rectangle('Position',[xB-0.05,yB-0.05,1,1],'Curvature',[1,1],'FaceColor','r');
Cpos = rectangle('Position',[xC-0.05,yC-0.05,1,1],'Curvature',[1,1],'FaceColor','r');
Dpos = rectangle('Position',[xD-0.05,yD-0.05,0.5,0.5],'Curvature',[1,1],'FaceColor','r');
Epos = rectangle('Position',[xE-0.05,yE-0.05,0.5,0.5],'Curvature',[1,1],'FaceColor','r');
Fpos = rectangle('Position',[xF-0.05,yF-0.05,0.5,0.5],'Curvature',[1,1],'FaceColor','r');
timedisplay = text(1,1,num2str(t(1)));

for i = 1:1:N
set(drag,'xdata',[xE xD],'ydata',[yE yD]);  
set(shock_strut,'xdata',[0 xF],'ydata',[0 yF]);
set(link2,'xdata',[xA xD],'ydata',[yA yD]);
set(link3,'xdata',[xB xC],'ydata',[yB yC]);
set(link4,'xdata',[xC xD],'ydata',[yC yD]);
set(link5,'xdata',[xA1 xA],'ydata',[yA1 yA]);
set(link6,'xdata',[xB1 xB],'ydata',[yB1 yB]);
set(Apos,'Position',[xA-0.05,yA-0.05,0.1,0.1]);
set(Bpos,'Position',[xB-0.05,yB-0.05,0.1,0.1]);
set(Cpos,'Position',[xC-0.05,yC-0.05,0.1,0.1]);
set(Dpos,'Position',[xD-0.05,yD-0.05,0.1,0.1]);
set(Epos,'Position',[xE-0.05,yE-0.05,0.1,0.1]);
set(Fpos,'Position',[xF-0.05,yF-0.05,0.1,0.1]);
set(timedisplay,'Position',[1 1],'string',num2str(t(i)));


drawnow();
    phi1 = pcoordsall(3,i);
    phi2 = pcoordsall(6,i);
    phi3 = pcoordsall(9,i);
    phi4 = pcoordsall(12,i);
    phi5 = pcoordsall(15,i);
    x1 = pcoordsall(1,i);y1 = pcoordsall(2,i);
    x2 = pcoordsall(4,i);y2 = pcoordsall(5,i);
    x3 = pcoordsall(7,i);y3 = pcoordsall(8,i);
    r1 = [x1 y1]';r2 = [x2 y2]';r3 = [x3 y3]';
    s_1_B = [l13*cos(zhi13) -l13*sin(zhi13)]';s_1_F = [-l1/2 0]';s_1_B1 = [l13*cos(zhi13) 0]';
    s_1_A = [l12*cos(zhi12) l12*sin(zhi12)]';s_1_A1 = [l12*cos(zhi12) 0]';
    s_2_D = [-l2/2 0]';s_3_C = [-l3/2 0]';
    S_1_A =  Rot(phi1)*s_1_A;S_1_A1 =  Rot(phi1)*s_1_A1;
    S_1_B =  Rot(phi1)*s_1_B;S_1_B1 =  Rot(phi1)*s_1_B1;
    S_1_F =  Rot(phi1)*s_1_F;
    S_2_D =  Rot(phi2)*s_2_D;
    S_3_C =  Rot(phi3)*s_3_C;
    xA = x1+S_1_A(1,1);xA1 = x1 +S_1_A1(1,1); 
    yA = y1+S_1_A(2,1);yA1 = y1 +S_1_A1(2,1);
    xB = x1+S_1_B(1,1);xB1 = x1 +S_1_B1(1,1);
    yB = y1+S_1_B(2,1);yB1 = y1 +S_1_B1(2,1);
    xC = x3+S_3_C(1,1);
    yC = y3+S_3_C(2,1);
    xD = x2+S_2_D(1,1);
    yD = y2+S_2_D(2,1);
    xF = x1+S_1_F(1,1);
    yF = y1+S_1_F(2,1);
    xE = -7.589 ;yE = 3.975;
    pause(0.001)
    
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)

function output = Rot(phi)
output = [cos(phi) -sin(phi);sin(phi) cos(phi)];
end

