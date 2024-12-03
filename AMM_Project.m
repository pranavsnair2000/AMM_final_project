%% start
clc
clear
% close all

sim_coppelia=0;  % open 5_link_manipulator.ttt and 

%% Phase 1

l1=2; l2=2.5; l3=1; l4=1.5; l5=1;

T=20;
T_sim=T*1;
Tspan=linspace(0,T_sim,360);

% ellipse
a=3; b=1;
C_e = [1;1];
r_t = pi/12;
R = [ cos(r_t), -sin(r_t) ; sin(r_t), cos(r_t)];

% Square
s=4.5;
C_s=[-0.5;1];


% did not end up using pcrt for the whole project, but used it for initial visualization:

% pcrt plot
L1 = Link('d', 0, 'a', l1, 'alpha', 0);
L2 = Link('d', 0, 'a', l2, 'alpha', 0);
L3 = Link('d', 0, 'a', l3, 'alpha', 0);
L4 = Link('d', 0, 'a', l4, 'alpha', 0);
L5 = Link('d', 0, 'a', l5, 'alpha', 0);

bot = SerialLink([L1 L2 L3 L4 L5], 'name', '5 Link Manipulator');

theta=[pi/6,pi/6,pi/6,pi/6,pi/6];
% theta=pi/4 * ones(1,5);

% figure(1)
% bot.plot(theta)



%% Phase 2: ellipse
T_sim=T*3;
Tspan=linspace(0,T_sim,360*3);


% initial desired position for phi=0
phi=0;
Xd_init =  C_e + R*[a*cos(phi) ; b*sin(phi)];
td_init=fsolve_5R(Xd_init,[0,0,0,0,0]);

Y_init= [td_init';phi];


% i

function Y_dot=i_trad_opnlp_der_ellipse(t,y)
    
    [J,Xdot] = get_J_Xdot_opnlp_ellipse(t,y);
    
    w=2*pi/20;
    Jp=pinv(J);
    
    Td_dot = Jp*Xdot; % no homogenous part for traditional pseudoinverse solution
    
    Y_dot=[Td_dot;w];
end

% ODE
[tout,yout] = ode45(@i_trad_opnlp_der_ellipse, Tspan, Y_init);

T1=yout(:,1);
T2=yout(:,2);
T3=yout(:,3);
T4=yout(:,4);
T5=yout(:,5);

% plot_ellipse(1,'Phase 2 (i) traditional pseudoinverse solution: ellipse',T1,T2,T3,T4,T5)


% disp('CoppeliaSim: traditional pseudoinverse solution for ellipse')
% coppelia(td_init,[T1,T2,T3,T4,T5],Tspan)


% ===============================================================================================

% ii

% ii a
function Y_dot=ii_a_opnlp_der_ellipse(t,y)
    [J,Xdot] = get_J_Xdot_opnlp_ellipse(t,y);

    w=2*pi/20;
    Jp=pinv(J);
    
    % particular
    Td_dot_p = Jp*Xdot;

    t4dot = Td_dot_p(4);
    t5dot = Td_dot_p(5);

    Tz=[1,0,0,t4dot,t5dot]';
    
    % homogenous
    Td_dot_h = (eye(5) - Jp*J)*Tz;
    
    % solution
    Td_dot = Td_dot_p + Td_dot_h;

    Y_dot=[Td_dot;w];
end

% ODE
[tout,yout] = ode45(@ii_a_opnlp_der_ellipse, Tspan, Y_init);

T1=yout(:,1);
T2=yout(:,2);
T3=yout(:,3);
T4=yout(:,4);
T5=yout(:,5);

% plot_ellipse(1,'Phase 2 (ii) a : ellipse',T1,T2,T3,T4,T5)

% disp('CoppeliaSim: Phase 2 (ii) a : ellipse')
% coppelia(td_init,[T1,T2,T3,T4,T5],Tspan)


% ============================================================

% ii b
function Y_dot=ii_b_opnlp_der_ellipse(t,y)
    [J,Xdot] = get_J_Xdot_opnlp_ellipse(t,y);

    w=2*pi/20;
    Jp=pinv(J);
    
    % particular
    Td_dot_p = Jp*Xdot;
    
    t1dot = Td_dot_p(1);
    t4dot = Td_dot_p(4);
    t5dot = Td_dot_p(5);

    % Tz=[t1dot, -t1dot, 0, t4dot, t5dot]';
    Tz=[1, -1, 0, t4dot, t5dot]';
    
    % homogenous
    Td_dot_h = (eye(5) - Jp*J)*Tz;
    
    % solution
    Td_dot = Td_dot_p + Td_dot_h;

    Y_dot=[Td_dot;w];
end

% b
% ODE
[tout,yout] = ode45(@ii_b_opnlp_der_ellipse, Tspan, Y_init);

T1=yout(:,1);
T2=yout(:,2);
T3=yout(:,3);
T4=yout(:,4);
T5=yout(:,5);

% plot_ellipse(1,'Phase 2 (ii) b : ellipse',T1,T2,T3,T4,T5)


% disp('CoppeliaSim: Phase 2 (ii) b : ellipse')
% coppelia(td_init,[T1,T2,T3,T4,T5],Tspan)

% =========================================================================

% iii

function Y_dot=iii_art_potent_opnlp_der_ellipse(t,y)
    [J,Xdot] = get_J_Xdot_opnlp_ellipse(t,y);
    
    w=2*pi/20;
    Jp=pinv(J);
    
    % particular
    Td_dot_p = Jp*Xdot;
    
    % V = (t1-pi/6)^2 + 0.25(t2-pi/2)^2 + 0.66(t3-pi/3)^2
    t1=y(1); t2=y(2); t3=y(3);
    Vd = [2*(t1-pi/6), 0.5*(t2-pi/2), 1.32*(t3-pi/3)];
    
    t4dot = Td_dot_p(4);
    t5dot = Td_dot_p(5);
    
    Tz=[-Vd,t4dot,t5dot]';
    
    % homogenous
    Td_dot_h = (eye(5) - Jp*J)*Tz;
    
    % solution
    Td_dot = Td_dot_p + Td_dot_h;
    
    Y_dot=[Td_dot;w];
end

% ODE
[tout,yout] = ode45(@iii_art_potent_opnlp_der_ellipse, Tspan, Y_init);

T1=yout(:,1);
T2=yout(:,2);
T3=yout(:,3);
T4=yout(:,4);
T5=yout(:,5);

% plot_ellipse(1,'Phase 2 (iii) : ellipse',T1,T2,T3,T4,T5)

% disp('CoppeliaSim: Phase 2 (iii) : ellipse')
% coppelia(td_init,[T1,T2,T3,T4,T5],Tspan)


%% Phase 2: square

% start at the bottom right corner
Xd_init =  C_s + [s/2;-s/2];
td_init=fsolve_5R(Xd_init,[0,0,0,0,0]);

Y_init= [td_init'];


% i

function Y_dot=i_trad_opnlp_der_square(t,y)
    
    Xdot= get_Xdot_square(t);

    J = get_J(y);
    Jp=pinv(J);

    Td_dot = Jp*Xdot; % no homogenous part for traditional pseudoinverse solution
    Y_dot=[Td_dot];
end


% ODE
[tout,yout] = ode45(@i_trad_opnlp_der_square, Tspan, Y_init);


T1=yout(:,1);
T2=yout(:,2);
T3=yout(:,3);
T4=yout(:,4);
T5=yout(:,5);

% plot_square(1,'Phase 2 (i) traditional pseudoinverse solution: square',T1,T2,T3,T4,T5)

% disp('CoppeliaSim: traditional pseudoinverse solution for square')
% coppelia(td_init,[T1,T2,T3,T4,T5],Tspan)

% ===============================================================================================


% ii a

function Y_dot=ii_a_opnlp_der_square(t,y)
    % y = [t1, t2, t3, t4, t5]
    
    Xdot=get_Xdot_square(t);
    
    J = get_J(y);
    Jp=pinv(J);
    
    % particular
    Td_dot_p = Jp*Xdot;

    t4dot = Td_dot_p(4);
    t5dot = Td_dot_p(5);

    Tz=[1,0,0,t4dot,t5dot]';
    
    % homogenous
    Td_dot_h = (eye(5) - Jp*J)*Tz;
    
    % solution
    Td_dot = Td_dot_p + Td_dot_h;

    Y_dot=Td_dot;

end
% ===============================


% ODE
[tout,yout] = ode45(@ii_a_opnlp_der_square, Tspan, Y_init);

T1=yout(:,1);
T2=yout(:,2);
T3=yout(:,3);
T4=yout(:,4);
T5=yout(:,5);

% plot_square(1,'Phase 2 (ii) a : square',T1,T2,T3,T4,T5)


% disp('CoppeliaSim: Phase 2 (ii) a : square')
% coppelia(td_init,[T1,T2,T3,T4,T5],Tspan)


% ===============================================================================================


% ii b
function Y_dot=ii_b_opnlp_der_square(t,y)
    % y = [t1, t2, t3, t4, t5]

    Xdot = get_Xdot_square(t);
    
    J = get_J(y);
    Jp=pinv(J);
    
    % particular
    Td_dot_p = Jp*Xdot;

    t4dot = Td_dot_p(4);
    t5dot = Td_dot_p(5);
    
    Tz=[1,-1,0,t4dot,t5dot]';
    
    % homogenous
    Td_dot_h = (eye(5) - Jp*J)*Tz;
    
    % solution
    Td_dot = Td_dot_p + Td_dot_h;

    Y_dot=Td_dot;

end
% ===============================



% ODE
[tout,yout] = ode45(@ii_b_opnlp_der_square, Tspan, Y_init);

T1=yout(:,1);
T2=yout(:,2);
T3=yout(:,3);
T4=yout(:,4);
T5=yout(:,5);

% plot_square(1,'Phase 2 (ii) b : square',T1,T2,T3,T4,T5)


% disp('CoppeliaSim: Phase 2 (ii) b : square')
% coppelia(td_init,[T1,T2,T3,T4,T5],Tspan)


% ===============================================================================================



% ii b
% ===============================
function Y_dot=iii_art_potent_opnlp_der_square(t,y)
    
    Xdot = get_Xdot_square(t);
    
    J = get_J(y);
    Jp=pinv(J);
    
    % particular
    Td_dot_p = Jp*Xdot;
    
    % V = (t1-pi/6)^2 + 0.25(t2-pi/2)^2 + 0.66(t3-pi/3)^2
    t1=y(1); t2=y(2); t3=y(3);
    Vd = [2*(t1-pi/6), 0.5*(t2-pi/2), 1.32*(t3-pi/3)];
    
    t4dot = Td_dot_p(4);
    t5dot = Td_dot_p(5);
    
    Tz=[-Vd,t4dot,t5dot]';

    
    % homogenous
    Td_dot_h = (eye(5) - Jp*J)*Tz;
    
    % solution
    Td_dot = Td_dot_p + Td_dot_h;

    Y_dot=Td_dot;

end
% ===============================


% ODE
[tout,yout] = ode45(@iii_art_potent_opnlp_der_square, Tspan, Y_init);

T1=yout(:,1);
T2=yout(:,2);
T3=yout(:,3);
T4=yout(:,4);
T5=yout(:,5);

% plot_square(1,'Phase 2 (iii)  : square',T1,T2,T3,T4,T5)

% disp('CoppeliaSim: Phase 2 (ii) b : square')
% coppelia(td_init,[T1,T2,T3,T4,T5],Tspan)


% ===============================================================================================



%% Phase 3 ellipse
T_sim=T*3;
Tspan=linspace(0,T_sim,360*3);


% ii closed loop joint space control: ellipse
% initial desired position for phi=0
Xd_init =  C_e + R*[a*cos(phi) ; b*sin(phi)];
td_init=fsolve_5R(Xd_init,[0,0,0,0,0]);

% actual initial position
X_init=Xd_init;
% X_init=[3;2];
t_init=fsolve_5R(X_init,[0,0,0,0,0]);


Y_init= [t_init';td_init';phi];

% ODE
[tout,yout] = ode45(@der_clsdlp_jntspc_ellipse, Tspan, Y_init);

T1=yout(:,1);
T2=yout(:,2);
T3=yout(:,3);
T4=yout(:,4);
T5=yout(:,5);

% plot_ellipse(1,'Phase 3 joint space controller : ellipse',T1,T2,T3,T4,T5)

% disp('CoppeliaSim: Phase 3 joint space controller : ellipse')
% coppelia(t_init,[T1,T2,T3,T4,T5],Tspan)


% ==================================================================================

% ii closed loop task space control: ellipse

% initial desired position for phi=0
Xd_init =  C_e + R*[a*cos(phi) ; b*sin(phi)];

% td_init=fsolve_RR(X_init,[0,0],l1,l2);

% actual initial position
% X_init=Xd_init;
X_init=[1;4];
t_init=fsolve_5R(X_init,[0,0,0,0,0]);


Y_init= [X_init;Xd_init;t_init';phi];

[tout,yout] = ode45(@der_clsdlp_tskspc_ellipse, Tspan, Y_init);

T1=yout(:,5);
T2=yout(:,6);
T3=yout(:,7);
T4=yout(:,8);
T5=yout(:,9);

% plot_ellipse(1,'Phase 3 task space controller : ellipse',T1,T2,T3,T4,T5)

% disp('CoppeliaSim: Phase 3 task space controller : ellipse')
% coppelia(t_init,[T1,T2,T3,T4,T5],Tspan)



%% Phase 3: square
% closed loop joint space control: square


% start at the bottom right corner
Xd_init =  C_s + [s/2;-s/2];
td_init=fsolve_5R(Xd_init,[0,0,0,0,0]);

% Y_init= [td_init'];

% actual initial position
% X_init=Xd_init;
X_init=[2;-1];
t_init=fsolve_5R(X_init,[0,0,0,0,0]);

Y_init= [t_init';td_init'];


% ODE
[tout,yout] = ode45(@der_clsdlp_jntspc_square, Tspan, Y_init);


T1=yout(:,1);
T2=yout(:,2);
T3=yout(:,3);
T4=yout(:,4);
T5=yout(:,5);

% plot_square(1,'Phase 3 joint space controller : square act',T1,T2,T3,T4,T5)
% plot_square(12,'Phase 3 joint space controller : square des',yout(:,6),yout(:,7),yout(:,8),yout(:,9),yout(:,10))

% disp('CoppeliaSim: Phase 3 joint space controller : square')
% coppelia(t_init,[T1,T2,T3,T4,T5],Tspan)


% ===========================================================================================================

% closed loop task space control: square

% start at the bottom right corner
Xd_init =  C_s + [s/2;-s/2]
td_init=fsolve_5R(Xd_init,[0,0,0,0,0]);

% Y_init= [td_init'];

t_side=length(Tspan)/4; 

% actual initial position
% X_init=Xd_init;
X_init=[3;2];
t_init=fsolve_5R(X_init,[0,0,0,0,0]);

Y_init= [X_init;Xd_init;t_init'];


% ODE
[tout,yout] = ode45(@der_clsdlp_tskspc_square, Tspan, Y_init);

T1=yout(:,5);
T2=yout(:,6);
T3=yout(:,7);
T4=yout(:,8);
T5=yout(:,9);

% plot_square(1,'Phase 3 task space controller : square',T1,T2,T3,T4,T5)

% disp('CoppeliaSim: Phase 3 task space controller : square')
% coppelia(t_init,[T1,T2,T3,T4,T5],Tspan)

%% Intermediate functions

function [J,Xdot] = get_J_Xdot_opnlp_ellipse(t,y)
    % y = [t1, t2, t3, t4, t5, phi]
    a=3; b=1;
    % C_e = [1;1];
    r_t=pi/12;
    R = [ cos(r_t), -sin(r_t) ; sin(r_t), cos(r_t)];
    
    l1=2; l2=2.5; l3=1; l4=1.5; l5=1;
    w=2*pi/20;
    phi=y(6);
    
    % X_dot
    xdot= a * -sin(phi) * w;
    ydot= b * cos(phi) * w;
    Xdot = R*[xdot;ydot];
    
    J=get_J(y(1:5));

end 

% ==========================================================================================
function J = get_J(y)
    % y = [t1, t2, t3, t4, t5]

    l1=2; l2=2.5; l3=1; l4=1.5; l5=1;

    t1=y(1); t2=y(2); t3=y(3); t4=y(4); t5=y(5);
    
    abs_t1 = t1;
    abs_t2 = abs_t1 + t2; 
    abs_t3 = abs_t2 + t3; 
    abs_t4 = abs_t3 + t4; 
    abs_t5 = abs_t4 + t5;
    
    J_t1= [- l1*sin(abs_t1) - l2*sin(abs_t2) - l3*sin(abs_t3) - l4*sin(abs_t4) - l5*sin(abs_t5);
           + l1*cos(abs_t1) + l2*cos(abs_t2) + l3*cos(abs_t3) + l4*cos(abs_t4) + l5*cos(abs_t5)];
    J_t2= [- l2*sin(abs_t2) - l3*sin(abs_t3) - l4*sin(abs_t4) - l5*sin(abs_t5);
           + l2*cos(abs_t2) + l3*cos(abs_t3) + l4*cos(abs_t4) + l5*cos(abs_t5)];
    J_t3= [- l3*sin(abs_t3) - l4*sin(abs_t4) - l5*sin(abs_t5);
           + l3*cos(abs_t3) + l4*cos(abs_t4) + l5*cos(abs_t5)];
    J_t4= [- l4*sin(abs_t4) - l5*sin(abs_t5);
           + l4*cos(abs_t4) + l5*cos(abs_t5)];
    J_t5= [- l5*sin(abs_t5);
           + l5*cos(abs_t5)];
    
    J=[J_t1 J_t2 J_t3 J_t4 J_t5];

end 
% ==========================================================================================

function Xdot= get_Xdot_square(t)
    T=20;
    s=4.5; 
    v=4*s/T;
    
    % right
    if(mod(t,T)<1*T/4)
        xdot = 0;
        ydot = v;
    
    % top
    elseif(mod(t,T)<2*T/4)
        xdot = -v;
        ydot = 0;
    
    % left
    elseif(mod(t,T)<3*T/4)
        xdot = 0;
        ydot = -v;
    
    % bottom
    elseif(mod(t,T)<4*T/4)
        xdot = v;
        ydot = 0;
    
    else
        disp('in else')
        disp(mod(t,T))
    end
    
    Xdot = [xdot;ydot];

end



%% Derivative functions


function Ydot = der_clsdlp_jntspc_ellipse(t,y)
    % y = [t1, t2, t3, t4, t5, td1, td2, td3, td4, td5, phi]
    a=3; b=1;
    % C_e = [1;1];
    r_t=pi/12;
    R = [ cos(r_t), -sin(r_t) ; sin(r_t), cos(r_t)];
    
    l1=2; l2=2.5; l3=1; l4=1.5; l5=1;
    w=2*pi/20;
    phi=y(11);
    
    
    % X_dot
    xdot= a * -sin(phi) * w;
    ydot= b * cos(phi) * w;
    Xdot = R*[xdot;ydot];
    

    J= get_J(y(1:5));
    
    T_e = y(6:10) - y(1:5); % desired - actual
    K = 3*eye(5);   
    
    Jp=pinv(J);
    
    Td_dot = Jp*Xdot;
    T_dot = Jp*Xdot + K*T_e;
    
    Ydot=[T_dot;Td_dot;w];
end


function Ydot = der_clsdlp_tskspc_ellipse(t,y)
    % y = [x, y, xd, yd, t1, t2, t3, t4, t5, phi]
    a=3; b=1;
    % C_e = [1;1];
    r_t=pi/12;
    R = [ cos(r_t), -sin(r_t) ; sin(r_t), cos(r_t)];
    
    l1=2; l2=2.5; l3=1; l4=1.5; l5=1;
    w=2*pi/20;
    phi=y(10);
    
    X_e = y(3:4) - y(1:2);
    K = [5, 0;
         0, 10];
    
    J=get_J(y(5:9));
        
    Jp=pinv(J);

    % dX_des
    xd_dot= a * -sin(phi) * w;
    yd_dot= b * cos(phi) * w;
    Xd_dot = R*[xd_dot;yd_dot];
    
    % dX_act
    X_dot = Xd_dot + K*X_e;
    
    % dTheta_act
    T_dot = Jp*X_dot;
    
    Ydot=[X_dot;Xd_dot;T_dot;w];
    
end

function Y_dot=der_clsdlp_jntspc_square(t,y)
    % y = [t1, t2, t3, t4, t5, td1, td2, td3, td4, td5]
    
    T=20;
    s=4.5; 
    v=4*s/T;
    
    % right
    if(mod(t,T)<1*T/4)
        xdot = 0;
        ydot = v;
    
    % top
    elseif(mod(t,T)<2*T/4)
        xdot = -v;
        ydot = 0;
    
    % left
    elseif(mod(t,T)<3*T/4)
        xdot = 0;
        ydot = -v;
    
    % bottom
    elseif(mod(t,T)<4*T/4)
        xdot = v;
        ydot = 0;
    
    else
        disp('in else')
        disp(mod(t,T))
    end
    
    Xd_dot = [xdot;ydot];
    
    T_e = y(6:10) - y(1:5); 
    K = 3*eye(5);
    
    % actual
    J = get_J(y(1:5));
    Jp=pinv(J);
    
    % desired
    J_d = get_J(y(6:10));
    Jp_d=pinv(J_d);
    
    Td_dot = Jp_d * Xd_dot;
    T_dot = Jp*Xd_dot + K*T_e;
    
    Y_dot=[T_dot;Td_dot];
    
end

function Ydot = der_clsdlp_tskspc_square(t,y)
    % y = [x, y, xd, yd, t1, t2, t3, t4, t5]
    
    % disp(t)

    T=20;
    s=4.5; 
    v=4*s/T;
    
    % right
    if(mod(t,T)<1*T/4)
        xdot = 0;
        ydot = v;

    % top
    elseif(mod(t,T)<2*T/4)
        xdot = -v;
        ydot = 0;
    
    % left
    elseif(mod(t,T)<3*T/4)
        xdot = 0;
        ydot = -v;
    
    % bottom
    elseif(mod(t,T)<4*T/4)
        xdot = v;
        ydot = 0;
    
    else
        disp('in else')
        disp(mod(t,T))
    end
    
    Xd_dot = [xdot;ydot];

    J = get_J(y(5:9));
    Jp=pinv(J);
    
    X_e = y(3:4) - y(1:2);
    K = [5, 0;
         0, 10];
        
    % dX_act
    X_dot = Xd_dot + K*X_e;
    
    % dTheta_act
    T_dot = Jp*X_dot;
    
    Ydot=[X_dot;Xd_dot;T_dot];
end


%% helper functions

function [solution, fval, exitflag, output]= fsolve_5R(pos, guess)
    l1=2; l2=2.5; l3=1; l4=1.5; l5=1;
    
    function F = F_Equations(pos,guess)

        x=pos(1);
        y=pos(2);

        t1=guess(1);
        t2=guess(2);
        t3=guess(3);
        t4=guess(4);
        t5=guess(5);

        % F equations
        f1= x- (l1*cos(t1) + l2*cos(t1+t2) + l3*cos(t1+t2+t3) + l4*cos(t1+t2+t3+t4) + l5*cos(t1+t2+t3+t4+t5));
        f2= y- (l1*sin(t1) + l2*sin(t1+t2) + l3*sin(t1+t2+t3) + l4*sin(t1+t2+t3+t4) + l5*sin(t1+t2+t3+t4+t5));

        F = [f1;f2];
    end

    % guess should be the only variable
    func_F=@(guess) F_Equations(pos,guess);
    [solution, fval, exitflag, output]= fsolve(func_F, guess);

end

function plot_ellipse(Fig,Title,T1,T2,T3,T4,T5)
    
    % link lengths
    l1=2; l2=2.5; l3=1; l4=1.5; l5=1;
    
    % ellipse center
    C_e = [1;1];
    a=3; b=1;
    r_t = pi/12;
    R = [ cos(r_t), -sin(r_t) ; sin(r_t), cos(r_t)];


    x=  l1*cos(T1) + l2*cos(T1+T2) + l3*cos(T1+T2+T3) + l4*cos(T1+T2+T3+T4) + l5*cos(T1+T2+T3+T4+T5);
    y=  l1*sin(T1) + l2*sin(T1+T2) + l3*sin(T1+T2+T3) + l4*sin(T1+T2+T3+T4) + l5*sin(T1+T2+T3+T4+T5);
    
    % EE positions
    X=[x,y];

    clear x y
    for i=[1:360]
        phi=deg2rad(i);
        x(i) = a*cos(phi);
        y(i) = b*sin(phi);
    end

    X_e = C_e + R* [x;y];
    
    % plot
    figure(Fig)
    title(Title)
    for i=1:size(T1,1)
        % disp(i)
        cla;
        hold on;
        xlim([-3,3]); ylim([-1.5,4.5]); 
        axis equal; grid on;
        
        % absolute angles
        t1=T1(i); 
        t2=t1+T2(i);
        t3=t2+T3(i);
        t4=t3+T4(i);
        t5=t4+T5(i);
        
        % joint positions (x5 is EE)
        x1=    l1*cos(t1);     y1=    l1*sin(t1);
        x2= x1+l2*cos(t2);     y2= y1+l2*sin(t2);
        x3= x2+l3*cos(t3);     y3= y2+l3*sin(t3);
        x4= x3+l4*cos(t4);     y4= y3+l4*sin(t4);
        x5= x4+l5*cos(t5);     y5= y4+l5*sin(t5);
    
        % arm
        plot([0, x1],[0, y1],"-b");
        plot([x1,x2],[y1,y2],"-b");
        plot([x2,x3],[y2,y3],"-b");
        plot([x3,x4],[y3,y4],"-b");
        plot([x4,x5],[y4,y5],"-b");
        
        % EE path
        plot(X(1:i,1),X(1:i,2),"-c")
        
        % center of ellipse
        plot(C_e(1),C_e(2),".r");
        
        % desired ellipse
        plot(X_e(1,:),X_e(2,:),"--m")
        
        hold off;
    
        % pause(0.05)
        drawnow;
    end
end


function plot_square(Fig,Title,T1,T2,T3,T4,T5)
    
    % link lengths
    l1=2; l2=2.5; l3=1; l4=1.5; l5=1;
    
    % square
    C_s = [-0.5;1];
    s = 4.5;

    x=  l1*cos(T1) + l2*cos(T1+T2) + l3*cos(T1+T2+T3) + l4*cos(T1+T2+T3+T4) + l5*cos(T1+T2+T3+T4+T5);
    y=  l1*sin(T1) + l2*sin(T1+T2) + l3*sin(T1+T2+T3) + l4*sin(T1+T2+T3+T4) + l5*sin(T1+T2+T3+T4+T5);
    
    % EE positions
    X=[x,y];
    
    % plot
    figure(Fig)
    title(Title)
    for i=1:size(T1,1)
        % disp(i)
        cla;
        hold on;
        xlim([-3,3]); ylim([-2,4]); 
        axis equal; grid on;
        
        % absolute angles
        t1=T1(i); 
        t2=t1+T2(i);
        t3=t2+T3(i);
        t4=t3+T4(i);
        t5=t4+T5(i);
        
        % joint positions (x5 is EE)
        x1=    l1*cos(t1);     y1=    l1*sin(t1);
        x2= x1+l2*cos(t2);     y2= y1+l2*sin(t2);
        x3= x2+l3*cos(t3);     y3= y2+l3*sin(t3);
        x4= x3+l4*cos(t4);     y4= y3+l4*sin(t4);
        x5= x4+l5*cos(t5);     y5= y4+l5*sin(t5);
    
        % arm
        plot([0, x1],[0, y1],"-b");
        plot([x1,x2],[y1,y2],"-b");
        plot([x2,x3],[y2,y3],"-b");
        plot([x3,x4],[y3,y4],"-b");
        plot([x4,x5],[y4,y5],"-b");
        
        % EE path
        plot(X(1:i,1),X(1:i,2),"-c")
        
        % center of ellipse
        plot(C_s(1),C_s(2),".r");
        
        % desired square
        tr= C_s + [s/2;s/2];
        tl= C_s + [-s/2;s/2];
        bl= C_s + [-s/2;-s/2];
        br= C_s + [s/2;-s/2];

        plot([br(1), tr(1)],[br(2), tr(2)],"--m");
        plot([tr(1), tl(1)],[tr(2), tl(2)],"--m");
        plot([tl(1), bl(1)],[tl(2), bl(2)],"--m");
        plot([bl(1), br(1)],[bl(2), br(2)],"--m");
        
    
        hold off;
    
        % pause(0.02)
        drawnow;
    end
end