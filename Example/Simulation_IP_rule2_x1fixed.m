%% 240815 4 by 4 Inverted pendulum rule 2
% by Seungyong Han
% x1 == 10*pi/180, x2 \in [-b, b]
clear
clc
close all
%%
global f0 f1 m M l J g
m = 0.22;
M = 1.3282;
f0 = 22.915;
f1 = 0.007056;
l = 0.304;
J = 0.004963;
g = 9.8;
n_r = 2;
K{1}= [33.90392585  3.36177302  1.14239616 25.31916985];
K{2}= [33.90392736  3.33854645  1.14239633 25.31917022];
%%
% load('h_inf_TS_gain_4by4.mat')

tf=20;  % final time
ti=0.001;  % runge kutta sample time  fail :0.00001
tspan=0:ti:tf;
sample_size = size(tspan,2);
deg =   -52; % -52 ~ 52
x(:,1) = [deg*(pi/180);0;0;0];
a = 10*pi/180;
b = 1;
pre_x1 = [a];
pre_x3 = [-b, b];
U=0*x(:,1);
u_temp(:,1)=U;
for i=1:sample_size-1
    W(:,i) = 0*exp(-0.1*i*ti)*sin(i*ti);
     
    z1 = x(1,i);
    z2 = x(2,i);
    
 
    Mf{1} = 1;
    if z2 < -b
        N{1} = 1;
    elseif z2 > b
        N{1} = 0;
    else
        N{1} = ((z2-b)/(-2*b));
    end
    
    if z2 < -b
        N{2} = 0;
    elseif z2 > b
        N{2} = 1;
    else
        N{2} = ((z2+b)/(2*b));
    end
    
    
%% normalized membership function
    for k1 = 1:size(pre_x1,2)
        for k2 = 1:size(pre_x3,2)
            w{k1,k2} = Mf{k1}*N{k2};
        end
    end
    
    sum_w = 0;
    for kk = 1:n_r
        sum_w = sum_w + sum(w{kk}); % if designer did not consider min, max of premise variable, then summation rule is importnat
    end
    
    num_k = 1;
    for k1 = 1:size(pre_x1,2)
        for k2 = 1:size(pre_x3,2)
            h{num_k} = w{k1,k2}/sum_w;
            num_k  = num_k + 1;
        end
    end
       
    gain = 0;
    for kk = 1:n_r
        gain = gain + h{kk}*K{kk};
    end
    
    U= gain* x(:,i);

    
    x(:,i+1)=rk5_2(x(:,i),U+W(:,i),ti); 
    u_temp(:,i+1)=U;


    end
figure(1)
plot(tspan,x(1,:),'r');
hold on
plot(tspan,x(2,:),'b');
plot(tspan,x(3,:),'m');
plot(tspan,x(4,:),'g');
xlabel('Time')
legend('x_{1}','x_{2}','x_{3}','x_{4}')
grid on

figure(2)
plot(tspan,u_temp(1,:),'r')
grid on
legend('u_{1}')
xlabel('Time');
% axis([0 10 -10 10])


