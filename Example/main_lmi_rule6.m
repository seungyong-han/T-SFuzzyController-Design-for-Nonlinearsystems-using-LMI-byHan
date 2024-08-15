%% 240815 Inverted pendulum for rule 6
% By Seungyong Han
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
n_r = 6;
a = pi/3;  
b = 1;        
pre_x1 = [-a 0 a];
pre_x3 = [-b, b];

K{1} = [160.00237196  41.94243039   1.74557348  57.68074125];
K{2} = [161.17644047  42.72651613   1.77527554  58.37913544];
K{3} = [607.6785384  182.70637164   6.40284217 159.81816265];
K{4} = [613.97890776 178.49036058   6.49907699 160.62573939];
K{5} = [169.80604397  47.83264729   1.98638712  63.67989567];
K{6} = [169.40514833  47.6316946    1.98509784  63.6627364 ];

%%
tf=60;  % final time
ti=0.001;  % runge kutta sample time  fail :0.00001
tspan=0:ti:tf;
sample_size = size(tspan,2);
deg =   -71; % -71 ~ 71
x(:,1) = [deg*(pi/180);0;0;0];

U=0*x(:,1);
u_temp(:,1)=U;
    for i=1:sample_size-1
    W(:,i) = 0*exp(-0.1*i*ti)*sin(i*ti);
    
    z1 = x(1,i);
    z2 = x(2,i);
    
    if z1 < -a
        Mf{1} = 1;
    elseif z1 > 0
        Mf{1} = 0;
    else
        Mf{1} = (z1/-a);
    end
    
    if z1 < -a || z1 > a
        Mf{2} = 0;
    else
        Mf{2} = ((z1+a)/a)*((z1-a)/-a);
    end
    
    if z1 < 0
        Mf{3} = 0;
    elseif z1 > a
        Mf{3} = 1;
    else
        Mf{3} = z1/a;
    end
    
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


