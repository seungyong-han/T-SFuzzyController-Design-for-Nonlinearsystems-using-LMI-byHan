%% 240813 Nonlinear Dynamics Model, 2 rules
% by Seungyong Han
clear
clc
close all
%%
global a b
a = 1.4;
b = 0.7;
n_r = 2;
pre_x1 = [a];
pre_x3 = [0 b];

%% 
K{1}= [-0.46175703 -0.69774825 -0.61291852  0.04828783];
K{2}= [-0.4629239  -0.69790545 -0.58677045  0.04121685];

%%
tf=20;  % final time
ti=0.01;  % runge kutta sample time 
tspan=0:ti:tf;
sample_size = size(tspan,2);

% Theoretical range
% x(1) \in [-a a] = [-1.4 1.4]
% x(3) \in [-b b] = [-0.7, 0.7]

x(:,1) = [1 0.5 0.7 -0.6]
% simulation success upto x(1,1) \in [-20 18]

U=0;
u_temp(:,1)=U;

    for i=1:sample_size-1
    W(:,i) = 0*exp(-0.1*i*ti)*sin(1*i*ti);
        
%     if x(1,i) > -a && x(1,i) < a 
    M{1} = (x(1,i)^2)/(a^2); % x1=0인 지점에서 선형화하였으나 모든 x1의 범위에 대해서 sub system을 적용한다는 의미
%     else  
        
    if x(3,i) ==0
        N{1} = 1;
    else
        N{1} = (b*sin(x(3,i)-x(3,i)*sin(b)))/(x(3,i)*(b-sin(b)));
    end
        N{2} = 1 - N{1};

 
    %% normalized membership function
    for k1 = 1:size(pre_x1,2)
        for k2 = 1:size(pre_x3,2)
            w{k1,k2} = M{k1}*N{k2};
        end
    end
    
    sum_w = 0;
    for kk = 1:n_r
        sum_w = sum_w + sum(w{kk});
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
    
    x(:,i+1)=rk5_4(x(:,i),U+W(:,i),ti); 
    u_temp(:,i+1)=U;

    end
figure()
plot(tspan,x(1,:),'r');
hold on
plot(tspan,x(2,:),'b');
plot(tspan,x(3,:),'m');
plot(tspan,x(4,:),'g');
xlabel('Time')
legend('x_{1}','x_{2}','x_{3}','x_{4}')
grid on

% figure()
% plot(tspan,u_temp(1,:),'r')
% grid on
% legend('u_{1}')
% xlabel('Time');
% axis([0 10 -10 10])


