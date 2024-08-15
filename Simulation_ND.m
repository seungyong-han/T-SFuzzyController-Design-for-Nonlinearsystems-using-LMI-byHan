%% 240815 Nonlinear Dynamics Model
% by Seungyong Han
% NOTE
% 1. M(fuzzy set) should be constructed in the premise variable dependent form, not system state.
% 2. If designer did not consider min, max of premise variable, then summation rule is importnat.

clear
clc
close all
%%
global a b
a = 1.4;
b = 0.7;
n_r = 12;
pre_x1 = [0, a/2, a];
pre_x3 = [-b -b/2, b/2, b];

%% 
K{1}= [-1.32917195 -2.76154173 -2.56703181  0.42378244];
K{2}= [-1.34340337 -2.73687678 -2.60494711  0.40034077];
K{3}= [-1.34340337 -2.73687678 -2.60494711  0.40034077];
K{4}= [-1.31834486 -3.75673021 -3.21218563  0.72266578];
K{5}= [-0.89800759 -2.43220558 -2.15645429  0.49362814];
K{6}= [-0.89611137 -2.4212017  -2.19928298  0.51755707];
K{7}= [-0.89611137 -2.4212017  -2.19928298  0.51755707];
K{8}= [-0.89800759 -2.43220558 -2.15645429  0.49362814];
K{9}= [-0.42260042 -1.0382448  -1.08037547  0.24956538];
K{10}= [-0.42299112 -1.04067044 -1.10125     0.24208813];
K{11}= [-0.42299112 -1.04067044 -1.10125     0.24208813];
K{12}= [-0.42260042 -1.0382448  -1.08037547  0.24956538]; 
%%
tf=20;  % final time
ti=0.01;  % runge kutta sample time 
tspan=0:ti:tf;
sample_size = size(tspan,2);

% Theoretical range
% x(1) \in [-a a] = [-1.4 1.4]
% x(3) \in [-b b] = [-0.7, 0.7]

x(:,1) = [-30 0.5 0.7 -0.6];  
% simulation success upto x(1,1) \in [-30 28]

U=0;
u_temp(:,1)=U;

    for i=1:sample_size-1
    W(:,i) = 0*exp(-0.1*i*ti)*sin(1*i*ti);
    
    z1 = x(1,i)^2; % premise variable!, M should be constructed in the premise variable dependent form, not system state
    z2 = x(3,i);
    
        if z1 < 0
        M{1} = 1;
    elseif z1 > (a/2)^2
        M{1} = 0;
    else
        M{1} = (z1 - (a/2)^2)/(0 - (a/2)^2);
    end
    
    if z1 < 0
        M{2} = 0;
    elseif z1 >  a^2
        M{2} = 0;
    else
        M{2} = ((z1 - 0)/((a/2)^2 - 0))*((z1 - a^2)/((a/2)^2 - a^2));
    end
    
    if z1 < (a/2)^2
        M{3} = 0;
    elseif z1 > a^2
        M{3} = 1;
    else
        M{3} = (z1 - (a/2)^2)/(a^2 - (a/2)^2);
    end
    
   
    if z2 < -b
        N{1} = 1;
    elseif z2 > -b/2
        N{1} = 0;
    else
        N{1} = (z2+b/2)/(-b+b/2);
    end
    
    if z2 < -b || z2 > b/2
        N{2} = 0;
    else
        N{2} = ((z2+b)/(-b/2 + b))*((z2-b/2)/(-b/2 -b/2));
    end
    
    if z2 < -b/2 || z2 > b
        N{3} = 0;
    else
        N{3} = ((z2+b/2)/(b/2 + b/2))*((z2-b)/(b/2 -b));
    end
 
   if z2 < b/2
        N{4} = 0;
   elseif z2 > b
        N{4} = 1;
    else
        N{4} = ((z2-b/2)/(b - b/2));
    end
 
    %% normalized membership function
    for k1 = 1:size(pre_x1,2)
        for k2 = 1:size(pre_x3,2)
            w{k1,k2} = M{k1}*N{k2};
        end
    end
    
    sum_w = 0;
    for kk = 1:n_r
        sum_w = sum_w + sum(w{kk}); % if designer did not consider min, max of premise variable, then summation rule is importnat
    end
    
    
%     if isnan(sum_w)==1
%         disp('error!!!!!!!!!!!')
%     end
    num_k = 1;
    sum_h = 0;
    for k1 = 1:size(pre_x1,2)
        for k2 = 1:size(pre_x3,2)
            h{num_k} = w{k1,k2}/sum_w;
            sum_h = sum_h + h{num_k};
            num_k  = num_k + 1;
            
        end
    end
%     sum_h
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


