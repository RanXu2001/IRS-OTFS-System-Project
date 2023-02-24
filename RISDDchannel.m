close all;clc;
%% parameters
max_UE_speed_kmh = 1000;
max_UE_speed = max_UE_speed_kmh*(1000/3600);

max_RIS_speed_kmh = 0;
max_RIS_speed = max_RIS_speed_kmh*(1000/3600);

% number of Doppler bins (time slots) 
N=16; 
% number of delay bins (subcarriers) 
M=64;

% subcarrier spacing 
delta_f=15e3;

% block duration
T=1/delta_f;

% carrier frequency 
fc=4e9;

delays_EVA = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510]*10^(-9); 
pdp_EVA = [0.0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];

delays_EPA = [0, 30, 70, 90, 110, 190, 410]*10^(-9); 
pdp_EPA = [0.0, -1.0, -2.0, -3.0, -8.0, -17.2, -20.8];

delays_ETU = [0, 50, 120, 200, 230, 500, 1600, 2300, 5000]*10^(-9);
pdp_ETU = [-1.0 -1.0 -1.0 0.0 0.0 0.0 -3.0 -5.0 -7.0];

%speed of light
c=299792458;

% OTFS grid delay and Doppler resolution 
delay_resolution = 1/(M*delta_f);
Doppler_resolution = 1/(N*T);


%% direct beam
%maximum Doppler spread(one side)
nu_d_max = (max_UE_speed*fc)/(c);
        
% maximum normalized Doppler spread (one-sided) 
k_d_max = nu_d_max/Doppler_resolution;

delays_d = delays_EVA;
pdp_d=pdp_EVA;
pdp_d_linear = 10.^(pdp_d/10);       
pdp_d_linear = pdp_d_linear/sum(pdp_d_linear);
taps_d=length(pdp_d);
a_i = sqrt(pdp_d_linear).*(sqrt(1/2) * (randn(1,taps_d)+1i*randn(1,taps_d)));
%delay tapes
l_d_i=  round(delays_d./delay_resolution);  %integer delay
k_d_i = round(k_d_max-k_d_max*rand(1,taps_d));%integer Doppler
%DD域响应转换成矩阵
a = zeros(N,M);
  for ith = 1:length(a_i)
          delay = l_d_i(ith)+1;
          Doppler = k_d_i(ith)+1;
          a(Doppler,delay) = a_i(ith);
  end
        a_abs = abs(a);

%% for each reflecting element:
L = 20;%element number

cascadeChannelPerElement = cell(L,1);

theta = zeros(1,L);
phaseFactor = zeros(1,L);

% 多组theta值
theta1 = zeros(1,L);
phaseFactor1 = zeros(1,L);
theta2 = zeros(1,L);
phaseFactor2 = zeros(1,L);
theta3 = zeros(1,L);
phaseFactor3 = zeros(1,L);
theta4 = zeros(1,L);
phaseFactor4 = zeros(1,L);
theta_list = [0,pi*2/8,pi*4/8,pi*6/8,pi*8/8,pi*10/8,pi*12/8,pi*14/8];





for i = 1:L
    theta(i) = theta_list(randi(length(theta_list)));%从theta list里面随机取一个
    phaseFactor(i) = cos(theta(i))+1i*sin(theta(i));
    % 多组theta值
    theta1(i) = theta_list(randi(length(theta_list)));
    phaseFactor1(i) = cos(theta1(i))+1i*sin(theta1(i));
    theta2(i) = theta_list(randi(length(theta_list)));
    phaseFactor2(i) = cos(theta2(i))+1i*sin(theta2(i));
    theta3(i) = theta_list(randi(length(theta_list)));
    phaseFactor3(i) = cos(theta3(i))+1i*sin(theta3(i));
    theta4(i) = theta_list(randi(length(theta_list)));
    phaseFactor4(i) = cos(theta4(i))+1i*sin(theta4(i));
    
    %% RIS-User link
        %maximum Doppler spread(one side)
        nu_max = (max_UE_speed*fc)/(c);
        
        % maximum normalized Doppler spread (one-sided) 
        k_max = nu_max/Doppler_resolution;
        
        delays=delays_EVA; 
        pdp=pdp_EVA;
        pdp_linear = 10.^(pdp/10);
        pdp_linear = pdp_linear/sum(pdp_linear);
        taps=length(pdp);
        g_i = sqrt(pdp_linear).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));
        %delay tapes
        l_i=  round(delays./delay_resolution);  %integer delay
        %l_i= delays./delay_resolution;
        %doppler tapes
        %k_i = round(k_max*cos(2*pi*rand(1,taps)));%integer Doppler Jake spectrum
        k_i = round(k_max-k_max*rand(1,taps));%integer Doppler uniform spectrum, all positive.
        %k_i = (k_max*cos(2*pi*rand(1,taps)));
        
        
        %% BS-RIS link
        %maximum Doppler spread(one side)
        nu1_max = (max_RIS_speed*fc)/(c);
        % maximum normalized Doppler spread (one-sided) 
        k1_max = nu1_max/Doppler_resolution;
        
        delays1 = delays_EVA;
        pdp1 = pdp_EVA;
        pdp1_linear = 10.^(pdp1/10);
        pdp1_linear = pdp1_linear/sum(pdp1_linear);
        taps1=length(pdp1);
        %gain
        f_i = sqrt(pdp1_linear).*(sqrt(1/2) * (randn(1,taps1)+1i*randn(1,taps1)));
        %delay tapes
        l1_i=  round(delays1./delay_resolution);  %integer delay
        %doppler tapes 问题：怎么确定两个link的Doppler shift的max。
        k1_i = round(k1_max*cos(2*pi*rand(1,taps1)));%integer Doppler Jake spectrum
        
        %% 得到DD域信道冲激响应（矩阵形式，横坐标为Doppler，纵坐标为delay）
        %RIS-User
        %g_i,l_i,k_i
        g = zeros(N,M);
        
        for ith = 1:length(g_i)
            delay = l_i(ith)+1;
            Doppler = k_i(ith)+1;
            g(Doppler,delay) = g_i(ith);
        end
        g_abs = abs(g);
        
        %BS-RIS
        f = zeros(N,M);
        for ith = 1:length(f_i)
            delay = l1_i(ith)+1;
            Doppler = k1_i(ith)+1;
            f(Doppler,delay) = f_i(ith);
        end
        f_abs = abs(f);

        h = conv2(f,g); %h是RIS单反射面级联信道冲激

        %将单反射面级联信道储存
        cascadeChannelPerElement(i) = {h};

        %生成h_tot
        h = phaseFactor(i)*h;
       
        % 多组theta值
        h1 = phaseFactor1(i)*h;
        h2 = phaseFactor2(i)*h;
        h3 = phaseFactor3(i)*h;
        h4 = phaseFactor4(i)*h;

        if i==1
            h_tot = h;
        else
            h_tot = h_tot+h;
        end

        if i==1
            h_tot1 = h1;
        else
            h_tot1 = h_tot1+h1;
        end
        if i==1
            h_tot2 = h2;
        else
            h_tot2 = h_tot2+h2;
        end
        if i==1
            h_tot3 = h3;
        else
            h_tot3 = h_tot3+h3;
        end
        if i==1
            h_tot4 = h4;
        else
            h_tot4 = h_tot4+h4;
        end
        %P:总的path数
        

end

%% 关于RIS 优化：
%得到全部P条path的delay和Doppler shift，以及对应的直射路的增益和过每个RIS element的级联信道的增益
theta_optimized = zeros(1,L);

a_extended = padarray(a,[N-1 M-1],0,'post');

P = length(find(h_tot~=0));
for ith = 1:P
    [l_tot,k_tot]= find(h_tot);%l_tot, k_tot是RIS+LOS的所有时延和Doppler
end
%问题：要在生成所有的信道冲激响应矩阵之后将他们保存。（保存每一个element的级联矩阵）

%选出the path that maximizes the upper bound on the effective channel gain:
target = zeros(P,1);
for ith = 1:P %i 是径
    target(ith) =abs(a_extended(l_tot(ith),k_tot(ith))) ;
    for jth = 1:L%j是element
        h_temp = cascadeChannelPerElement{jth,1};
        target(ith) = target(ith) +abs(h_temp(l_tot(ith),k_tot(ith))) ;
    end
    target(ith) = target(ith)^2;
end
[targetVal,p] = max(target);
%得到优化后的调制角度和最优化后的信道响应矩阵
h_tot_optimized = a_extended;
for ith = 1:L
    h_temp = cascadeChannelPerElement{ith,1};
    theta_optimized(ith) = angle(a_extended(l_tot(p),k_tot(p))) - angle(h_temp(l_tot(p),k_tot(p)));
    h_tot_optimized = h_tot_optimized + (cos(theta_optimized(ith))+1i*sin(theta_optimized(ith)))*h_temp;
end


%% 

%加上直射链路
h_tot = h_tot + a_extended;
h_tot1 = h_tot1 + a_extended;
h_tot2 = h_tot2 + a_extended;
h_tot3 = h_tot3 + a_extended;
h_tot4 = h_tot4 + a_extended;

h_abs = abs(h);

h_tot_abs = abs(h_tot);
h_tot1_abs = abs(h_tot1);
h_tot2_abs = abs(h_tot2);
h_tot3_abs = abs(h_tot3);
h_tot4_abs = abs(h_tot4);

%经过RIS反射之后的信道冲激响应
figure;
bar3(h_tot_abs);
xlabel('delay');
ylabel('Doppler');
zlabel('Channel coeficient(abs)');
title('Random RIS phase channel impulse response 1');
figure;
bar3(h_tot1_abs);
xlabel('delay');
ylabel('Doppler');
zlabel('Channel coeficient(abs)');
title('Random RIS phase channel impulse response 2');
figure;
bar3(h_tot2_abs);
xlabel('delay');
ylabel('Doppler');
zlabel('Channel coeficient(abs)');
title('Random RIS phase channel impulse response 3');
figure;
bar3(h_tot3_abs);
xlabel('delay');
ylabel('Doppler');
zlabel('Channel coeficient(abs)');
title('Random RIS phase channel impulse response 4');
figure;
bar3(h_tot4_abs);
xlabel('delay');
ylabel('Doppler');
zlabel('Channel coeficient(abs)');
title('Random RIS phase channel impulse response 5');

figure;
bar3(h_abs);
xlabel('delay');
ylabel('Doppler');
zlabel('Channel coeficient(abs)');

figure;
bar3(abs(a_extended));
xlabel('delay');
ylabel('Doppler');
zlabel('Channel coeficient(abs)');
title('Direct path');

figure;
bar3(abs(h_tot_optimized));
xlabel('delay');
ylabel('Doppler');
zlabel('Channel coeficient(abs)');
title('Optimized channel impulse response');

