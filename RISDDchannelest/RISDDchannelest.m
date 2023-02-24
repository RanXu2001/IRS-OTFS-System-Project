%% 代码目标：
%1. 留下一个random phase和一个optimized。
%2. 生成随机信息比特，用OTFS调制。
%3. 信号分别经过random phase 和optimized RIS，在接收端信道估计
%4. 对于每一个SNR，做出channel est
%5. threshold 方法：找到适合的threshold
%6. 只做信道估计，则对于两个信道（random和optimize）和每一个SNR来说，都只需要传输一个frame
%7. 问题：不能估计出LOS和RIS信道分别的参数。
close all;clc;
%% parameters
max_UE_speed_kmh = 100;
max_UE_speed = max_UE_speed_kmh*(1000/3600);

max_RIS_speed_kmh = 0;
max_RIS_speed = max_RIS_speed_kmh*(1000/3600);

% number of Doppler bins (time slots) 
N=16; 
% number of delay bins (subcarriers) 
M=32;

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

%% 生成比特以及调制。
N_fram = 1;%传输1个frame(channel estimation)
M_mod = 4; %4QAM
M_bits = log2(M_mod);%一个symbol有两个bits
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));   %符号的能量
% number of symbols per frame
N_syms_perfram = N*M;   %一帧中包含的符号数=单个载波包含的码元数*载波个数
% number of bits per frame
N_bits_perfram = N*M*M_bits;    %一帧中的符号数*一个符号包含的2bit信息

SNR_dB = 0:5:25;    %信噪比dB从0~15，间隔5dB取样，length为4
SNR = 10.^(SNR_dB/10);  %dB和10进制换算
noise_var_sqrt = sqrt(1./SNR);  %噪声
sigma_2 = abs(eng_sqrt*noise_var_sqrt).^2;

z = exp(1i*2*pi/(M*N));%这个z用于信道估计
for noisei = 1:length(SNR)
data_info_bit = randi([0,1],N_bits_perfram,1); %生成0或1的随机数，生成了一个帧结构中包含的bit数
data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));%二进制数转换成10进制数，reshape为重组数组，语法为reshape(A,a,b)将矩阵A重构为a*b的矩阵
x = qammod(data_temp,M_mod, 'gray'); %输出使用正交幅度调制4-QAM消息信号X的复包络,gray格雷码编码
x = reshape(x,N,M); 
x_ran = x;%经过random信道的frame
x_opt = x;%经过opt信道的frame


%% 生成信道： direct beam
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

%% 生成信道：RIS for each reflecting element:
L = 20;%element number

cascadeChannelPerElement = cell(L,1);

theta = zeros(1,L);
phaseFactor = zeros(1,L);

theta_list = [0,pi*2/8,pi*4/8,pi*6/8,pi*8/8,pi*10/8,pi*12/8,pi*14/8];


for i = 1:L
    theta(i) = theta_list(randi(length(theta_list)));%从theta list里面随机取一个
    phaseFactor(i) = cos(theta(i))+1i*sin(theta(i));
    
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

        %生成h_tot:h_tot是random phase RIS 的信道响应
        h = phaseFactor(i)*h;
       
        if i==1
            h_tot = h;
        else
            h_tot = h_tot+h;
        end

        %P:总的path数
        

end

%% 关于RIS 优化：
%得到全部P条path的delay和Doppler shift，以及对应的直射路的增益和过每个RIS element的级联信道的增益
theta_optimized = zeros(1,L);

a_extended = padarray(a,[N-1 M-1],0,'post');

P = length(find(h_tot~=0));
for ith = 1:P
    %[l_tot,k_tot]= find(h_tot);%l_tot, k_tot是RIS+LOS的所有时延和Doppler
    [k_tot,l_tot]= find(h_tot);
end
%问题：要在生成所有的信道冲激响应矩阵之后将他们保存。（保存每一个element的级联矩阵）

%选出the path that maximizes the upper bound on the effective channel gain:
target = zeros(P,1);
for ith = 1:P %i 是径
    %target(ith) =abs(a_extended(l_tot(ith),k_tot(ith))) ;
    target(ith) =abs(a_extended(k_tot(ith),l_tot(ith))) ;
    for jth = 1:L%j是element
        h_temp = cascadeChannelPerElement{jth,1};
        %target(ith) = target(ith) +abs(h_temp(l_tot(ith),k_tot(ith))) ;
        target(ith) = target(ith) +abs(h_temp(k_tot(ith),l_tot(ith))) ;
    end
    target(ith) = target(ith)^2;
end
[targetVal,p] = max(target);
%得到优化后的调制角度和最优化后的信道响应矩阵
h_tot_optimized = a_extended;
for ith = 1:L
    h_temp = cascadeChannelPerElement{ith,1};
    %theta_optimized(ith) = angle(a_extended(l_tot(p),k_tot(p))) - angle(h_temp(l_tot(p),k_tot(p)));
    theta_optimized(ith) = angle(a_extended(k_tot(p),l_tot(p))) - angle(h_temp(k_tot(p),l_tot(p)));
    h_tot_optimized = h_tot_optimized + (cos(theta_optimized(ith))+1i*sin(theta_optimized(ith)))*h_temp;
end
h_tot_optimized_abs = abs(h_tot_optimized);


%% 

%加上直射链路
h_tot = h_tot + a_extended;
h_abs = abs(h);
h_tot_abs = abs(h_tot);

%% 从以上，得到h_tot和h_tot_optimized.
[ki_ran,li_ran] = find(h_tot);
kmax_ran = max(ki_ran);
lmax_ran = max(li_ran);
mp_ran = lmax_ran+1;
np_ran = floor(N/2);
x_ran((np_ran-2*kmax_ran):(np_ran+2*kmax_ran),(mp_ran-lmax_ran):(mp_ran+lmax_ran))=0;
x_ran(np_ran,mp_ran) = 1;
x_ran_abs = abs(x_ran);

[ki_opt,li_opt] = find(h_tot_optimized);
kmax_opt = max(ki_opt);
lmax_opt = max(li_opt);
mp_opt = lmax_opt+1;
np_opt = floor(N/2);
x_opt((np_opt-2*kmax_opt):(np_opt+2*kmax_opt),(mp_opt-lmax_opt):(mp_opt+lmax_opt))=0;
x_opt(np_opt,mp_opt) = 1;
x_opt_abs = abs(x_opt);

%% OTFS modulation
s_ran = OTFS_modulation(N,M,x_ran);
s_opt = OTFS_modulation(N,M,x_opt);

%% OTFS channel output
taps_ran = length(li_ran);
chan_coef_ran = zeros(1,taps_ran);
for i = 1:taps_ran
    chan_coef_ran(i) = h_tot(ki_ran(i),li_ran(i));
end
r_ran = OTFS_channel_output(N,M,taps_ran,li_ran,ki_ran,chan_coef_ran,sigma_2(noisei),s_ran);

taps_opt = length(li_opt);
chan_coef_opt = zeros(1,taps_opt);
for i = 1:taps_opt
    chan_coef_opt(i) = h_tot_optimized(ki_opt(i),li_opt(i));
end
r_opt = OTFS_channel_output(N,M,taps_opt,li_opt,ki_opt,chan_coef_opt,sigma_2(noisei),s_opt);

%%  OTFS demodulation
y_ran = OTFS_demodulation(N,M,r_ran);%y:N*M
y_ran_abs = abs(y_ran);

y_opt = OTFS_demodulation(N,M,r_opt);%y:N*M
y_opt_abs = abs(y_opt);

 %% RIS-OTFS channel estimation
 h_est_ran = zeros((2*N-1),(2*M-1));
 for l = 1:(lmax_ran+1)
      for k = 1:(kmax_ran+1)
          h_est_ran(k,l) = y_ran((np_ran+k),(mp_ran+l))/(1*z^(k*mp_ran));
      end
 end
 h_est_ran_abs = abs(h_est_ran);

 h_est_opt = zeros((2*N-1),(2*M-1));
 for l = 1:(lmax_opt+1)
      for k = 1:(kmax_opt+1)
          h_est_opt(k,l) = y_opt((np_opt+k),(mp_opt+l))/(1*z^(k*mp_opt));
      end
 end
 h_est_opt_abs = abs(h_est_opt);

%% 画图
%经过RIS反射之后的信道冲激响应
% figure;
% bar3(h_tot_abs);
% xlabel('delay');
% ylabel('Doppler');
% zlabel('Channel coeficient(abs)');
% title('Random RIS phase channel impulse response 1');
% 
% 
% figure;
% bar3(h_abs);
% xlabel('delay');
% ylabel('Doppler');
% zlabel('Channel coeficient(abs)');
% 
% figure;
% bar3(abs(a_extended));
% xlabel('delay');
% ylabel('Doppler');
% zlabel('Channel coeficient(abs)');
% title('Direct path');
% 
% figure;
% bar3(abs(h_tot_optimized));
% xlabel('delay');
% ylabel('Doppler');
% zlabel('Channel coeficient(abs)');
% title('Optimized channel impulse response');
figure;
bar3(abs(h_est_opt));
title('Estimated optimized channel impulse response');
figure;
bar3(abs(h_tot_optimized));
title('Optimized channel impulse response');

% figure;
% bar3(abs(h_est_ran));
% title('Estimated random phase channel impulse response');
% figure;
% bar3(abs(h_tot));
% title('Random phase channel impulse response');
end