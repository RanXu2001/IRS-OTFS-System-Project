%% DD域
max_UE_speed_kmh = 1000;
max_UE_speed = max_UE_speed_kmh*(1000/3600);

N=16; 
M=64;

delta_f=15e3;
T=1/delta_f;
fc=4e9;

delays_EVA = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510]*10^(-9); 
pdp_EVA = [0.0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];

c=299792458;

delay_resolution = 1/(M*delta_f);
Doppler_resolution = 1/(N*T);

nu_max = (max_UE_speed*fc)/(c);
k_max = nu_max/Doppler_resolution;

delays=delays_EVA; 
pdp=pdp_EVA;       
pdp_linear = 10.^(pdp/10);
pdp_linear = pdp_linear/sum(pdp_linear);
taps=length(pdp);
g_i = sqrt(pdp_linear).*(sqrt(1/2) * (randn(1,taps)+1i*randn(1,taps)));

l_i=  round(delays./delay_resolution);
k_i = round(k_max-k_max*rand(1,taps));

%% 得到DD域信道冲激响应（矩阵形式，横坐标为Doppler，纵坐标为delay）
g = zeros(M,N);

for ith = 1:length(g_i)
    delay = l_i(ith)+1;
    Doppler = k_i(ith)+1;
    g(Doppler,delay) = g_i(ith);
end
g_abs = abs(g);
figure;
bar3(g_abs);
title('Channel response in delay-Doppler domain');
xlabel('delay');
ylabel('Doppler shift');
zlabel('Channel coeficient(abs)');

%% 得到TD域冲激响应
g_td = zeros(M,N);
for ith = 1:M
    g_td(ith,:) = fft(g(ith,:));
end
figure;
g_td_abs = abs(g_td);
surf(g_td_abs);
title('Channel response in delay-time domain');
xlabel('delay');
ylabel('time');
zlabel('Channel coeficient(abs)');
%% 得到TF域冲激响应
g_tf = zeros(M,N);
for ith = 1:M
    g_tf(ith,:) = fft(g(ith,:));
end
for ith = 1:N
    g_tf(:,ith) = fft(g_tf(:,ith));
end
figure;
g_tf_abs = abs(g_tf);
surf(g_tf_abs);
title('Channel response in time-frequency domain');
xlabel('frequency');
ylabel('time');
zlabel('Channel coeficient(abs)');
