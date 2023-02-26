%
% Copyright (c) 2018, Raviteja Patchava, Yi Hong, and Emanuele Viterbo, Monash University
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%    - Latest version of this code may be downloaded from: https://ecse.monash.edu/staff/eviterbo/
%    - Freely distributed for educational and research purposes
%%

clc
clear all
close all
tic
%% OTFS parameters%%%%%%%%%%
% number of symbol
N = 32; %单个载波包含的码元数为8
% number of subcarriers
M = 16;  %载波数为8
% size of constellation
M_mod = 4;  %进制数为4，那么一码元就包含2比特
M_bits = log2(M_mod);
% average energy per data symbol
eng_sqrt = (M_mod==2)+(M_mod~=2)*sqrt((M_mod-1)/6*(2^2));   %符号的能量
% number of symbols per frame
N_syms_perfram = N*M;   %一帧中包含的符号数=单个载波包含的码元数*载波个数
% number of bits per frame
N_bits_perfram = N*M*M_bits;    %一帧中的符号数*一个符号包含的2bit信息

z = exp(1i*2*pi/(M*N));

SNR_dB = 0:5:15;    %信噪比dB从0~15，间隔5dB取样，length为4
SNR = 10.^(SNR_dB/10);  %dB和10进制换算
noise_var_sqrt = sqrt(1./SNR);  %噪声
sigma_2 = abs(eng_sqrt*noise_var_sqrt).^2;
%%
rng(1)
N_fram = 500;
err_ber = zeros(length(SNR_dB),1); %生成0矩阵，length返回 X 中最大数组维度的长度，此时为4,即此命令生成4*1维的0矩阵

fprintf('Start to do OTFS simulation. Please wait...\n');
for iesn0 = 1:length(SNR_dB)    %1~4
    SNR_temp = SNR_dB(iesn0);   %SNR_temp取0,5,10,15
    SNR_temp 
    for ifram = 1:N_fram 
        %% random input bits generation and 4QAM modulation%%%%%
        data_info_bit = randi([0,1],N_bits_perfram,1); %生成0或1的随机数，且维数是128*1，也就是生成了一个帧结构中包含的bit数
        data_temp = bi2de(reshape(data_info_bit,N_syms_perfram,M_bits));%二进制数转换成10进制数，reshape为重组数组，语法为reshape(A,a,b)将矩阵A重构为a*b的矩阵
        %reshape(128,64,2)，即将128*1的data_info_bit重组成64*2的矩阵，
        %x = qammod(data_temp,M_mod,0,'gray');
        x = qammod(data_temp,M_mod, 'gray'); %输出使用正交幅度调制4-QAM消息信号X的复包络,gray格雷码编码
        x = reshape(x,N,M); %将64*1的x重组成8*8数组


        %% OTFS channel generation%%%%
        [taps,delay_taps,Doppler_taps,chan_coef] = OTFS_channel_gen(N,M);
        g = zeros(M,N);% channel response
        for ith = 1:length(chan_coef)
            delay = delay_taps(ith)+1;
            Doppler = Doppler_taps(ith)+1;
            g(Doppler,delay) = chan_coef(ith);
        end
        g_abs = abs(g);

        %% 加导频：在第m_p行, n_p列添加导频，周围一圈是0，作为guard。
        lmax = max(delay_taps);
        kmax = max(Doppler_taps);
        m_p = lmax+1;
        n_p =floor(N/2);
        x((m_p-lmax):(m_p+lmax),(n_p-2*kmax):(n_p+2*kmax))=0;
        x(m_p,n_p) = 1;
        x_abs = abs(x);

        %% OTFS modulation%%%%
        s = OTFS_modulation(N,M,x);
        
        
        %% OTFS channel output%%%%%
        r = OTFS_channel_output(N,M,taps,delay_taps,Doppler_taps,chan_coef,sigma_2(iesn0),s);
        
        %% OTFS demodulation%%%%
        y = OTFS_demodulation(N,M,r);%y:M*N
        y_abs = abs(y);
        %% OTFS channel estimation
        g_est = zeros(M,N);
        for l = 1:(lmax+1)
            for k = 1:(kmax+1)
                g_est(l,k) = y((m_p+l-1),(n_p+k-1))/(1*z^(k*m_p));
            end
        end
        g_est_abs = abs(g_est);
        
        %% message passing detector%%%%
        x_est = OTFS_mp_detector(N,M,M_mod,taps,delay_taps,Doppler_taps,chan_coef,sigma_2(iesn0),y);
        
        %% output bits and errors count%%%%%
       
        data_demapping = qamdemod(x_est,M_mod,'gray');
        data_info_est = reshape(de2bi(data_demapping,M_bits),N_bits_perfram,1);
        errors = sum(xor(data_info_est,data_info_bit));
        err_ber(iesn0) = errors + err_ber(iesn0);        
    end
    err_ber_fram_temp = err_ber(iesn0) / N_bits_perfram / N_fram;
    err_ber_fram_temp      %表达式加;可以隐藏输出
    err_ber_fram(iesn0) = err_ber_fram_temp;
end
%err_ber_fram = err_ber/N_bits_perfram./N_fram;这个一直注释
err_ber_fram;
semilogy(SNR_dB, err_ber_fram,'-*','LineWidth',2);
title(sprintf('OTFS'))
ylabel('BER'); xlabel('SNR in dB');grid on

figure;
        bar3(g_est_abs);
        title('g_est_abs');
        figure;
        bar3(g_abs);
        title( 'g_abs');
