# IRS-OTFS-System-Project
## Doubly-selective channel in different domains
###channel_in_different_domains.m
## RIS-OTFS channel generation and RIS optimization
###RISDDchannel.m
## Channel estimation(The effective channel response)
###RISDDchannelest.m
## Channel estimation(The RIS-reflected channel and LOS channel)
###RISDDchannelest_separated.m
## RIS优化前后的信道冲激响应对比图，优化前后的概率密度函数(用曲线图体现)
## DD域信道估计在不同SNR下，估计的效果。将原信道矩阵和估计信道矩阵求相关系数
## OTFS系统的系统容量曲线，RIS-OTFS系统的系统容量曲线
## 关于OTFS仿真：
1.	体现OTFS在高速信道下的表现很好（根据速度不同跑不同的曲线）
2.	和OFDM的对比，不同的速度下，用OTFS和OFDM跑下来的曲线
3.	OTFS的信道估计方法：数学原理和MATLAB实现
4.	不同检测方法下OTFS的性能（可有可无）
5.	Zak变换和普通方法的差别
