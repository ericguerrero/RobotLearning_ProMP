clc;close all;clear all; format compact;

%% Demostrations
num = 100; % # of demostrations
position = [1, 2.5, 3]; 
time = [0 30 60];  
[trajT, trajX] = trajGeneration(num,position,time);

%% 
n = 100;