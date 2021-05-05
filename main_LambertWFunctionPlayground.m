clc;
close all;
clear all;

syms x y
f = lambertw(2,x + 1i*y);
interval = [-100 100 -100 100];
fmesh(real(f),interval,'ShowContours','On')