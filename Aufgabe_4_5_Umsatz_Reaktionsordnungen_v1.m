clearvars; close all; clc
%% Aufgabe 4.5
%% a.) 
%{
    Reatkion 1: 1. Ordnung, da ln(c_A) entspricht der Linearisierung
    A -> B

    Reatkion 2: 2. Ordnung, da 1/c_A entspricht der Linearisierung
    2A -> B

    Reatkion 3: 2. Ordnung, da 1/c_A entspricht der Linearisierung
    A + B -> C + 2D
    c_A = c_B
    
    Reaktion 3: 2. Ordnung
    A + B -> C + 2D
    c_A != c_B
%}
%% c.)
%Reaktion 1
c_A_0 = 1.25; %[mol/L]
c_A_ln = -1.8; %[ln(mol/L)]
t = 400; %[s]

k_1 = (c_A_ln - log(c_A_0)) / -t %[1/s]
%% 
%Reaktion 2
c_A_0 = 1.25; %[mol/L]
c_A_U1 = 2.8; %[L/mol]
t_2 = 200; %[s]

k_2 = (c_A_U1 - (1 / c_A_0)) / (2 * t_2) %[(L / mol) * (1 / s)]
%% 
%Reaktion 3.1
c_A_0 = 1.25; %[mol/L]
c_A_U2 = 2.8; %[L/mol]
t = 400; %[s]

k_2_AB = (c_A_U2 - (1 / c_A_0)) / t %[(L / mol) * (1 / s)]
%%
%Reaktion 3.2
c_A_0 = 1.25; %[mol/L]
c_B_0 = 1.5; %[mol/L]
AB_01 = -0.5; %[L/mol]
t = 400; %[s]

k_2_AB2 = AB_01 / ((c_A_0 - c_B_0) * t) %[(L / mol) * (1 / s)]
