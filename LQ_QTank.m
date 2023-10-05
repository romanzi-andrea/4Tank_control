clc
clear all
close all

addpath utils;
load('Data/pump_params.mat');
load('Data/valve_opening.mat');

S = (4.44 / 2)^2 * pi;
y_target = [10; 10];

sys = struct('A', AA, 'S', S, 'g', 981, 'k1_fcn', k1_fcn, 'k2_fcn', k2_fcn, ...
                'gamma1_fcn', gamma1_fcn, 'gamma2_fcn', gamma2_fcn, ...
                'u_min', 4, 'u_max', 12, 'x_min', 0.1, 'x_max', 30);

[linsys, x_bar, u_bar] = linearize_4tanks(sys, y_target);

h0 = [ 2, 2, 2, 2 ].';

% Approximate the maps -> They are implemented via Lookup tables

xv = 0:0.05:18;
yk1 = k1_fcn(xv);
yk2 = k2_fcn(xv);
ygamma1 = gamma1_fcn(xv);
ygamma2 = gamma2_fcn(xv);

A=linsys.A;
B=linsys.B;
C=linsys.C;
D=linsys.D;
%% sys enlargement
%conditions

n = size(A);
n = n(1);
p = size(C);
p = p(1);
m = size(B);
m = m(2);

% p>=i p=2 so we can bring to 0 2 errors
zero(tf(ss(A,B,C(1,:),D(1,:))))    %no zeros in s=0
zero(tf(ss(A,B,C(2,:),D(2,:))))    %no zeros in s=0
i = 2;

A_tilde = [A, zeros(4,2); 
           -C, zeros(2)];
B_tilde = [B; zeros(2)];
C_tilde = [C, zeros(2)];
D_tilde = D;

en_sys = ss(A_tilde,B_tilde,C_tilde,D_tilde);
%% LQ control
%conditions 
Q = diag([10,10,10,10,.8,.7]);
R = 5*eye(m);
Cq = sqrt(Q);
rank(obsv(A_tilde,Cq))   %rank 6 ok
rank(ctrb(A_tilde, B_tilde))   %rank 6 ok
a = 0;
rank(ctrb(a*eye(n+i)+A_tilde, B_tilde))
rank(obsv(a*eye(n+i)+A_tilde,Cq)) 
K_lq = lqr(a*eye(n+i)+A_tilde,B_tilde, Q, R);
eig(A_tilde-B_tilde*K_lq)
tzero(ss(A_tilde-B_tilde*K_lq,B_tilde, C_tilde, D_tilde))

sigma(tf(ss(A,B,C,D)))
hold on
sigma(tf(ss(A_tilde-B_tilde*K_lq, B_tilde, C_tilde, D_tilde)))
legend OL_sv CL_sv
grid on

