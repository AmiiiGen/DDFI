function [A, Bu, C, Du, K] = modelExample
%--modelExample(.)--/Generating an observable realization of 
%  state-space  matrices corresponding to a predefined discrete-time 
%  transfer function matrix.
%
%         [       ]
%  G(z) = |  TFs  |   to    x(k+1) = A*x(k) + B*u(k) + K*e(k)
%         [       ]           y(k) = C*x(k) + D*u(k) + e(k)  
%         
%  ----------------
%  @author: Amin Sheikhi
% 

%% Discrete-time System Matrices

% denominator in descending powers of z
a = poly([0.8 0.38 0.9 0.5]);    

% numerator in descending powers of z
b11 = poly([-0.59 0.95 0.4 -0.7 ]);
b21 = poly([-0.59 0.2 -0.45 0.63]);
b31 = poly([-0.59 0.95 -0.18 0.83]);

G11 = tf(b11,a,-1);
G21 = tf(b21,a,-1);
G31 = tf(b31,a,-1);

G = [G11; G21; G31];

[A,Bu,C,Du] = tf2ss([b11;b21;b31],a);
sys_ss = ss(A,Bu,C,Du,1);

sys_obs = compreal(sys_ss, 'o');

A_o = sys_obs.A;
Bu_o = sys_obs.B;
C_o = sys_obs.C;
Du_o = sys_obs.D;

n = size(A_o,2);
nu = size (Bu_o,2);
ny = size(C_o,1);

%--->>> process and measurement noises:
var_w = 0.25;       % process noise variance
var_v = 0.25;       % measurement noise variance
mu = 0;             % zero-mean noises
Qw = var_w*eye(n);
Rv = var_v*eye(ny);

% the innovation form
[P,K,~] = idare(A',C',Qw,Rv,[],[]);  
K = K';                              % steady-state Kalman gain 
cov_e = C*P*C' + Rv;                 % Covariance of innovation

end