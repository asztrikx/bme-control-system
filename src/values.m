clear;

m = 5;
M = 1.8;
J = 1.4e-4;
P = 0.04;
L0 = 1.5;
y0 = 0.9;

g = 9.81;
a = sqrt((m+M)*g / (m*L0));
omega0 = a/2;
xi = sqrt(2)/2;
scinf = -1.5*a;
soinf = -3*a;
Tdstart = 4/a;
Tdstart = 4;

dnorm = 1;
dnorm = 0.09/(0.3-0.09) * dnorm;

W1 = [
   m+M m*L0;
   m   m*L0
];
O1 = [
    0 0;
    0 -m*g
];
Q1 = [
    1;
    0
];
W2 = [m+J/(P^2)];
O2 = [0];
Q2 = [-1/P];

cA1 = [
    zeros(2,2), eye(2);
    inv(W1)*O1, zeros(2,2)
];
cB1 = [
    zeros(2,1);
    inv(W1)*Q1
];
cC1 = [1, zeros(1, 3)];
cD1 = 0;

cA2 = [
    0           1;
    inv(W2)*O2, 0
];
cB2 = [
    0;
    inv(W2)*Q2
];
cC2 = [1 0];
cD2 = 0;

sys1 = ss(cA1, cB1, cC1, cD1);
sys2 = ss(cA2, cB2, cC2, cD2);
T = 1/(5*abs(soinf));
Dsys1 = c2d(sys1, T, 'zoh');
Dsys2 = c2d(sys2, T, 'zoh');
[phi1, gamma1, C1, D1] = ssdata(Dsys1);
[phi2, gamma2, C2, D2] = ssdata(Dsys2);

disp("Stabilitás");
eig(phi1) % érintik a kört => labilis
eig(phi2) % érintik a kört => labilis

Mc1 = ctrb(phi1, gamma1);
assert(rank(Mc1) == 4);
Mc2 = ctrb(phi2, gamma2);
assert(rank(Mc2) == 2);

sDominant1 = -omega0*xi + omega0*sqrt(1-xi^2)*1i;
sDominant2 = conj(sDominant1);
zDominant1 = exp(sDominant1*T);
zDominant2 = exp(sDominant2*T);
zoinf = exp(soinf*T);
zcinf = exp(scinf*T);

K1 = acker(phi1, gamma1, [zDominant1, zDominant2, zcinf, zcinf]);
K2 = acker(phi2, gamma2, [zDominant1, zDominant2]);
disp("Zárt kör_1,2 kapott SÉ-i");
eig(phi1-gamma1*K1);
eig(phi2-gamma2*K2);
%initial(ss(phi1-gamma1*K1,[],[1 0 0 0],[], T), [10 0 0 0])
%initial(ss(phi2-gamma2*K2,[],[1 0    ],[], T), [10 0])

N1 = inv([
    phi1 - eye(4), gamma1;
    C1,            0
]) * [
    zeros(4,1);
    1
]; % [1; 0; 0; 0] as r is the first state variable
N2 = inv([
    phi2 - eye(2), gamma2;
    C2,              0
]) * [
    zeros(2,1);
    1
]; % [1; 0] as l ...
Nx1 = N1(1:4);
Nx2 = N2(1:2);
Nu1 = N1(end);
Nu2 = N2(end);

Mo1 = obsv(phi1, C1);
Mo2 = obsv(phi2, C2);
assert(rank(Mo1) == 4);
assert(rank(Mo2) == 2);

Gt1 = acker(phi1', phi1'*C1', [zoinf zoinf zoinf zoinf]); %phi1'*C1'
Gt2 = acker(phi2', phi2'*C2', [zoinf zoinf]);
G1 = Gt1';
G2 = Gt2';
F1 = phi1 - G1*C1*phi1;
F2 = phi2 - G2*C2*phi2;
H1 = gamma1 - G1*C1*gamma1;
H2 = gamma2 - G2*C2*gamma2;

phiT1 = [
    phi1,       gamma1;
    zeros(1,4), 1
];
phiT2 = [
    phi2,       gamma2;
    zeros(1,2), 1
];
gammaT1 = [
    gamma1;
    0
];
gammaT2 = [
    gamma2;
    0
];
CT1 = [C1 0];
CT2 = [C2 0];
GTt1 = acker(phiT1', phiT1'*CT1', [zoinf zoinf zoinf zoinf zoinf]);
GTt2 = acker(phiT2', phiT2'*CT2', [zoinf zoinf zoinf]);
GT1 = GTt1';
GT2 = GTt2';
FT1 = phiT1-GT1*CT1*phiT1;
FT2 = phiT2-GT2*CT2*phiT2;
HT1 = gammaT1-GT1*CT1*gammaT1;
HT2 = gammaT2-GT2*CT2*gammaT2;