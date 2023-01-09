%31 is entering CF2 and state 34 is obtained after throttling state 30 %

%    state=[h(kj/kg),temp(C),pressure(bar)]    %

state=[
    3390,540,177.5;
    3059.3,339.96,43.723;
    3537.1,540,39.351;
    3362.1,452.1,21.724;
    3165.6,353.50,10.674;
    2984.9,261.12,5.063;
    2746.9,136.76,1.393;
    2639.8,90.43,0.713;
    2515,68.5,.292;
    2382.2,45.24,0.097;
    189.4,45.24,0.097;
    191.2,45.68,0.951;
    194.0,46.36,20;
    271.1,64.33,19.4;
    361.9165,86.08,18.9;
    369.9378,88,18.4;
    439.9882,104.64,18.7;
    630.6,149.4,17;
    770,181.4,20;
    792.9,184.41,204;
    909.2,211.02,204;
    1095.1,251.83,204;
    1.1171e+03,256.49,200;
    939.1767,219.02,25;
    818.516,192.41,17.5;
    249.7747,59.64,1.7;
    2.6748e+03,94.08,.2;
    2.6652e+03,89.08,.2;
    395.5183,94.08,18.8;
    2.6250e+03,67.33,.11;
    2.6752e+03,99.4,0.95;
    2.6839e+03,98.21,.11;
    2.6270e+03,68.3,0.1;
    2.9280e+03,452.51,21.073;
    2.6250e+03,67.33,.01;
    ];

for i=1:34
    state(i,1)=XSteam('h_pT',state(i,3),state(i,2));
end


state(35,1)=XSteam('h_pT',state(30,3),state(30,2));

% entropy for each state using XSteam function kJ/(kg.celsius) %
entropy=[
         XSteam('s_ph',state(1,3),state(1,1));
         XSteam('s_ph',state(2,3),state(2,1));
         XSteam('s_ph',state(3,3),state(3,1));
         XSteam('s_ph',state(4,3),state(4,1));
         XSteam('s_ph',state(5,3),state(5,1));
         XSteam('s_ph',state(6,3),state(6,1));
         XSteam('s_ph',state(7,3),state(7,1));
         XSteam('s_ph',state(8,3),state(8,1));
         XSteam('s_ph',state(9,3),state(9,1));
         XSteam('s_ph',state(10,3),state(10,1));
         XSteam('s_ph',state(11,3),state(11,1));
         XSteam('s_ph',state(12,3),state(12,1));
         XSteam('s_ph',state(13,3),state(13,1));
         XSteam('s_ph',state(14,3),state(14,1));
         XSteam('s_ph',state(15,3),state(15,1));
         XSteam('s_ph',state(16,3),state(16,1));
         XSteam('s_ph',state(17,3),state(17,1));
         XSteam('s_ph',state(18,3),state(18,1));
         XSteam('s_ph',state(19,3),state(19,1));
         XSteam('s_ph',state(20,3),state(20,1));
         XSteam('s_ph',state(21,3),state(21,1));
         XSteam('s_ph',state(22,3),state(22,1));
         XSteam('s_ph',state(23,3),state(23,1));
         XSteam('s_ph',state(24,3),state(24,1));
         XSteam('s_ph',state(26,3),state(26,1));
         XSteam('s_ph',state(26,3),state(26,1));
         XSteam('s_ph',state(27,3),state(27,1));
         XSteam('s_ph',state(28,3),state(28,1));
         XSteam('s_ph',state(29,3),state(29,1));
         XSteam('s_ph',state(30,3),state(30,1));
         XSteam('s_ph',state(31,3),state(31,1));
         XSteam('s_ph',state(32,3),state(32,1));
         XSteam('s_ph',state(33,3),state(33,1));
         XSteam('s_ph',state(34,3),state(34,1));
         XSteam('s_ph',state(35,3),state(35,1));
         ];

% exergy y values for all states 
y=rand(34,1);
for i=1:35
    y(i,1)=state(i,1)-(298.15*entropy(i,1));
end

%initializing mass for all 34 states, where all mass flow rates are taken
%in kilograms/second
k1=270.278;
k31=0.204;
%first law equations to find mass flow rates in the power plant %

syms k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12 k13 k14 k15 k16 k17 k18 k19 k20 k21 k22 k23 k24 k25 k26 k27 k28 k29 k30 k32 k33 k34 kx

eqn1 = (k1)*(state(1,1)-state(2,1))+(k3)*(state(3,1)-state(4,1))+(k3-k4)*(state(4,1)-state(5,1))+(k3-k4-k5)*(state(5,1)-state(6,1))+(kx)*(state(6,1)-state(7,1))+(kx-k7)*(state(7,1)-state(8,1))+(kx-k7-k8)*(state(8,1)-state(9,1))+(kx-k7-k8-k9)*(state(9,1)-state(10,1))-((XSteam('v_pT',state(11,3),state(11,2)))*(state(12,3)-state(11,3)))-((XSteam('v_pT',state(28,3),state(28,2)))*(state(29,3)-state(28,3)))-((XSteam('v_pT',state(19,3),state(19,2)))*(state(20,3)-state(19,3)))-330000000==0;

eqn2 = (k12*state(12,1))+(k31*state(31,1))==(k32*state(32,1))+(k30*state(30,1));

eqn3 = (k13*state(13,1))+(k9*state(9,1))==(k14*state(14,1))+(k30*state(30,1));

eqn4 = (k14*state(14,1))+(k8*state(8,1))+(k27*state(27,1))==(k28*state(28,1))+(k15*state(15,1));

eqn5 = (k7*state(7,1))+(k16*state(16,1))+(k26*state(26,1))==(k17*state(17,1))+(k27*state(27,1));

eqn6 = (k6*state(6,1))+(k17*state(17,1))==(k18*state(18,1))+(k26*state(26,1));

eqn7 = (k34*state(34,1))+(k20*state(20,1))+(k24*state(24,1))==(k21*state(21,1))+(k25*state(25,1));

eqn8 = (k21*state(21,1))+(k2*state(2,1))==(k24*state(24,1))+(k22*state(22,1));

eqn9 = (k22*state(22,1))+(k4*state(4,1))==(k34*state(34,1))+(k23*state(23,1));

eqn10= (k5*state(5,1))+(k25*state(25,1))+(k18*state(18,1))==(k19*state(19,1));

eqn11= (k32*state(32,1))+(k30*state(30,1))==(k33*state(33,1));

eqn12= (k29*state(29,1))+(k15*state(15,1))==(k16*state(16,1));

eqn13=k19==k1;

eqn14=k20==k1;

eqn15=k21==k1;

eqn16=k22==k1;

eqn17=k23==k1;

eqn18=k32==k31;

eqn19=kx==k3-k4-k5-k6;

eqn20=k10==kx-k7-k8-k9;

eqn21=k25==k24+k34;

eqn22=k19==k5+k18+k25;

eqn23=k16==k17;

eqn24=k17==k18;

eqn25=k6==k26;

eqn26=k29==k28;

eqn27=k11==k10+k33;

eqn28=k3==k1-k2;

eqn29=k28==k8+k27;

eqn30=k30==k9;

eqn31=k11==k12;

eqn32=k27==k7+k26;

eqn33=k33==k30+k32;


[A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4,eqn5, eqn6, eqn7, eqn8, eqn9, eqn10, eqn11, eqn12, eqn13, eqn14, eqn15, eqn16, eqn17, eqn18, eqn19, eqn20, eqn21, eqn22, eqn23, eqn24, eqn25, eqn26, eqn27, eqn28, eqn29, eqn30, eqn31, eqn32, eqn33], [k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29, k30, k32, k33, k34, kx]);
X = linsolve(A,B);

k2=23.619;
k3=246.659;
k4=13.638;
k5=12.1;
k6=15.887;
k7=6.42;
k8=7.65;
k9=6.36;
k10=183.597;
k11=190.358;
k12=190.283;
k13=190.283;
k14=190.283;
k15=190.283;
k16=220.919;
k17=220.919;
k18=220.919;
k19=270.278;
k20=270.278;
k21=270.278;
k22=270.278;
k23=270.278;
k24=23.619;
k25=37.258;
k26=16.236;
k27=22.986;
k28=7.647;
k29=30.633;
k30=6.547;
k32=.204;
k33=6.760;
k34=13.638;
kx=208.97;


state(33,1)=(k32*state(32,1)+k30*state(30,1))/k33;

entropy(33,1)=XSteam('s_ph',state(33,3),state(33,1));

y(33,1)=state(33,1)-(298.15*entropy(33,1));

% energy %
w_hp=(k1)*(state(1,1)-state(2,1))

w_ip=(k3)*(state(3,1))-k4*(state(4,1))-k5*(state(5,1))-(kx+k6)*(state(6,1))

w_lp=(kx)*(state(6,1)-state(7,1))+(kx-k7)*(state(7,1)-state(8,1))+(kx-k7-k8)*(state(8,1)-state(9,1))+(kx-k7-k8-k9)*(state(9,1)-state(10,1))


% second law exergetic efficiencies %
n_hp=w_hp/(k1*(y(1,1)-y(2,1)))

n_ip=w_ip/(k3*(y(3,1)-y(4,1))+(k3-k4)*(y(4,1)-y(5,1))+(k3-k4-k5)*(y(5,1)-y(6,1)))

n_lp=w_lp/(kx*(y(6,1)-y(7,1))+(kx-k7)*(y(7,1)-y(8,1))+(kx-k7-k8)*(y(8,1)-y(9,1))+(kx-k7-k8-k9)*(y(9,1)-y(10,1)))

n_dj1=(k1*(state(1,1)-state(2,1)))/(k1*(y(1,1)-y(2,1)))

n_dj2=todecimal(((k3)*(state(3,1)-state(4,1))+(k3-k4)*(state(4,1)-state(5,1))+(k3-k4-k5)*(state(5,1)-state(6,1)))/(k3*(y(3,1)-y(4,1))+(k3-k4)*(y(4,1)-y(5,1))+(k3-k4-k5)*(y(5,1)-y(6,1))))

n_dj3=todecimal((kx)*(state(6,1)-state(7,1))+(kx-k7)*(state(7,1)-state(8,1))+(kx-k7-k8)*(state(8,1)-state(9,1))+(kx-k7-k8-k9)*(state(9,1)-state(10,1))/(kx)*(y(6,1)-y(7,1))+(kx-k7)*(y(7,1)-y(8,1))+(kx-k7-k8)*(y(8,1)-y(9,1))+(kx-k7-k8-k9)*(y(9,1)-y(10,1)))

n_dj4=(k18*(y(18,1)-y(19,1)))/(k6*(y(6,1)-y(26,1)))

n_mx1=(k33*(y(31,1)-y(32,1)))/(k30*(y(33,1)-y(35,1)))

n_mx2=(k15*(y(16,1)-y(15,1)))/(k29*(y(29,1)-y(26,1)))

n_deaerator=(k18*(y(19,1)-y(18,1)))/(k5*(y(5,1)-y(19,1))+k25*(y(25,1)-y(19,1)))

n_gj6=(k20*(y(21,1)-y(20,1)))/(k34*(y(34,1)-y(25,1))+k24*(y(24,1)-y(25,1)))

n_gj7=(k21*(y(22,1)-y(21,1)))/(k2*(y(2,1)-y(24,1)))

n_cf2=(k12*(y(13,1)-y(12,1)))/(k31*(y(31,1)-y(32,1)))

n_ex1=(k22*(y(23,1)-y(22,1)))/(k4*(y(4,1)-y(34,1)))
