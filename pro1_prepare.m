data = xlsread('附件3-弹性模量与压力.xlsx','A2:B402');
p=data(:,1);
E=data(:,2);
%拟合弹性模量与压力的关系
[fitresult, ~] = createFit(p, E);
rhospan = [0.85 0.9];
p0 = 100;
[rho,p] = ode15s(@(rho,x) (fitresult.p1*x^6 + fitresult.p2*x^5 + fitresult.p3*x^4 +...
    fitresult.p4*x^3 + fitresult.p5*x^2 + fitresult.p6*x + fitresult.p7)/rho, rhospan, p0);

%拟合密度和压力的关系
 [fitresult2, ~] = createFit1(rho, p);

pspan=[100 160];
rho0=0.85;
[p1 rhogao] = ode15s(@(x,rho) rho/(fitresult.p1*x^6 + fitresult.p2*x^5 + fitresult.p3*x^4 +...
    fitresult.p4*x^3 + fitresult.p5*x^2 + fitresult.p6*x + fitresult.p7),pspan, rho0);

rhogao = rhogao(length(rhogao));%高压侧密度

pspan=[100 150];
rho0=0.85;
[p1 rho_150] = ode15s(@(x,rho) rho/(fitresult.p1*x^6 + fitresult.p2*x^5 + fitresult.p3*x^4 +...
    fitresult.p4*x^3 + fitresult.p5*x^2 + fitresult.p6*x + fitresult.p7),pspan, rho0);
rho_150 = rho_150(length(rho_150));
