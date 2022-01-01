data = xlsread('����3-����ģ����ѹ��.xlsx','A2:B402');
p=data(:,1);
E=data(:,2);
%��ϵ���ģ����ѹ���Ĺ�ϵ
[fitresult, ~] = createFit(p, E);
rhospan = [0.85 0.9];
p0 = 100;
[rho,p] = ode15s(@(rho,x) (fitresult.p1*x^6 + fitresult.p2*x^5 + fitresult.p3*x^4 +...
    fitresult.p4*x^3 + fitresult.p5*x^2 + fitresult.p6*x + fitresult.p7)/rho, rhospan, p0);

%����ܶȺ�ѹ���Ĺ�ϵ
 [fitresult2, ~] = createFit1(rho, p);

pspan=[100 160];
rho0=0.85;
[p1 rhogao] = ode15s(@(x,rho) rho/(fitresult.p1*x^6 + fitresult.p2*x^5 + fitresult.p3*x^4 +...
    fitresult.p4*x^3 + fitresult.p5*x^2 + fitresult.p6*x + fitresult.p7),pspan, rho0);

rhogao = rhogao(length(rhogao));%��ѹ���ܶ�

pspan=[100 150];
rho0=0.85;
[p1 rho_150] = ode15s(@(x,rho) rho/(fitresult.p1*x^6 + fitresult.p2*x^5 + fitresult.p3*x^4 +...
    fitresult.p4*x^3 + fitresult.p5*x^2 + fitresult.p6*x + fitresult.p7),pspan, rho0);
rho_150 = rho_150(length(rho_150));
