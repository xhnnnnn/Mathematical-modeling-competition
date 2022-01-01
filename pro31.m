clear
data = xlsread('附件3-弹性模量与压力.xlsx','A2:B402');
p=data(:,1);
E=data(:,2);
%拟合弹性模量与压力的关系
[fitresult, gof] = createFit(p, E);
%求密度从0.85到0.9之间的压力的数值解

rhospan = [0.85 0.9];
p0 = 100;
[rho,p] = ode15s(@(rho,x) (fitresult.p1*x^6 + fitresult.p2*x^5 + fitresult.p3*x^4 +...
    fitresult.p4*x^3 + fitresult.p5*x^2 + fitresult.p6*x + fitresult.p7)/rho, rhospan, p0);
%用求出的数值解拟合压力与密度的关系
[beta1, ~] = createFit1(rho, p);

data = xlsread('附件1-凸轮边缘曲线.xlsx','A2:B629');
L = data(:,2);
theta = data(:,1);
%拟合极径与极角的关系
 [beta2, ~] = createFit2(theta, L);

%上升段的距离随时间变化
yundong1 = xlsread('附件2-针阀运动曲线.xlsx','A2:B46');
penyouzui_t1 = yundong1(:,1);
penyouzui_h1 = yundong1(:,2);
%拟合上升段距离时间变化的多项式
[h_beta1,~] = juli_t1(penyouzui_t1,penyouzui_h1);
%下降段的距离随时间变化
yundong3 = xlsread('附件2-针阀运动曲线.xlsx','D2:E46');
penyouzui_t3 = yundong3(:,1);
penyouzui_h3 = yundong3(:,2);
[h_beta3,~] = juli_t3(penyouzui_t3,penyouzui_h3);

Lmax = max(L);%初始位置的极经 上止点
Lmin = min(L);%下止点

%得到上止点液体体积和下止点液体体积
 V_min = 20;
 V_max = 20 + pi*2.5^2*(Lmax-Lmin);
 Pmin = 0.5;
 
 %求0.5压力下的密度 也就是下止点的密度
pspan = [100 0.5];
rho0 = 0.85;
[p1,rho_min] = ode15s(@(x,rho) rho/(fitresult.p1*x^6 + fitresult.p2*x^5 + fitresult.p3*x^4 +...
    fitresult.p4*x^3 + fitresult.p5*x^2 + fitresult.p6*x + fitresult.p7),pspan, rho0);
 rho_min = rho_min(length(rho_min));%0.5MPA下的密度

m = rho_min*V_max;%高压油泵的质量
S_youbeng = pi*2.5^2;%高压油泵的面积

for t_penyoujiu = 2.45:0.3:97.55 
for w = 0.01:0.01:1
jishi = 0;%计时器

total_P2 = [];%管道压力
total_rho2 = [];    %管道密度
    
T1 = roundn(6.27/w,-2);%活塞周期
T2 = 100;

rho2_pre = 0.85;
i=1;

tic

 while jishi < 1000
    
     T1_shengyu = roundn(mod(jishi,T1+10),-2);
     T2_shengyu = roundn(mod(jishi,T2),-2);
     
    if mod(roundn(T1_shengyu,-2),T1+10)==0
        T1_shengyu=0;
    end
    if mod(roundn(T2_shengyu,-2),T2)==0
        T2_shengyu=0;
    end
    t_shengyu = 0;
     
    if T1_shengyu < roundn(T1/2,-2) && T2_shengyu < roundn(0.33,-2)
        t_shengyu = min(roundn(T1/2,-2)-T1_shengyu,0.33-T2_shengyu);
        tspan = roundn(jishi:0.01:jishi+t_shengyu,-2);
        tspan = mod(tspan,T2);
        [t1,rho] = ode15s(@(t,x) -0.85*pi/(pi*5^2*500)*(((7.892+(h_beta1.p1*t^4+h_beta1.p2*t^3+...
            h_beta1.p3*t^2+h_beta1.p4*t+h_beta1.p5))*tand(9))^2-1.25^2)*sqrt(2*x*...
        ((beta1.p1*x^4 +beta1.p2*x^3 +beta1.p3*x^2 +beta1.p4*x + beta1.p5)-0.101)),tspan,rho2_pre);
        rho(1) =[];
        total_rho2 = [total_rho2;rho];
        rho2_pre = rho(length(rho));%更新当前管道密度
        jishi = jishi+t_shengyu;
        T1_shengyu = roundn(mod(jishi,T1),-2);
        T2_shengyu = roundn(mod(jishi,T2),-2);
        if (roundn(T1/2,-2)-T1_shengyu) == min(roundn(T1/2,-2)-T1_shengyu,2.12-T2_shengyu)
            rho1_pre = rho_min;
            while t_pengyouzui<100
            end
        end
    elseif T1_shengyu < roundn(T1/2,-2) && T2_shengyu < roundn(2.12,-2)
        t_shengyu = min(roundn(T1/2,-2)-T1_shengyu,2.12-T2_shengyu);
        tspan = roundn(jishi:0.01:jishi+t_shengyu,-2);
        tspan = mod(tspan,T2);
        
        [t1,rho] = ode15s(@(t,x) -0.85*pi/(pi*5^2*500)*0.7^2*sqrt(2*x*...
        ((beta1.p1*x^4 +beta1.p2*x^3 +beta1.p3*x^2 +beta1.p4*x + beta1.p5)-0.101)),tspan,rho2_pre);
        rho(1) =[];
        total_rho2 = [total_rho2;rho];
        rho2_pre = rho(length(rho));%更新当前管道密度
        jishi = jishi+t_shengyu;
        T1_shengyu = roundn(mod(jishi,T1),-2);
        T2_shengyu = roundn(mod(jishi,T2),-2);
        if (roundn(T1/2,-2)-T1_shengyu) == min(roundn(T1/2,-2)-T1_shengyu,2.46-T2_shengyu)
            rho1_pre = rho_min;
        end

    elseif T1_shengyu < roundn(T1/2,-2) && T2_shengyu < roundn(2.46,-2)
        t_shengyu = min(roundn(T1/2,-2)-T1_shengyu,2.46-T2_shengyu);
        tspan = roundn(jishi:0.01:jishi+t_shengyu,-2);
        tspan = mod(tspan,T2);
        [t1,rho] = ode15s(@(t,x) -0.85*pi/(pi*5^2*500)*(((7.892+(h_beta3.p1*t^4+h_beta3.p2*t^3+h_beta3.p3*t^2+...
            h_beta3.p4*t+h_beta3.p5))*tand(9))^2-1.25^2)*sqrt(2*x*...
            ((beta1.p1*x^4 +beta1.p2*x^3 +beta1.p3*x^2 +beta1.p4*x + beta1.p5)-0.101)),tspan,rho2_pre);
        rho(1) =[];
        total_rho2 = [total_rho2;rho];
        rho2_pre = rho(length(rho));%更新当前管道密度
        jishi = jishi+t_shengyu;
        T1_shengyu = roundn(mod(jishi,T1),-2);
        T2_shengyu = roundn(mod(jishi,T2),-2);
        if (roundn(T1/2,-2)-T1_shengyu) == min(roundn(T1/2,-2)-T1_shengyu,100-T2_shengyu)
            rho1_pre = rho_min;
        end
        
     elseif T1_shengyu < roundn(T1/2,-2) && T2_shengyu < roundn(T2,-2)
        t_shengyu = min(roundn(T1/2,-2)-T1_shengyu,100-T2_shengyu);
        tspan = roundn(jishi:0.01:jishi+t_shengyu,-2);
        tspan = mod(tspan,T2);
        tspan(tspan == 0) = T2;
        [t1,rho] = ode15s(@(t,x) 0,tspan,rho2_pre);
        rho(1) =[];       
        rho2_pre = rho(length(rho));%更新当前管道密度
        total_rho2 = [total_rho2;rho];

        jishi = jishi+t_shengyu;
        T1_shengyu = roundn(mod(jishi,T1),-2);
        T2_shengyu = roundn(mod(jishi,T2),-2);
        if (roundn(T1/2,-2)-T1_shengyu) == min(roundn(T1/2,-2)-T1_shengyu,100-T2_shengyu)
            rho1_pre = rho_min;
        end       
    elseif T1_shengyu <roundn(T1,-2)  && T2_shengyu < roundn(0.33,-2)   
        t_shengyu = min(roundn(T1,-2)-T1_shengyu,0.33-T2_shengyu);
        tspan = 0:0.01:0.01;
        tspan = roundn(tspan+jishi,-2);
        if roundn(rho1_pre,-3)>roundn(rho2_pre,-3)
            [t1,rho] = ode15s(@(t,x) odefun1(t,x,w,beta1,beta2,h_beta1,T1),tspan,[rho1_pre rho2_pre]);
        else
            [t1,rho] = ode15s(@(t,x) odefun5(t,x,w,beta1,beta2,h_beta1,T1),tspan,[rho1_pre rho2_pre]);
        end
        rho = real(rho);

        rho1 = rho(:,1);
        rho2 = rho(:,2);
        rho2(1) =[];
        rho2_pre = rho2(length(rho2));        
        total_rho2 = [total_rho2;rho2_pre];
        rho1_pre = rho1(length(rho1));%更新当前的密度
        jishi = jishi+0.01;
        
    elseif T1_shengyu <roundn(T1,-2)  && T2_shengyu < roundn(2.12,-2)
         t_shengyu = min(roundn(T1,-2)-T1_shengyu,2.12-T2_shengyu);
        tspan = 0:0.01:0.01;
        tspan = roundn(tspan+jishi,-2);

        if roundn(rho1_pre,-3)>roundn(rho2_pre,-3)
            [t1,rho] = ode15s(@(t,x) odefun2(t,x,w,beta1,beta2,T1),tspan,[rho1_pre rho2_pre]);
        else
             [t1,rho] = ode15s(@(t,x) odefun6(t,x,w,beta1,beta2,T1),tspan,[rho1_pre rho2_pre]);
        end
        rho = real(rho);
        
        rho1 = rho(:,1);
        rho2 = rho(:,2);
        rho2(1) =[];
        rho2_pre = rho2(length(rho2));
        total_rho2 = [total_rho2;rho2_pre];
        rho1_pre = rho1(length(rho1));%更新当前的密度
        jishi = jishi+0.01;
        
     elseif T1_shengyu <roundn(T1,-2) && T2_shengyu < roundn(2.46,-2)
         t_shengyu = min(roundn(T1,-2)-T1_shengyu,2.46-T2_shengyu);
        tspan = 0:0.01:0.01;
        tspan = roundn(tspan+jishi,-2);
        
        if roundn(rho1_pre,-3)>roundn(rho2_pre,-3)
            [t1,rho] = ode15s(@(t,x) odefun3(t,x,w,beta1,beta2,h_beta3,T1),tspan,[rho1_pre rho2_pre]);
        else
            [t1,rho] = ode15s(@(t,x) odefun7(t,x,w,beta1,beta2,h_beta3,T1),tspan,[rho1_pre rho2_pre]);
        
        end
        rho = real(rho);
        
        rho1 = rho(:,1);
        rho2 = rho(:,2);
        rho2(1) =[];
        rho2_pre = rho2(length(rho2));
        total_rho2 = [total_rho2;rho2_pre];
        rho1_pre = rho1(length(rho1));%更新当前的密度
        jishi = jishi+0.01;
        
     elseif T1_shengyu <roundn(T1,-2) && T2_shengyu < roundn(100,-2)
         t_shengyu = min(roundn(T1,-2)-T1_shengyu,100-T2_shengyu);
        tspan = 0:0.01:0.01;
        tspan = roundn(tspan+jishi,-2);

        if roundn(rho1_pre,-3)>roundn(rho2_pre,-3)
            [t1,rho] = ode15s(@(t,x) odefun4(t,x,w,beta1,beta2,T1),tspan,[rho1_pre rho2_pre]);
        else
            [t1,rho] = ode15s(@(t,x) odefun8(t,x,w,beta1,beta2,T1),tspan,[rho1_pre rho2_pre]);
        end
        rho = real(rho);
        
        rho1 = rho(:,1);
        rho2 = rho(:,2);
        rho2(1) =[];
        rho2_pre = rho2(length(rho2));
        total_rho2 = [total_rho2;rho2_pre];
        rho1_pre = rho1(length(rho1));%更新当前的密度
        jishi = jishi+0.01;
    elseif T1_shengyu < roundn(T1+10,-2) && T2_shengyu < roundn(0.33,-2)
        t_shengyu = min(roundn(T1+10,-2)-T1_shengyu,0.33-T2_shengyu);
        tspan = roundn(jishi:0.01:jishi+t_shengyu,-2);
        tspan = mod(tspan,T2);
        [t1,rho] = ode15s(@(t,x) -0.85*pi/(pi*5^2*500)*(((7.892+(h_beta1.p1*t^4+h_beta1.p2*t^3+...
            h_beta1.p3*t^2+h_beta1.p4*t+h_beta1.p5))*tand(9))^2-1.25^2)*sqrt(2*x*...
        ((beta1.p1*x^4 +beta1.p2*x^3 +beta1.p3*x^2 +beta1.p4*x + beta1.p5)-0.101)),tspan,rho2_pre);
        rho(1) =[];
        total_rho2 = [total_rho2;rho];
        rho2_pre = rho(length(rho));%更新当前管道密度
        jishi = jishi+t_shengyu;

    elseif T1_shengyu < roundn(T1+10,-2) && T2_shengyu < roundn(2.12,-2)
        t_shengyu = min(roundn(T1+10,-2)-T1_shengyu,2.12-T2_shengyu);
        tspan = roundn(jishi:0.01:jishi+t_shengyu,-2);
        tspan = mod(tspan,T2);
        
        [t1,rho] = ode15s(@(t,x) -0.85*pi/(pi*5^2*500)*0.7^2*sqrt(2*x*...
        ((beta1.p1*x^4 +beta1.p2*x^3 +beta1.p3*x^2 +beta1.p4*x + beta1.p5)-0.101)),tspan,rho2_pre);
        rho(1) =[];
        total_rho2 = [total_rho2;rho];
        rho2_pre = rho(length(rho));%更新当前管道密度
        jishi = jishi+t_shengyu;

    elseif T1_shengyu < roundn(T1+10,-2) && T2_shengyu < roundn(2.46,-2)
        t_shengyu = min(roundn(T1+10,-2)-T1_shengyu,2.46-T2_shengyu);
        tspan = roundn(jishi:0.01:jishi+t_shengyu,-2);
        tspan = mod(tspan,T2);
        [t1,rho] = ode15s(@(t,x) -0.85*pi/(pi*5^2*500)*(((7.892+(h_beta3.p1*t^4+h_beta3.p2*t^3+h_beta3.p3*t^2+...
            h_beta3.p4*t+h_beta3.p5))*tand(9))^2-1.25^2)*sqrt(2*x*...
            ((beta1.p1*x^4 +beta1.p2*x^3 +beta1.p3*x^2 +beta1.p4*x + beta1.p5)-0.101)),tspan,rho2_pre);
        rho(1) =[];
        total_rho2 = [total_rho2;rho];
        rho2_pre = rho(length(rho));%更新当前管道密度
        jishi = jishi+t_shengyu;
        T1_shengyu = roundn(mod(jishi,T1),-2);
        T2_shengyu = roundn(mod(jishi,T2),-2);
        if (roundn(T1/2,-2)-T1_shengyu) == min(roundn(T1/2,-2)-T1_shengyu,100-T2_shengyu)
            rho1_pre = rho_min;
        end
        
     elseif T1_shengyu < roundn(T1+10,-2) && T2_shengyu < roundn(T2,-2)
        t_shengyu = min(roundn(T1+10,-2)-T1_shengyu,100-T2_shengyu);
        tspan = roundn(jishi:0.01:jishi+t_shengyu,-2);
        tspan = mod(tspan,T2);
        tspan(tspan == 0) = T2;
        [t1,rho] = ode15s(@(t,x) 0,tspan,rho2_pre);
        rho(1) =[];       
        rho2_pre = rho(length(rho));%更新当前管道密度
        total_rho2 = [total_rho2;rho];

        jishi = jishi+t_shengyu;
        T1_shengyu = roundn(mod(jishi,T1),-2);
        T2_shengyu = roundn(mod(jishi,T2),-2);
        if (roundn(T1/2,-2)-T1_shengyu) == min(roundn(T1/2,-2)-T1_shengyu,100-T2_shengyu)
            rho1_pre = rho_min;
        end              
        
    end
    
 end
 
    toc
    total_P = beta1.p1*total_rho2.^4+beta1.p2*total_rho2.^3+beta1.p3*total_rho2.^2+beta1.p4*total_rho2+beta1.p5;
    all_P(i,:) = total_P(1:100000);
    pianyi_sum(i) = sum(abs(total_P-100));
    pianyi_max(i) = max(abs(total_P-100));
    i=i+1;
end
end 
[min_1 mindex] = min(all_P)