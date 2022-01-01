clear
pro1_prepare

all_P = zeros(200,500000);
pianyi_max = zeros(1,200);
pianyi_sum = zeros(1,200);
i=1;
for T = 0.1:0.01:2
jishi = 0;%计时器
total_P=[];
rho0=rho_150;%燃油目前的密度
P0=150;
tic
while jishi<=5000
    t1 = mod(jishi,T+10);
    t2 = mod(jishi,100);
    if roundn(t1,-3)==T+10
        t1=0;
    end
    if roundn(t2,-3)==100
        t2=0;
    end
    t_shengyu = 0;

    if roundn(t1,-3)<roundn(T,-3)&&roundn(t2,-3)<roundn(0.2,-3)
        t_shengyu=min(T-t1,0.2-t2);%当前状态的剩余时
        tspan = 0:0.01:roundn(t_shengyu,-3);
        if  (fitresult2.p1*rho0^4 + fitresult2.p2*rho0^3 +...
                fitresult2.p3*rho0^2 +fitresult2.p4*rho0 + fitresult2.p5) < 160
            [t1,rho1] = ode15s(@(t,x) 0.85 * (pi*0.7^2)* sqrt(2*rhogao*(160-(fitresult2.p1*x^4 + fitresult2.p2*x^3 + fitresult2.p3*x^2 +...
            fitresult2.p4*x + fitresult2.p5)))/(500*pi*5^2)-100*t2*x/(500*pi*5^2), tspan, rho0);
        else
            [t1,rho1] = ode15s(@(t,x) -100*t2*x/(500*pi*5^2), tspan, rho0);
        end
        rho0 = rho1(length(rho1));%燃油目前的密度
        P1=fitresult2.p1*rho1.^4 + fitresult2.p2*rho1.^3 + fitresult2.p3*rho1.^2 +fitresult2.p4*rho1 + fitresult2.p5;
        P1(1)=[];
        total_P=[total_P,P1'];
        P0 = total_P(length(total_P));%燃油目前的压力
        jishi = jishi + t_shengyu;        
    elseif roundn(t1,-3)<roundn(T,-3)&&roundn(t2,-3)<roundn(2.2,-3)
        t_shengyu=min(T-t1,2.2-t2);%当前状态的剩余时
        tspan = 0:0.01:roundn(t_shengyu,-3);
        if  (fitresult2.p1*rho0^4 + fitresult2.p2*rho0^3 +...
                fitresult2.p3*rho0^2 +fitresult2.p4*rho0 + fitresult2.p5) < 160
            [t1,rho1] = ode15s(@(t,x) 0.85 * (pi*0.7^2)* sqrt(2*rhogao*(160-(fitresult2.p1*x^4 + fitresult2.p2*x^3 + fitresult2.p3*x^2 +...
            fitresult2.p4*x + fitresult2.p5)))/(500*pi*5^2)-20*x/(500*pi*5^2), tspan, rho0);
        else
            [t1,rho1] = ode15s(@(t,x) -20*x/(500*pi*5^2), tspan, rho0);
        end
        rho0 = rho1(length(rho1));%燃油目前的密度
        P1=fitresult2.p1*rho1.^4 + fitresult2.p2*rho1.^3 + fitresult2.p3*rho1.^2 +fitresult2.p4*rho1 + fitresult2.p5;
        P1(1)=[];
        total_P=[total_P,P1'];
        P0 = total_P(length(total_P));%燃油目前的压力
        jishi = jishi + t_shengyu;  
    elseif roundn(t1,-3)<roundn(T,-3)&&roundn(t2,-3)<roundn(2.4,-3)
        t_shengyu=min(T-t1,2.4-t2);%当前状态的剩余时
        tspan = 0:0.01:roundn(t_shengyu,-3);
        if  (fitresult2.p1*rho0^4 + fitresult2.p2*rho0^3 +...
                fitresult2.p3*rho0^2 +fitresult2.p4*rho0 + fitresult2.p5) < 160            
            [t1,rho1] = ode15s(@(t,x) 0.85 * (pi*0.7^2)* sqrt(2*rhogao*(160-(fitresult2.p1*x^4 + fitresult2.p2*x^3 + fitresult2.p3*x^2 +...
            fitresult2.p4*x + fitresult2.p5)))/(500*pi*5^2)-(480-20*t2)*x/(500*pi*5^2), tspan, rho0);
        else
            [t1,rho1] = ode15s(@(t,x) -(480-20*t2)*x/(500*pi*5^2), tspan, rho0);
        end
        rho0 = rho1(length(rho1));%燃油目前的密度
        P1=fitresult2.p1*rho1.^4 + fitresult2.p2*rho1.^3 + fitresult2.p3*rho1.^2 +fitresult2.p4*rho1 + fitresult2.p5;
        P1(1)=[];
        total_P=[total_P,P1'];
        P0 = total_P(length(total_P));%燃油目前的压力
        jishi = jishi + t_shengyu;  
    elseif roundn(t1,-3)<roundn(T,-3)&&roundn(t2,-3)<roundn(100,-3)
        t_shengyu=min(T-t1,100-t2);%当前状态的剩余时
        tspan = 0:0.01:roundn(t_shengyu,-3);
        if  (fitresult2.p1*rho0^4 + fitresult2.p2*rho0^3 +...
                fitresult2.p3*rho0^2 +fitresult2.p4*rho0 + fitresult2.p5) < 160            
            [t1,rho1] = ode15s(@(t,x) 0.85 * (pi*0.7^2)* sqrt(2*rhogao*(160-(fitresult2.p1*x^4 + fitresult2.p2*x^3 + fitresult2.p3*x^2 +...
            fitresult2.p4*x + fitresult2.p5)))/(500*pi*5^2), tspan, rho0);
        else
            [t1,rho1] = ode15s(@(t,x) 0, tspan, rho0);
        end
        rho0 = rho1(length(rho1));%燃油目前的密度
        P1=fitresult2.p1*rho1.^4 + fitresult2.p2*rho1.^3 + fitresult2.p3*rho1.^2 +fitresult2.p4*rho1 + fitresult2.p5;
        P1(1)=[];
        total_P=[total_P,P1'];
        P0 = total_P(length(total_P));%燃油目前的压力
        jishi = jishi + t_shengyu;  
    elseif roundn(t1,-3)<roundn(T+10,-3)&&roundn(t2,-3)<roundn(0.2,-3)
        t_shengyu=min(T+10-t1,0.2-t2);%当前状态的剩余时
        tspan = 0:0.01:roundn(t_shengyu,-3);
        [t1,rho1] = ode15s(@(t,x) -(100*t2)*x/(500*pi*5^2), tspan, rho0);
        rho0 = rho1(length(rho1));%燃油目前的密度
        P1=fitresult2.p1*rho1.^4 + fitresult2.p2*rho1.^3 + fitresult2.p3*rho1.^2 +fitresult2.p4*rho1 + fitresult2.p5;
        P1(1)=[];
        total_P=[total_P,P1'];
        P0 = total_P(length(total_P));%燃油目前的压力
        jishi = jishi + t_shengyu;  
    elseif roundn(t1,-3)<roundn(T+10,-3)&&roundn(t2,-3)<roundn(2.2,-3)
        t_shengyu=min(T+10-t1,2.2-t2);%当前状态的剩余时
        tspan = 0:0.01:roundn(t_shengyu,-3);
        [t1,rho1] = ode15s(@(t,x) -20*x/(500*pi*5^2), tspan, rho0);
        rho0 = rho1(length(rho1));%燃油目前的密度
        P1=fitresult2.p1*rho1.^4 + fitresult2.p2*rho1.^3 + fitresult2.p3*rho1.^2 +fitresult2.p4*rho1 + fitresult2.p5;
        P1(1)=[];
        total_P=[total_P,P1'];
        P0 = total_P(length(total_P));%燃油目前的压力
        jishi = jishi + t_shengyu;  
    elseif roundn(t1,-3)<roundn(T+10,-3)&&roundn(t2,-3)<roundn(2.4,-3)
        t_shengyu=min(T+10-t1,2.4-t2);%当前状态的剩余时
        tspan = 0:0.01:roundn(t_shengyu,-3);
        [t1,rho1] = ode15s(@(t,x) -(240-120*t2)*x/(500*pi*5^2), tspan, rho0);
        rho0 = rho1(length(rho1));%燃油目前的密度
        P1=fitresult2.p1*rho1.^4 + fitresult2.p2*rho1.^3 + fitresult2.p3*rho1.^2 +fitresult2.p4*rho1 + fitresult2.p5;
        P1(1)=[];
        total_P=[total_P,P1'];
        P0 = total_P(length(total_P));%燃油目前的压力
        jishi = jishi + t_shengyu;  
   elseif roundn(t1,-3)<roundn(T+10,-3)&&roundn(t2,-3)<roundn(100,-3)
        t_shengyu=min(T+10-t1,100-t2);%当前状态的剩余时
        tspan = 0:0.01:roundn(t_shengyu,-3);
        [t1,rho1] = ode15s(@(t,x) 0, tspan, rho0);
        rho0 = rho1(length(rho1));%燃油目前的密度
        P1=fitresult2.p1*rho1.^4 + fitresult2.p2*rho1.^3 + fitresult2.p3*rho1.^2 +fitresult2.p4*rho1 + fitresult2.p5;
        P1(1)=[];
        total_P=[total_P,P1'];
        P0 = total_P(length(total_P));%燃油目前的压力
        jishi = jishi + t_shengyu; 
    end
end

all_P(i,:)=total_P(1:500000);
pianyi_sum(i)=sum(abs(total_P(1:500000)-150));
pianyi_max(i)=max(abs(total_P(1:500000)-150));
i=i+1;
toc
end
