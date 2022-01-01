function dYdt = odefun7(t,Y,w,beta1,beta2,h_beta3,T1)
   %两个都开启的状态
    dYdt = zeros(2,1);
    %高压油泵的密度变化
    t1 = roundn(mod(t,T1),-2);
    t2 = roundn(mod(t,100),-2);


   dYdt(1) = 1/(20+(7.239-((beta2.p1*(w*t1)^4+beta2.p2*(w*t1)^3+ beta2.p3*(w*t1)^2 +...
        beta2.p4*(w*t1)+beta2.p5)))*pi*2.5^2)*...
        (pi*2.5^2*Y(1)*...
        (4*beta2.p1*w^4*t1^3+3*beta2.p2*w^3*t1^2+2*beta2.p3*w^2*t1 + beta2.p4*w));
    %油管的密度变化
    dYdt(2) = -0.85*pi/(pi*5^2*500)*(((7.892+(h_beta3.p1*t2^4+h_beta3.p2*t2^3+h_beta3.p3*t2^2+h_beta3.p4*t2+h_beta3.p5))*tand(9))^2-...
        1.25^2)*sqrt(2*Y(2)*...
        ((beta1.p1*Y(2)^4 +beta1.p2*Y(2)^3 +beta1.p3*Y(2)^2 +beta1.p4*Y(2) + beta1.p5)-...
        0.101));
    
    
end
