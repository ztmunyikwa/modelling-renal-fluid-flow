% function called by ode solver -- u(1) = Dsa, u(2) = Asa, u(3) = Dla, u(4) = Ala
% DAall2
function z = arterioleODE(t,u,Pin,Qnum,Pdrop,multiplier,capdensity,CMD_DA,td,Pac,...
    Pend,MyoOn_dist,TGFOn_dist,ndistc,mu_dist,Leff_dist,Qdistc,Pdistc,Params_dist,CLc,ta)
%CMD_DA in mol/cm^3

% load FiveCompConsMRSRCR2.mat Pac Pend MyoOn_dist TGFOn_dist
% load CStateValsMRSRCR2.mat

% td = 1;
% ta = 4.8;

Cdistold = (pi*(u(1)^4)*ndistc)/(128*mu_dist*Leff_dist);

Qden = 1/Cdistold;
Qdistnew = Qdistc*((Pin-Pend)/(Pac-Pend))*(Qnum/Qden);
Qdistnewtest = (Pin-Pend)/(Qden);  %this is equivalent to the previous line

Pdistnew = 0.5*(Pin + Pend); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Key for autoregulation curves!   - next section, just pressure remains normal in each case!!
        
if MyoOn_dist == 1
    Pstartdist = Pdistnew;
else
    Pstartdist = Pdistc;  %control state assumption
%     Pstartdist = Pdistnew;  %completely zeroed out
%     Params_dist(6) = 0;     %completely zeroed out
end

if TGFOn_dist == 1
    CMD_DA = CMD_DA;
else
    CMD_DA = CLc;   %control state assumption
%     Params_dist(7) = 0;  %completely zero out
end


%z = [dDdt; dAdt];
% %%main one I'm using.
z = [1/td*2/Pdistc*(Pstartdist*u(1)/2 - ...
    (Params_dist(1)*exp(Params_dist(2)*((u(1)/Params_dist(10)) - 1)) + ...
    (Params_dist(3)*exp(-1*(((u(1)/Params_dist(10))-Params_dist(4))/Params_dist(5))^2))*u(2)));...
    1/ta*(1/(1+exp(-(Params_dist(6)*Pstartdist*u(1)/2 + Params_dist(7)/(1+exp(-Params_dist(8)*(CMD_DA-CLc))) - Params_dist(9))))-u(2))];
% save DAallvals2

