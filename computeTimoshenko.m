% This file is used to find the BC's for a Timoshenko beam (Sandwich)

function beam = computeTimoshenko(model)

plotcurves = model.plotcurves;

names = fieldnames(model);
for i=1:numel(names)
    eval([names{i},'=model.',names{i},';'])
end

x=0:L/1e3:L; %Convert scalar length to vector form

switch loadcase
    case 'simple_pt'
        
        a = xP;
        b = 1-xP;
        
        %Shear force distribution
        Tx = P*(1-xP)*ones(size(x));
        Tx(x>xP) = -P*xP;
        
        %Bending moment distribution
        Mx = P*b*x;
        Mx(x>xP) = P*a*(L-x(x>xP));

        
        %Twist angle distribution
        
        tb(x<xP) = -b*P*((b^2-1)*L^2 + 3*x(x<xP).^2)/(6*D);
        tb(x==xP) = 0;
        tb(x>xP) = a*P*((a^2+2)*L^2 -6*L*x(x>xP) + 3*x(x>xP).^2)/(6*D);   
        
        ts(x<xP) = P*b/S;
        ts(x==xP) = 0;
        ts(x>xP) = -P*a/S;
        
        ttot=tb+ts;
        
        %t = P*L^3*b/(6*D)*( (1-b^2)/L -3*x.^2./(L^3) ) + P*b/S;
        %t(x>xP) = P*L^3*a/(6*D)*((a^2-1)/L + 3*(L-x(x>xP)).^2./(L^3)) -P*a/S;
        
        %Deflections due to bending
        wb = P*L^3*b/(6*D)*((1-b^2)*(x./L)-(x./L).^3);
        wb(x>xP) = P*L^3*a/(6*D)*((1-a^2)*(L-x(x>xP))./L-((L-x(x>xP))./L).^3);
        
        %Deflections due to shear (Timoshenko)
        ws = P*b*x/S;
        ws(x>xP) = P*a*(L-x(x>xP))/S; 
        %ws = Tx*b/S;
        
        %Total deflection
        w = wb+ws;
        %w = ws;
        w(x>xP) = wb(x>xP) + ws(x>xP);
        
        
    case 'simple_dist'
        %Assuming that the input load is in Newton, the distribution is
        %created as 
        q = P/L; %[Force/length]
        
        %Reaction forces
        RL = q*L/2;
        RR = RL;
        
        %Shear force distribution
        Tx = RL-q*x;
        
        %Bending moment distribution
        Mx = RL*x-0.5*q*x.^2;
        
        %Twist angle distribution
        tb = q*L^4/(24*D)*( 4*x.^3/(L^4) -6.*x.^2/(L^3) + 1/L);
        ts = q/(2*S)*(L-2.*x);
        ttot = tb + ts;
        
        %Deflections due to bending
        wb = q*L^4/(24*D)*((x/L).^4 -2*(x/L).^3 + (x/L));% + q/(2*S).*(L*x-x.^2);
        %Deflections due to shear (Timoshenko)
        ws = q/(2*S)*(L*x-x.^2);
        
        %Total deflection
        w = wb + ws;
        
    case 'cantilever_pt'
        %Force input is in Newtons, located at the free end of the
        %cantilever beam.
        
        b = xP;
        a = 1-b;
        eta = x/L;

        %Shear force distribution
        Tx = -P*ones(size(x));
        Tx(x>xP) = 0;
        
        %Bending moment distribution
        Mx = Tx.*(xP*L-x);
        
        %Twist angle distribution
        ttot = P*L^2/(2*D)*(b^2-(eta-a).^2);
        ttot(eta>a) = P*L^2/(2*D)*b^2;
      
        %Deflections due to bending
        wb = P*L^3/(6*D)*(-b^3 + 3*b^2*(1-eta));
        wb(eta>a) = P*L^3/(6*D)*((eta(eta>a)-a).^3-3*b^2*(eta(eta>a)-a)+2*b^3);
        wb = flip(wb);
        
        %Deflections due to shear (Timoshenko)
        ws = Tx.*x/S;
        ws(x>xP) = ws(x==xP);
        
        %Total deflections
        w = wb + ws;

    case 'cantilever_dist'

        %Force input in Newton. 
        Q = P/L;
        
        %Shear force distribution
        %Tx = Q*x;
        Tx = Q.*(x-L);
        
        %Moment distribution
        %Mx = -Q*0.5*x.^2;
        Mx = -Q*L/(2*L)*(L-x).^2;
        
        %Twist angle distribution due to bending
        tb = Q*L^4/(24*D)*( 4*x.^3/(L^4) -12*x.^2/(L^3) +12*x/(L^2) );
        
        %Twist angle distribution due to shear
        ts = Tx./S;
        
        %Total angle distribution
        ttot = tb + ts;
        
        %Deflections due to bending
        wb = Q*L^4/(24*D)*( (x/L).^4 -4*(x/L).^3 +6*(x/L).^2 );
        
        %wb = Q*L^3/(24*D)*(eta.^4-4*eta+3);
        
        %Deflections due to shear
        ws = Q*x./(2*S).*(2*L-x);
        %ws = Tx.*x/S;
        
        %Total deflection
        w = wb + ws;
        %w = flip(w);
        
    otherwise
        warning('case not defined')
end
    
switch plotcurves
    case 'on'
        subplot(2,2,1)
        cla reset
        hold on
        plot(x,w,'b','linewidth',1)
        plot(x,wb,'r--','linewidth',1)
        plot(x,ws,'b:','linewidth',1)
        xlabel('x'),ylabel('w')
        grid on
        title('Displacement')
        legend('Total','Bending','Shear')
        subplot(2,2,2)
        cla reset
        hold on
        plot(x,ttot,'b','linewidth',1)
        plot(x,tb,'r--','linewidth',1)
        plot(x,ts,'b:','linewidth',1)
        xlabel('x'),ylabel('\theta [rad]')
        grid on
        title('Rotation')
        legend('Total','Bending','Shear')
        subplot(2,2,3)
        cla reset
        plot(x,Mx,'b','linewidth',1)
        xlabel('x'),ylabel('M')
        grid on
        title('Bending moment')
        subplot(2,2,4)
        axis tight
        plot(x,Tx,'b','linewidth',1)
        xlabel('x'),ylabel('Q')
        grid on
        title('Shear force')
        
end

if exist('f')
    beam=struct('x',x,'w',w,'t',ttot,'M',Mx,'Q',Tx,'f',f);
else
    beam=struct('x',x,'w',w,'t',ttot,'M',Mx,'Q',Tx);
end

end