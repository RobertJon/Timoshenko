% This file is used to find the BC's for a Timoshenko beam (Sandwich)

function beam = computeTimoshenko(model,plotcurves)

if nargin < 2
    plotcurves = 0;
end

names = fieldnames(model);
for i=1:numel(names)
    eval([names{i},'=model.',names{i},';'])
end

x=0:L/1e3:L; %Convert scalar length to vector form

switch loadcase
    case 'simple_pt'
        
        a = xP;
        b = 1-xP;
        
%         Q = -P*(L - xP)/L*ones(size(x));
%         Q(x>xP) = P*xP/L;
%         M = -(P*x*(L - xP))/L;
%         M(x>xP) = -(P*xP*(L - x(x>xP)))/L;
%         t = -(P*(L - xP)*(3*x.^2 + xP^2 - 2*L*xP))/(6*E*I*L);
%         t(x>xP) = (P*xP*(2*L^2 - 6*L*x(x>xP) + 3*x(x>xP).^2 + xP^2))/(6*E*I*L);
%         w = -(P*x.*(L - xP).*(x.^2 + xP^2 - 2*L*xP))/(6*E*I*L);
%         w(x>xP) = -(P*xP*(L - x(x>xP)).*(x(x>xP).^2 - 2*L*x(x>xP) + xP^2))/(6*E*I*L);
%         beta = (1:5)*pi/L;
%         f=beta.^2*sqrt(model.E*model.I/(model.rho*model.A))/2/pi;
        
        %Shear force distribution
        Tx = P*(1-xP)*ones(size(x));
        Tx(x>xP) = -P*xP;
        
        %Bending moment distribution
        Mx = P*(1-xP)*x;
        Mx(x>xP) = P*xP*(L-x(x>xP));
        
        %Twist angle distribution
        t = P*L^3*b/(6*D)*( (1-b^2)/L -3*x.^2./(L^3) ) + P*b/S;
        t(x>xP) = P*L^3*a/(6*D)*( (a^2-1)/L + 3*(L-x(x>xP)).^2./(L^3)) -P*a/S;
        
        %Deflections due to bending
        wb = P*L^3*b/(6*D)*((1-b^2)*(x./L)-(x./L).^3);
        wb(x>xP) = P*L^3*a/(6*D)*((1-a^2)*(L-x(x>xP))./L-((L-x(x>xP))./L).^3);
        
        %Deflections due to shear (Timoshenko)
        ws = P*b*x/S;
        ws(x>xP) = P*a*(L-x(x>xP))/S; 
        
        %Total deflection
        w = wb+ws;
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
        t = q*L^4/(24*D)*( 4*x.^3/(L^4) -6.*x.^2/(L^3) + 1/L) + q/(2*S)*(L-2.*x);
        
        %Deflections due to bending
        wb = q*L^4/(24*D)*((x/L).^4 -2*(x/L).^3 + (x/L)) + q/(2*S).*(L*x-x.^2);
        
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
        t = P*L^2/(2*D)*(b^2-(eta-a).^2);
        t(eta>a) = P*L^2/(2*D)*b^2;
      
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
        
        eta = x/L;
        
        %Force input in Newton. 
        Q = P/L;
        
        %Shear force distribution
        Tx = Q*x;
        
        %Moment distribution
        Mx = -Q*0.5*x.^2;
        
        %Twist angle distribution
        t = Q*L^2/(6*D)*(1-eta.^3);
        t = flip(t);
        
        %Deflections due to bending
        wb = Q*L^3/(24*D)*(eta.^4-4*eta+3);
        
        %Deflections due to shear
        ws = Tx.*x/S;
        
        %Total deflection
        w = wb + ws;
        w = flip(w);
        
    otherwise
        warning('case not defined')
end
    
if plotcurves
    subplot(2,2,1)
    plot(x,w,'b','linewidth',1)
    xlabel('x'),ylabel('w')
    grid on
    subplot(2,2,2)
    plot(x,t,'b','linewidth',1)
    xlabel('x'),ylabel('\theta [rad]')
    grid on
    subplot(2,2,3)
    plot(x,Mx,'b','linewidth',1)
    xlabel('x'),ylabel('M')
    grid on
    subplot(2,2,4)
    plot(x,Tx,'b','linewidth',1)
    xlabel('x'),ylabel('Q')
    grid on
end

if exist('f')
    beam=struct('x',x,'w',w,'t',t,'M',Mx,'Q',Tx,'f',f);
else
    beam=struct('x',x,'w',w,'t',t,'M',Mx,'Q',Tx);
end

end