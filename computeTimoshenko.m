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
        

    case 'cantilever_ptend'
        %Force input is in Newtons, located at the free end of the
        %cantilever beam.
        
        a = xP;
        b = 1-a;
        eta = x/L;
        
        %Shear force distribution
        Tx = P*ones(size(x));
        
        %Bending moment distribution
        Mx = -P.*(L-x);
        
        %Twist angle distribution
        t = P*L^3/(2*D)*(b^2-(eta-a).^2);
        t(eta>xP) = P*L^2*b^2/(2*D);
        
        %Deflections due to bending
        wb = P*L^3/(3*D)*( (eta-a).^3 -3*b^2*(eta-a) +2*b^3);
        wb(x>xP) = P*L^3/(3*D)*(-b^3 + 3*b^2*(1-eta(eta>xP)));
        
        %Deflections due to shear (Timoshenko)
        ws = Tx.*x/S;
        
        %Total deflections
        w = wb + ws;
        
        t= ones(size(x));
        
    case 'cantilever_dist'
        %Force input in Newton. 
        
        
        
    otherwise
        warning('case not defined')
end
    
if plotcurves
    subplot(2,2,1)
    plot(x,w)
    xlabel('x'),ylabel('w')
    subplot(2,2,2)
    plot(x,t)
    xlabel('x'),ylabel('\theta [rad]')
    subplot(2,2,3)
    plot(x,Mx)
    xlabel('x'),ylabel('M')
    subplot(2,2,4)
    plot(x,Tx)
    xlabel('x'),ylabel('Q')
end

if exist('f')
    beam=struct('x',x,'w',w,'t',t,'M',Mx,'Q',Tx,'f',f);
else
    beam=struct('x',x,'w',w,'t',t,'M',Mx,'Q',Tx);
end

end