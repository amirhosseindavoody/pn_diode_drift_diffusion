% This code simulates the P-N diode in the thermal equilibrium state.
% 12/6/2010
% Amirhossein Davoody

% This code simulates the P-N diode under bias.
% 12/6/2010
% Amirhossein Davoody

Na = 10^16; %cm^-3
Nd = 10^17; %cm^-3
epss = 1.05*10^-12; %silicon permitivity F/cm
TL = 300; %K
kB = 8.617*10^-5; %eV/K
q = 1.602*10^-19; %coulombos
ni = 1.5*10^10;
x_min = 0;
x_max = 5*10^-4;
Nmax = max(Na, Nd);
VT = kB*TL; %Thermal Voltage [eV]
LDi = (epss*VT/(q*ni))^0.5; %intrinsic Debye length
LDe = (epss*VT/(q*Nmax))^0.5; %extrinsic Debye length
LDa = (epss*VT/(q*Na))^0.5;
LDd = (epss*VT/(q*Nd))^0.5;
delta = 0.9; %mesh size in normalized scale
D = 1; %diffusion coefficient cm^2/s
M = D/VT; %mobility scaling factor
Ti = LDi^2/D; %time
dx = delta*LDe; %actual mesh size in real space
s = 0; %flag for terminating the iterative solving loop
ss = 0; %flag for terminating the iterative solving loop
er = 2; %flag for terminating the iterative solving loop
err = 2;
%Vp = V(1);
%Vn = V(n);
V = x_min:dx:x_max; %Potential matrix
V = V.';
dx = dx/LDi;
X = V; %Space matrix
C = V; %Doping matrix
f = C; %Forcing matrix
y = V; % L*y = f
[n,m] = size(V); % n = number of total space meshes
L = eye(n); %Lower triangular matrix
U = eye(n); %Upper triangular matrix
A = zeros(n); % A*V = f
E = V; %electric field
NA = V;
ND = V;
mun = V;
mup = V;
G = V;
%--------------------Arora Model Constants
MU1N = 88;
MU1P = 54;
MU2N = 1252;
MU2P = 407;
ALPHAN = -0.57;
ALPHAP = -0.57;
BETAN = -2.33;
BETAP = -2.33;
GAMMAN = 2.546;
GAMMAP = 2.546;
NCRITN = 1.432*10^17;
NCRITP = 2.67*10^17;
BETANN = 2;
BETAPP = 1;
ALPHANN = 2.4*10^7;
ALPHAPP = 2.4*10^7;
THETAN = 0.8;
THETAP = 0.8;
TNOMN = 600;
TNOMP = 600;
VSATN = ALPHANN/(1+THETAN*exp(TL/TNOMN));
VSATP = ALPHAPP/(1+THETAP*exp(TL/TNOMP));
%--------------------Generation Rate Constants
TAUN0 = 1*10^-7;
TAUP0 = 1*10^-7;
NSRHN = 5.0*10^16;
NSRHP = 5.0*10^16;
AN = 1;
AP = 1;
BN = 1;
BP = 1;
CN = 0;
CP = 0;
EN = 0;
EP = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Solution for Equilibrium %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off;
%--------------------Setting C Matrix Elements and Initial Guess for
%--------------------Potential
for i = 1:n  
    if(i < n/2)    
        V(i,1) = -log(Na/ni);
        NA(i,1) = Na/ni;
        ND(i,1) = 0;
        C(i,1) = -Na/ni;
    else
        V(i,1) = log(Nd/ni);
        NA(i,1) = Na/ni;
        ND(i,1) = Nd/ni;
        C(i,1) = Nd/ni;
    end;
end;

%--------------------Iteration Loop for Possion's Equation Solver
while ((s <10)&&(er > 10^-5))
    
    Vold = V;
    
    %--------------------Setting Forcing Vector
    for i = 1:n
        if(i == 1)
            f(i,1) = -log(Na/ni);
        elseif (i < n)
            f(i,1) = -(exp(-V(i,1))-exp(V(i,1))+C(i,1))-(exp(-V(i,1))+exp(V(i,1)))*V(i,1);
        else
            f(i,1) = log(Nd/ni);
        end;
    end;
    
    %--------------------Setting A Matrix Elements
    for i = 1:n
        if(i == 1)
            A(i,i) = 1;
        elseif (i == n)
            A(i,i) = 1;
            A(i,i-1) = 0;
        else
            A(i,i) = -(2/dx^2+exp(-V(i,1))+exp(V(i,1)));
            A(i,i-1) = 1/dx^2;
            A(i,i+1) = 1/dx^2;
        end;
    end;
    
    V = LU_decomposition(A,f);
    
    %--------------------Calculating Error
    er = 0;
    for i = 1:n
        er = max(er,abs(V(i,1)-Vold(i,1)));
    end
    s = s+1;
end


N = exp(V);
P = exp (-V);
semilogy(X,N,'k');
%plot(X,V,'k');
hold on;
semilogy(X,P,'k');
Veq = V;
Neq = N;
Peq = P;
s = 0;
er = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Solution for Non-Equilibrium %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------Arora Model for Mobility
mun0_n = MU1N*(TL/300)^ALPHAN+(MU2N*(TL/300)^BETAN)/(1+Nd/(NCRITN*(TL/300)^GAMMAN));
mun0_p = MU1N*(TL/300)^ALPHAN+(MU2N*(TL/300)^BETAN)/(1+Na/(NCRITN*(TL/300)^GAMMAN));
mup0_n = MU1P*(TL/300)^ALPHAP+(MU2P*(TL/300)^BETAP)/(1+Nd/(NCRITP*(TL/300)^GAMMAP));
mup0_p = MU1P*(TL/300)^ALPHAP+(MU2P*(TL/300)^BETAP)/(1+Na/(NCRITP*(TL/300)^GAMMAP));

%--------------------Carrier Lifetime for Generation
tn_n = TAUN0/(AN+BN*Nd/NSRHN+CN*(Nd/NSRHN)^EN);
tn_p = TAUN0/(AN+BN*Na/NSRHN+CN*(Na/NSRHN)^EN);
tp_n = TAUP0/(AP+BP*Nd/NSRHP+CP*(Nd/NSRHP)^EN);
tp_p = TAUP0/(AP+BP*Na/NSRHP+CP*(Na/NSRHP)^EN);

%--------------------Diode Bias Voltage
Vb = 0.5/VT;
nvstep = 20;
for kk = 0:nvstep
    V(1) = Veq(1)+kk/nvstep*0.5/VT;
    V(n) = Veq(n);
    Vold = V*0;
    %--------------------Iterative Loop for Non-Equilibrium Solution
    while (err > 10^-4)
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Possion's Equation Solver%%%%%%%%%%%%%%%%%%
        %--------------------Setting Forcing Vector Elements
        for i = 1:n
            if ((i == 1) || (i == n))
                f(i) = V(i);
            else
                f(i) = N(i)-P(i)-C(i)-V(i)*(N(i)+P(i));
            end;
        end;
        
        %--------------------Setting A Matrix Elements
        A = A*0;
        for i = 1:n
            if ((i == 1) || (i == n))
                A(i,i) = 1;
            else
                A(i,i) = -(2/dx^2)-(N(i)+P(i));
                A(i,i-1) = 1/dx^2;
                A(i,i+1) = 1/dx^2;
            end;
        end;
        
        V = LU_decomposition (A,f);
        
        for i = 1:n-1
            if(i<n/2)
                mun(i) = mun0_p*(1/(1+(mun0_p*abs(V(i)-V(i+1))*VT/(VSATN*LDi*dx))^BETANN))^1/BETANN;
                mup(i) = mup0_p*(1/(1+(mup0_p*abs(V(i)-V(i+1))*VT/(VSATP*LDi*dx))^BETAPP))^1/BETAPP;
            else
                mun(i) = mun0_n*(1/(1+(mun0_n*abs(V(i)-V(i+1))*VT/(VSATN*LDi*dx))^BETANN))^1/BETANN;
                mup(i) = mup0_n*(1/(1+(mup0_n*abs(V(i)-V(i+1))*VT/(VSATP*LDi*dx))^BETAPP))^1/BETAPP;
            end;
        end;
        mun(n) = 0;
        mup(n) = 0;
        
        d_n = VT*mun;
        d_p = VT*mup;
        
        %%%%%%%%%%%%%%%%%%Electron Continiuty Equation Solver%%%%%%%%%%%%%%%%%%
        while (er>1)
            Nold = N;
            %--------------------Calculating Normalized Generation Rate
            for i = 1:n
                if(i < n/2)
                    G(i) = LDi^2*(1-N(i)*P(i))/(tp_p*(N(i)+1)+tn_p*(P(i)+1));
                else
                    G(i) = LDi^2*(1-N(i)*P(i))/(tp_n*(N(i)+1)+tn_n*(P(i)+1));
                end;
            end;
            
            %--------------------Setting Elements of Forcing Vector
            for i = 1:n
                if((i == 1)||(i == n))
                    f(i) = N(i);
                else
                    f(i) = G(i);
                end;
            end;
            
            %--------------------Setting Elements of A Matrix
            A = A*0;
            for i = 1:n
                if ((i == 1)||(i == n))
                    A(i,i) = 1;
                else
                    A(i,i-1) = 1/dx^2*d_n(i-1)*Bernouli(V(i-1)-V(i));
                    A(i,i) = -1/dx^2*(d_n(i-1)*Bernouli(V(i)-V(i-1))+d_n(i)*Bernouli(V(i)-V(i+1)));
                    A(i,i+1) = 1/dx^2*d_n(i)*Bernouli(V(i+1)-V(i));
                end;
            end;
            
            N = LU_decomposition(A,f);
            %semilogy(X,ni*N,'b');
            er = 0;
            for i = 1:n
                er = max(er,abs(N(i,1)-Nold(i,1)));
            end;
            s = s+1;
        end;
        
        er = 2;
        s = 0;
        
        %%%%%%%%%%%%%%%%%%Hole Continiuty Equation Solver%%%%%%%%%%%%%%%%%%%%%%
        while (er>1)
            Pold = P;
            %--------------------Calculating Normalized Generation Rate
            for i = 1:n
                if(i < n/2)
                    G(i) = LDi^2*(1-N(i)*P(i))/(tp_p*(N(i)+1)+tn_p*(P(i)+1));
                else
                    G(i) = LDi^2*(1-N(i)*P(i))/(tp_n*(N(i)+1)+tn_n*(P(i)+1));
                end;
            end;
            
            %--------------------Setting Elements of Forcing Vector
            for i = 1:n
                if((i == 1)||(i == n))
                    f(i) = P(i);
                else
                    f(i) = G(i);
                end;
            end;
            
            %--------------------Setting Elements of A Matrix
            A = A*0;
            for i = 1:n
                if ((i == 1)||(i == n))
                    A(i,i) = 1;
                else
                    A(i,i-1) = 1/dx^2*d_p(i-1)*Bernouli(V(i)-V(i-1));
                    A(i,i) = -1/dx^2*(d_p(i-1)*Bernouli(V(i-1)-V(i))+d_p(i)*Bernouli(V(i+1)-V(i)));
                    A(i,i+1) = 1/dx^2*d_p(i)*Bernouli(V(i)-V(i+1));
                end;
            end;
            
            P = LU_decomposition(A,f);
            
            er = 0;
            for i = 1:n
                er = max(er,abs(P(i,1)-Pold(i,1)));
            end;
            ss = ss+1;
        end;
        
        er = 2;
        s = 0;
        
        err = 0;
        for i = 1:n
            err = max(err,abs(V(i,1)-Vold(i,1)));
        end;
        s = s+1;
        Vold = V;
    end;
    semilogy(X,N,'b');
    semilogy(X,P,'b');
    %plot(X,V,'b');
    for i = 1:n
        NN(i,kk+21) = N(i);
        PP(i,kk+21) = P(i);
        VV(i,kk+21) = V(i);
    end;
    err = 2;
end;

V = Veq;
N = Neq;
P = Peq

for kk = 0:nvstep
    V(1) = Veq(1)-kk/nvstep*0.5/VT;
    V(n) = Veq(n);
    Vold = V*0;
    %--------------------Iterative Loop for Non-Equilibrium Solution
    while (err > 10^-4)
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Possion's Equation Solver%%%%%%%%%%%%%%%%%%
        %--------------------Setting Forcing Vector Elements
        for i = 1:n
            if ((i == 1) || (i == n))
                f(i) = V(i);
            else
                f(i) = N(i)-P(i)-C(i)-V(i)*(N(i)+P(i));
            end;
        end;
        
        %--------------------Setting A Matrix Elements
        A = A*0;
        for i = 1:n
            if ((i == 1) || (i == n))
                A(i,i) = 1;
            else
                A(i,i) = -(2/dx^2)-(N(i)+P(i));
                A(i,i-1) = 1/dx^2;
                A(i,i+1) = 1/dx^2;
            end;
        end;
        
        V = LU_decomposition (A,f);
        
        for i = 1:n-1
            if(i<n/2)
                mun(i) = mun0_p*(1/(1+(mun0_p*abs(V(i)-V(i+1))*VT/(VSATN*LDi*dx))^BETANN))^1/BETANN;
                mup(i) = mup0_p*(1/(1+(mup0_p*abs(V(i)-V(i+1))*VT/(VSATP*LDi*dx))^BETAPP))^1/BETAPP;
            else
                mun(i) = mun0_n*(1/(1+(mun0_n*abs(V(i)-V(i+1))*VT/(VSATN*LDi*dx))^BETANN))^1/BETANN;
                mup(i) = mup0_n*(1/(1+(mup0_n*abs(V(i)-V(i+1))*VT/(VSATP*LDi*dx))^BETAPP))^1/BETAPP;
            end;
        end;
        mun(n) = 0;
        mup(n) = 0;
        
        d_n = VT*mun;
        d_p = VT*mup;
        
        %%%%%%%%%%%%%%%%%%Electron Continiuty Equation Solver%%%%%%%%%%%%%%%%%%
        while (er>1)
            Nold = N;
            %--------------------Calculating Normalized Generation Rate
            for i = 1:n
                if(i < n/2)
                    G(i) = LDi^2*(1-N(i)*P(i))/(tp_p*(N(i)+1)+tn_p*(P(i)+1));
                else
                    G(i) = LDi^2*(1-N(i)*P(i))/(tp_n*(N(i)+1)+tn_n*(P(i)+1));
                end;
            end;
            
            %--------------------Setting Elements of Forcing Vector
            for i = 1:n
                if((i == 1)||(i == n))
                    f(i) = N(i);
                else
                    f(i) = G(i);
                end;
            end;
            
            %--------------------Setting Elements of A Matrix
            A = A*0;
            for i = 1:n
                if ((i == 1)||(i == n))
                    A(i,i) = 1;
                else
                    A(i,i-1) = 1/dx^2*d_n(i-1)*Bernouli(V(i-1)-V(i));
                    A(i,i) = -1/dx^2*(d_n(i-1)*Bernouli(V(i)-V(i-1))+d_n(i)*Bernouli(V(i)-V(i+1)));
                    A(i,i+1) = 1/dx^2*d_n(i)*Bernouli(V(i+1)-V(i));
                end;
            end;
            
            N = LU_decomposition(A,f);
            %semilogy(X,ni*N,'b');
            er = 0;
            for i = 1:n
                er = max(er,abs(N(i,1)-Nold(i,1)));
            end;
            s = s+1;
        end;
        
        er = 2;
        s = 0;
        
        %%%%%%%%%%%%%%%%%%Hole Continiuty Equation Solver%%%%%%%%%%%%%%%%%%%%%%
        while (er>1)
            Pold = P;
            %--------------------Calculating Normalized Generation Rate
            for i = 1:n
                if(i < n/2)
                    G(i) = LDi^2*(1-N(i)*P(i))/(tp_p*(N(i)+1)+tn_p*(P(i)+1));
                else
                    G(i) = LDi^2*(1-N(i)*P(i))/(tp_n*(N(i)+1)+tn_n*(P(i)+1));
                end;
            end;
            
            %--------------------Setting Elements of Forcing Vector
            for i = 1:n
                if((i == 1)||(i == n))
                    f(i) = P(i);
                else
                    f(i) = G(i);
                end;
            end;
            
            %--------------------Setting Elements of A Matrix
            A = A*0;
            for i = 1:n
                if ((i == 1)||(i == n))
                    A(i,i) = 1;
                else
                    A(i,i-1) = 1/dx^2*d_p(i-1)*Bernouli(V(i)-V(i-1));
                    A(i,i) = -1/dx^2*(d_p(i-1)*Bernouli(V(i-1)-V(i))+d_p(i)*Bernouli(V(i+1)-V(i)));
                    A(i,i+1) = 1/dx^2*d_p(i)*Bernouli(V(i)-V(i+1));
                end;
            end;
            
            P = LU_decomposition(A,f);
            
            er = 0;
            for i = 1:n
                er = max(er,abs(P(i,1)-Pold(i,1)));
            end;
            ss = ss+1;
        end;
        
        er = 2;
        s = 0;
        
        err = 0;
        for i = 1:n
            err = max(err,abs(V(i,1)-Vold(i,1)));
        end;
        s = s+1;
        Vold = V;
    end;
    semilogy(X,N,'b');
    semilogy(X,P,'b');
    %plot(X,V,'b');
    for i = 1:n
        NN(i,21-kk) = N(i);
        PP(i,21-kk) = P(i);
        VV(i,21-kk) = V(i);
    end;
    err = 2;
end;



