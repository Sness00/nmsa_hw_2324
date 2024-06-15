%% Data for the following PDE
% ro(x)*u_{tt} -  (mu(x)*u_{x})_x = f  in (a,b) x (0,T] 
% u(0) = u0;
% v(0) = v0;
% + b.c in {a,b} x (0,T]
function [Data] = DataTest(TestName)

if strcmp(TestName,'Q3a')
    Data.name = TestName;
    
    % Spatial domain
    Data.domain = [0, 1.5];
    
    % Time domain and discretization
    Data.T = 5;
    Data.dt = 0.01;
    
    % Parameters
    Data.ro = @(x) 1*(x <= 0.75) + 4*(x > 0.75);
    Data.mu = @(x) 4*(x <= 0.75) + 1*(x > 0.75);    
    Data.force = @(x,t) 0.*x + 0.*t;
    
    % B.C.
    Data.boundary = 'DN';
    Data.gD1  =  @(t) 4*pi*sqrt(exp(1))*(t-0.3)*exp(-8*pi^2*(t-0.3)^2);
    Data.gN2  =  @(t) 0.*t;
    
    % Initial conditions
    Data.u0 = @(x) 0.*x;
    Data.v0 = @(x) 0.*x;

elseif strcmp(TestName, 'Q3b')
    Data.name = TestName;

     % Spatial domain
    Data.domain = [0, 1.5];
    
    % Time domain and discretization
    Data.T = 5;
    Data.dt = 0.01;
    
    % Parameters
    Data.ro = @(x) (x<0.5) + (1 + 6*(x-0.5))*(x>=0.5 && x<=1) + 4*(x>1);
    Data.mu = @(x) 4*(x<0.5) + (4 + 3*(x-0.5))*(x>=0.5 && x<=1) + (x>1);    
    Data.force = @(x,t) 0.*x + 0.*t;
    
    % B.C.
    Data.boundary = 'DN';
    Data.gD1  =  @(t) 4*pi*sqrt(exp(1))*(t-0.3)*exp(-8*pi^2*(t-0.3)^2);
    Data.gN2  =  @(t) 0.*t;

    % I.C.
    Data.u0 = @(x) 0.*x;
    Data.v0 = @(x) 0.*x;

elseif strcmp(TestName, 'Q3c')
    Data.name = TestName;

    Data.name = TestName;

     % Spatial domain
    Data.domain = [0, 1.5];
    
    % Time domain and discretization
    Data.T = 5;
    Data.dt = 0.01;

    % Parameters
    Data.ro = @(x) 1 + 2*x; 
    Data.mu = @(x) 4 - 2*x;
    Data.force = @(x,t) 0.*x + 0.*t;
    
   % B.C.
    Data.boundary = 'DN';
    Data.gD1  =  @(t) 4*pi*sqrt(exp(1))*(t-0.3)*exp(-8*pi^2*(t-0.3)^2);
    Data.gN2  =  @(t) 0.*t;
    
    % I.C.
    Data.u0 = @(x) 0.*x;
    Data.v0 = @(x) 0.*x;

elseif strcmp(TestName, 'test')
    Data.name = TestName;
    
    % Exact solution
    Data.uex = @(x, t) sin(2*pi*x).*sin(2*pi*Data.c(x).*t);
    Data.uex_x = @(x, t) pi*t.*(Data.mu_x(x).*Data.ro(x) - Data.mu(x).*Data.ro_x(x)) ... 
        ./Data.ro(x).^2./Data.c(x).*sin(2*pi*x).*cos(2*pi*Data.c(x).*t) + ...
                         2*pi*cos(2*pi*x).*sin(2*pi*Data.c(x).*t);

    % Domain
    Data.domain = [0, 1];
    
    % Time
    Data.T = 1;
    Data.dt = 0.001;
    
    % Parameters
    Data.ro = @(x) 16 + 7.8.*x; 
    Data.mu = @(x) 21 + 3.3.*x;

    Data.ro_x = @(x) 7.8 + 0.*x;
    Data.mu_x = @(x) 3.3 + 0.*x;

    Data.c = @(x) sqrt(Data.mu(x)./Data.ro(x));    
    
    A = @(x, t) 4*pi^2*t./Data.c(x) .* (Data.mu_x(x).*Data.ro(x) - Data.mu(x).*Data.ro_x(x))./Data.ro(x).^2;
    B = @(x, t) -(pi*t./Data.c(x).*(Data.mu_x(x).*Data.ro(x) - Data.mu(x).*Data.ro_x(x))./Data.ro(x).^2).^2;
    C = @(x, t) -pi*t.*(2./Data.c(x).*Data.ro_x(x).*(Data.mu_x(x).*Data.ro(x) - Data.mu(x).*Data.ro_x(x))./Data.ro(x).^3 + ...
                        Data.c(x)/2.*((Data.mu_x(x).*Data.ro(x) - Data.mu(x).*Data.ro_x(x))./(Data.mu(x).*Data.ro(x))).^2);
    
    Data.force = @(x, t) -Data.mu_x(x).*Data.uex_x(x, t) - Data.mu(x).*(A(x, t).*cos(2*pi*x).*cos(2*pi*Data.c(x).*t) + ...
                                                                       B(x, t).*sin(2*pi*x).*sin(2*pi*Data.c(x).*t) + ...
                                                                       C(x, t).*sin(2*pi*x).*cos(2*pi*Data.c(x).*t));
    % B.C.
    Data.boundary = 'DD';    
    Data.gD1  =  @(t) 0.*t;
    Data.gD2  =  @(t) 0.*t;
    
    % I.C.
    Data.u0 = @(x) 0.*x;
    Data.v0 = @(x) 2*pi*Data.c(x).*sin(2*pi*x);
end