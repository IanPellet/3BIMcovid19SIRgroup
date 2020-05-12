function sol = run_SIR()
% RUN_SIR simulation of the SIR model

% Dynamical parameters 
N = 6.0e8;      % total population
mu = 3.4e-5;    % natural death rate 1/80 years = 1/80/365 per day 
Lambda = N*mu;  % birth rate: set to keep N constant (ignoring deaths from covid-19)
gam = 1/14.0;  % recovery rate 1/two weeks = 1/14 

% transmission rate based on R0 = 2.5
% R0 = beta/gam
% beta = R0*gam
R0 = 2.5;
beta = R0*gam;   % infection rate S -> I

% Integration parameters 
IC = [(N-1e4); 1e4; 0]; % seed a few infected individuals
tspan = [0 365]; % in days 

% simulations
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
sol = ode45(@sir,tspan,IC,options)

figure(2); clf;
semilogy(sol.x, sol.y(1:3,:));

disp(['Percentage of contaminated 100*(E+I+R)/N at day ' num2str(tspan(2)), ': ', ...
  num2str(100*sum(sol.y(2:3,end)/N))]);


    function dxdt = sir(~,xx)  % nested function 
        % SEIR is the ODE rhs 
        
        S = xx(1); I = xx(2); R = xx(3);
        
        dxdt = [ Lambda - mu*S - beta*I/N*S;
                 beta*I/N*S - (mu+gam)*I;
                 gam*I - mu*R ]; 
             
    end                         % end nested function sir

end                             % end main function run_SIR 
