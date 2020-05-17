function sol = model_SIR()
% MODEL_SIR simulation of the SIR model with n classes

% Dynamical parameters
n = 3;
N = popAge([25 70]);    % population
%mu = 3.4e-5;    % natural death rate 1/80 years = 1/80/365 per day 
%Lambda = N*mu;  % birth rate: set to keep N constant (ignoring deaths from covid-19)
gam = [1/7.0 ; 1/10.0 ; 1/14.0];  % recovery rate 1/two weeks = 1/14 
%R0 = 2.5;

%beta = eye(n).*R0.*gam;
%beta([2 3 4 6 7 8]) = [0.3 0.1 0.3 0.2 0.1 0.2];

% Contact matrix import
names = {'Age', '0-4y', '5-9y', '10-14y', '15-19y', '20-24y', '25-29y', '30-34y', '35-39y', '40-44y', '45-49y', '50-54y', '55-59y', '60-64y', '65-69y', '70y+'};
contact = readtable('./contact.txt');
contact.Properties.VariableNames = names;
contact = removevars(contact, {'Age'});
contact.Properties.RowNames = names;


% Integration parameters 
I0 = [0 ; 7 ; 0];
IC = [(N(1)-I0(1)) I0(1) 0 ; (N(2)-I0(2)) I0(2) 0 ; (N(3)-I0(3)) I0(3) 0 ]; % seed a few infected individuals
tspan = [46 365]; % in days 

% simulations
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
sol = ode45(@sir,tspan,IC,options);
soly = [sum(sol.y(1:3,:)) ;
        sum(sol.y(4:6,:)) ;
        sum(sol.y(7:9,:)) ];

figure(2); clf;
semilogy(sol.x, soly(1:3,:));
legend('S','I','R');

disp(['Percentage of contaminated 100*(E+I+R)/N at day ' num2str(tspan(2)), ': ', ...
  num2str(100*sum(sol.y(2:3,end)/N))]);


    function dxdt = sir(~,xx)  % nested function 
        % SEIR is the ODE rhs 
        S = xx(1:3);
        I = xx(4:6); 
        R = xx(7:9);
        
        
        dxdt = [ (- beta * I .* S./N);
                 (beta * I .* S./N - gam .* I);
                 (gam .* I) ]; 
             
    end                         % end nested function sir

end                             % end main function run_SIR 