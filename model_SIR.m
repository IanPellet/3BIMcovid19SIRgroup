% Saturday at 18:40
% modele SIR : Pauline , Maëlys, Nino , Ian, Alexandre Vi

function sol = model_SIR()
% MODEL_SIR simulation of the SIR model with n classes

% Dynamical parameters
n = 4;

N = popAge([18 40 70]);    % population

%gam = [1/7.0 ; 1/10.0 ; 1/14.0];  % recovery rate 1/two weeks = 1/14 
gam = [1/17.0 ; 1/.23 ; 1/35.0 ; 1/30.0];

R0 = 2.5;       %Avant confinement
%R0 = 0.5;       %Pendant confinement
%R0 = 0.9;       %Après confinement

%contact_de_base = 1.0;
%M = rand(n) * contact_de_base; % contact_de_base * 0.5 contact en moyenne
%M = M + diag(rand(n,1) + 3.0);
%M = M + diag(rand(n-2,1) + 3.0,2);
%M = M + diag(diag(M,2),-2);        contact aléatoire      taille de
%ppulation total variable?

C = [1 2 1.5 0 ; 2 4 3 1 ; 1.5 3 2.5 1 ; 0 1 1 0.5];        %contact possible avant confinement
%C = [0. 1.5 1 0 ; 1.5 2 1 0 ; 1 1 1 0.5 ; 0 0 0.5 0.5];        %contact possible pendant confinement
%C = [0. 1.5 1 0 ; 1.5 3 1.5 0 ; 1 1.5 2 0.5 ; 0 0 0.5 0.5];        %contact possible apres confinement

%beta = eye(n).*R0.*gam;
%beta([2 3 4 5 7 8 9 10 12 13 14 15]) = [0.3 0.1 0.01 0.3 0.2 0.15 0.1 0.2 0.4 0.01 0.15 0.4];

beta = C .*R0.*gam;     %taux de transmibilité

% Integration parameters 
%IC = [(N(1)-0.25e1) 0.25e1 0 ; (N(2)-0.33e1) 0.33e1 0 ; (N(3)-0.42e1) 0.42e1 0 ; (N(4)-0.62e1) 0.62e1 0]; % seed a few infected individuals
IC = [(N(1)-N(1)*0.00001) N(1)*0.00001 0 ; (N(2)-N(2)*0.00001) N(2)*0.00001 0 ; (N(3)-N(3)*0.00001) N(3)*0.00001 0  ; (N(4)-N(4)*0.00001) N(4)*0.00001 0  ];    %infecte au debut
%IC = [(N(1)-N(1)*0.2) N(1)*0.2 0 ; (N(2)-N(2)*0.15) N(2)*0.15 0 ; (N(3)-N(3)*0.25) N(3)*0.25 0  ; (N(4)-N(4)*0.2) N(4)*0.2 0  ];    %infecté au debut du confinement
%IC = [(N(1)-N(1)*0.000001) N(1)*0.000001 0 ; (N(2)-N(2)*0.1) N(2)*0.1 0 ; (N(3)-N(3)*0.15) N(3)*0.15 0  ; (N(4)-N(4)*0.1) N(4)*0.1 0  ];    %infecté a la fin du confinement

%tspan = [0 365]; % in days 
tspan = [0 100];

% simulations
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
sol = ode45(@sir,tspan,IC,options);
soly = [sum(sol.y(1:4,:)) ;
        sum(sol.y(5:8,:)) ;
        sum(sol.y(9:12,:));
        sum(sol.y)];

figure(2); clf;
semilogy(sol.x, soly(1:4,:));
legend('S','I','R','Total');

disp(['Percentage of contaminated 100*(I+R)/N at day ' num2str(tspan(2)), ': ', ...
  num2str(100*sum(sol.y(5:12,end)/N))]);


    function dxdt = sir(~,xx)  % nested function 
        % SIR is the ODE rhs 
        S = xx(1:4);
        I = xx(5:8); 
        R = xx(9:12);
        
        
        dxdt = [ (- beta * I .* S./N);
                 (beta * I .* S./N - gam .* I);
                 (gam .* I) ]; 
             
    end                         % end nested function sir

end                             % end main function run_SIR 