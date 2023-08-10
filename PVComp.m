function FF = PVComp(G, Ta, Rp, Fig)
%{
    NAME: ABDULLAH ESSA
    SID: 201256467
    ELEC5564M - Electric Power Generation by Renewable Sources
    Assignment (1), PV Simulation Script - Variable (G, Ta and Rp)
    Last Updated: 25 July 2023
%}

%% Inputs
% G -> irradiation
% Ta -> ambient temperature
% Rp -> PV parallel resistance
% Fig -> plotting figure number (integer)

%% Outputs
% FF -> Array of Fill Factor for each scenario iteration

%% Initialization
% Define the parameters for the PV cell
A = 1.72;  % diode ideality factor
k = 1.380658e-23;  % boltzmann constant.
Ior = 19.9693e-6;  % reverse saturation current.
ns = 40;  % number of series cells per module.
np = 2;  % number of parallel cells per module.
Tr = 301.18;  % reference temperature (stc).
q = 1.6e-19;  % electron charge.
Eg = 1.1;  % silicon bandgap energy.
Iscr = 3.3;  % short circuit current (reference, stc) at Tr.
Ki = 1.7e-3;  % Temperature coefficient of I_sc.
Rs = 5e-5;  % series resistance.
%Rp = 5e+5;  % parallel resistance.

% initialize plot color array
C = {'b','r',"#77AC30","#4DBEEE",'m'};

% initialize the array for Tc
Tc =  zeros(length(Ta),1);

% initialize fill factor array to contain elements equal to number of
% scenarios
size = length(Ta)*length(G);
FF = zeros(size,1);  
FF_Counter = 1;  % counter for fill factor array increment

% if Rp is being tested then make G constant otherwise Rp is constant
if(length(Rp)>1)
    length_ = length(Rp);
    G = G(1);  % protect against if G was an array with array of Rp
    G = ones(length(Rp),1).*G;
else
    length_ = length(G);
    Rp = Rp(1);  % protect against if G was an array with array of Rp
    Rp = ones(length(G),1).*Rp;
end
%% Loop
for i = 1:length_
    for ii = 1:length(Ta)
        % cell surface temperature
        Tc(ii) = Ta(ii) + 0.2*G(i) + 273.18;   

        % Diode leakage current
        Is = Ior*(Tc(ii)/Tr)^3*exp(q*Eg*((1/Tr-1/Tc(ii))/(k*A)));

        % Short circuit current
        Isc = (Iscr +(Ki*(Tc(ii)-Tr)))*G(i);

        % Defining the loop boundaries -> initial guesses
        % short circuit condition is the initial:
        I = Isc;  % current is equal to Isc at short circuit condition
        V = 0;  % voltage is zero at short circuit condition

        % Newton Raphson criterion initialization
        tol = 1e-3;  % initializing the set maximum error the root should be.
        error = 1;  % initializing the error to be used in loop.

        % Implementation of Newton Raphson to find the Voc 
        while I > 0 
            while error >= tol
                % equation of I equal to 0; f(I)=0
                funcI_1 = I - np*Isc + np*Is*(exp((q*(V/ns + I*Rs/np))/(A*k*Tc(ii)) ) - 1) + (V*np/ns + I*Rs)/Rp(i);

                % derivative of the equation of I equal to 0; f'(I)
                funcI_2 = 1 - (q*Is*Rs)/(A*k*Tc(ii))*exp((q*(V/ns + I*Rs/np))/(A*k*Tc(ii))) + (Rs/Rp(i));

                % compute error
                error = funcI_1/funcI_2;

                % compute new Y-axis value of (I)
                I = I - error;
            end
                    % error reset for the new data point
        error = 1;

        % increment voltage by the tolerance for accuracy
        V = V + tol;
    end

    % store the root as the open circuit voltage
    Voc = V;

    % plotting the current, voltage and power characteristics, 
    % labelling key points.

    % initializing the voltage array according to boundary conditions
    V = 0:tol:Voc;

    % initializing the current array 
    I = zeros(length(V),1);
    I(1,1) = Isc;  % initializing the short circuit condition

    % initializing the power output array
    P = zeros(length(V),1);

    % Computing the values of (I) & (P)
    for z = 1:1:length(V)   
        if(z < length(V))
            I(z) = np*Isc - np*Is*(exp((q*(V(z)/ns + I(z)*Rs/np))/(A*k*Tc(ii)) ) - 1) - (V(z)*np/ns + I(z)*Rs)/Rp(i);
        end

        % Protect the model from giving negative output current 
        % open circuit condition
        if ((I(z) < 0) || (z == length(V)))
            I(z) = 0;
        end

        % output power compute
        P(z) = V(z).*I(z);
    end

    %% Maximum Power Point Tracking (MPPT)& Fill Factor Computation
    P_MPP = max(P);
    for n = 1:1:length(V)
        if(P_MPP == P(n))
            V_MPP = V(n);
            I_MPP = I(n);
        end
    end

    FF(FF_Counter) = (I_MPP*V_MPP)/(Isc);
    FF_Counter = FF_Counter + 1;

    %% PLOTTING
    if(length(Rp)>1)
        % IV Characteristics
        figure (Fig);
        yyaxis left;
        % add a legend on every irradiance or Temperature iteration 
        if(length(Ta)>length(G))
            plot(V,I,'-','color',C{ii},'LineWidth',2,'DisplayName',['IV Characteristics (Rp = ',num2str(Rp(i)),')']);
        else
            plot(V,I,'-','color',C{i},'LineWidth',2,'DisplayName',['IV Characteristics (Rp = ',num2str(Rp(i)),')']);
        end
        hold on;
        scatter(V_MPP, I_MPP, 'r','DisplayName',['MPP (Rp = ',num2str(Rp(i)),')']);
        xlabel('Voltage (V)');
        ylabel('Current (A)');
        grid on;
        legend('-DynamicLegend','Location','east');
        % show the MPP coordinate
        t = text(V_MPP+0.05*V_MPP, I_MPP, ['(',num2str(V_MPP),',',num2str(I_MPP),')']);
        t.FontSize = 7;
        title(['PV & IV Characteristics - for Rp']);
    
        % PV Characteristics
        figure(Fig);
        yyaxis right;
        % add a legend on every irradiance or Temperature iteration 
        if(length(Ta)>length(G))
            plot(V,P,'--','color',C{ii},'LineWidth',2,'DisplayName', ['PV Characteristics (Rp = ',num2str(Rp(i)),')']);
        else
            plot(V,P,'--','color',C{i},'LineWidth',2,'DisplayName', ['PV Characteristics (Rp = ',num2str(Rp(i)),')']);
        end
        hold on;
        scatter(V_MPP, P_MPP, 'r','DisplayName',['MPP (Rp = ',num2str(Rp(i)),')']);
        xlabel('Voltage (V)');
        ylabel('Power (W)');
        grid on;
        legend('-DynamicLegend','Location','northwest');
        % show the MPP coordinate
        t = text(V_MPP+0.05*V_MPP, P_MPP, ['(',num2str(V_MPP),',',num2str(P_MPP),')']);
        t.FontSize = 7;
    else
        % IV Characteristics
        figure (Fig);
        yyaxis left;
        % add a legend on every irradiance or Temperature iteration 
        if(length(Ta)>length(G))
            plot(V,I,'-','color',C{ii},'LineWidth',2,'DisplayName',['IV Characteristics (G = ',num2str(G(i)),')&(Ta = ',num2str(Ta(ii)),')']);
        else
            plot(V,I,'-','color',C{i},'LineWidth',2,'DisplayName',['IV Characteristics (G = ',num2str(G(i)),')&(Ta = ',num2str(Ta(ii)),')']);
        end
        hold on;
        scatter(V_MPP, I_MPP, 'r','DisplayName',['MPP (G = ',num2str(G(i)),')&(Ta = ',num2str(Ta(ii)),')']);
        xlabel('Voltage (V)');
        ylabel('Current (A)');
        grid on;
        legend('-DynamicLegend','Location','east');
        % show the MPP coordinate
        t = text(V_MPP+0.05*V_MPP, I_MPP, ['(',num2str(V_MPP),',',num2str(I_MPP),')']);
        t.FontSize = 7;
        title(['PV & IV Characteristics - for Ta(',num2str(Ta),') & G =(', num2str(G),')']);
    
        % PV Characteristics
        figure(Fig);
        yyaxis right;
        % add a legend on every irradiance or Temperature iteration 
        if(length(Ta)>length(G))
            plot(V,P,'--','color',C{ii},'LineWidth',2,'DisplayName', ['PV Characteristics (G = ',num2str(G(i)),')&(Ta = ',num2str(Ta(ii)),')']);
        else
            plot(V,P,'--','color',C{i},'LineWidth',2,'DisplayName', ['PV Characteristics (G = ',num2str(G(i)),')&(Ta = ',num2str(Ta(ii)),')']);
        end
        hold on;
        scatter(V_MPP, P_MPP, 'r','DisplayName',['MPP (G = ',num2str(G(i)),')&(Ta = ',num2str(Ta(ii)),')']);
        xlabel('Voltage (V)');
        ylabel('Power (W)');
        grid on;
        legend('-DynamicLegend','Location','northwest');
        % show the MPP coordinate
        t = text(V_MPP+0.05*V_MPP, P_MPP, ['(',num2str(V_MPP),',',num2str(P_MPP),')']);
        t.FontSize = 7;
    end
end
end

    




