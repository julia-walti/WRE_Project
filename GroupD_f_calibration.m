function [Ksat_val, c_val, tsub_val, z_val, NS_val] = GroupD_f_calibration(Q_obs, phi, temperature, n_years, P, sw, s1, kc, n, Qb, tsup, A, dt, tsub1, c1, Ksat1,z1, doTest)
    %iteration until convergence NS>0.87 - can it be as a conditions in a loop?
    N_it=10000;                           %[-]    Number of iterations for the simulated annealing algorithm 
    %should be 10 000 but because of long running time i put it to 100
    cr=1/1200;                          %cooling rate
    Nb_NS_stored=1;                     %[-], accepted values of NS

    % Preallocation of vectors
    Ksat=zeros(N_it, 1);                 %[m/s]
    c=zeros(N_it,1);                    %[-]
    tsub=zeros(N_it,1);                 %[h]
    z=zeros(N_it,1);                    %[mm]
    NS=zeros(N_it,1);                   %[-]
    T_SA=zeros(N_it,1);                 %[-]

    % Initial values are the arbitrary avalues defined earlier 
    Ksat(1)=Ksat1;                       %[m/s]  Hydraulic conductivity for saturated soil
    c(1)=c1;                             %[-]    Exponent for power-law relation L(s)
    tsub(1)=tsub1;                       %[h]    Mean sub-superficial residence time
    z(1)=z1;                             %[mm]   Root zone depth

    % First NS index
    [Q, R, I, s, L, ET] = GroupD_f_hydromodel(phi, temperature, n_years, P, Ksat1, sw, s1, kc, n, Qb, tsup, tsub1, A, c1, dt, z1, doTest);
    Q=Q';
    NS_function=@(Q) 1-sum((Q_obs-Q).^2)/sum((Q_obs-mean(Q_obs)).^2); %[-]
    NS(1)=NS_function(Q);
    NS_old = NS(1);

    % Std
    sigma_Ksat=0.05*(10^(-5)-10^(-7));  %[m/s]
    sigma_c=0.05*(20-1);                %[-]
    sigma_tsub=0.05*(400-1);            %[h]
    sigma_z=0.05*(2000-1);              %[mm]

    for i=1:N_it
           T_SA(i) = exp(-cr*i);            %[-]
           Ksat_new=TruncNormRnd(Ksat(Nb_NS_stored),sigma_Ksat,10^(-7),10^(-5));   %[m/s]
           c_new=TruncNormRnd(c(Nb_NS_stored),sigma_c,1,20);                       %[-]
           tsub_new=TruncNormRnd(tsub(Nb_NS_stored),sigma_tsub,1,400);             %[h]
           z_new=TruncNormRnd(z(Nb_NS_stored),sigma_z,1,2000);                     %[mm]

           %run model
           [Q_new, R, I, soil_saturation, L, ET]= GroupD_f_hydromodel(phi,temperature,n_years,P,Ksat_new,sw,s1,kc,n,Qb,tsup,tsub_new,A,c_new,dt,z_new,doTest);

           %NS new
           %Q_new=Q_new';
           NS_new=NS_function(Q_new'); %[-]

           if NS_new>NS_old
                %Update index counting accepted values in preallocated vector
                Nb_NS_stored=Nb_NS_stored+1;
                %Save the current parameter set
                Ksat(Nb_NS_stored)=Ksat_new;
                c(Nb_NS_stored)=c_new;
                tsub(Nb_NS_stored)=tsub_new;
                z(Nb_NS_stored)=z_new;
                NS(Nb_NS_stored)=NS_new;
                NS_old=NS_new;
           else
               if rand<exp((NS_new-NS_old)/T_SA(i))
               % Save based on probability annealing condition
      %not sure how to structure the probability condition
                %Update index counting accepted values in preallocated vector
                Nb_NS_stored=Nb_NS_stored+1;
                %Save the current parameter set
                Ksat(Nb_NS_stored)=Ksat_new;
                c(Nb_NS_stored)=c_new;
                tsub(Nb_NS_stored)=tsub_new;
                z(Nb_NS_stored)=z_new;
                NS(Nb_NS_stored)=NS_new;
                NS_old=NS_new;
               end           
           end
    end

    %clean empty vectors
    Ksat(Nb_NS_stored+1:end)=[];
    c(Nb_NS_stored+1:end)=[];
    z(Nb_NS_stored+1:end)=[];
    tsub(Nb_NS_stored+1:end)=[];
    NS(Nb_NS_stored+1:end)=[];
    T_SA(Nb_NS_stored+1:end)=[];
    
    % Save the values of the parameters during the calibration
    writematrix([Ksat, c, z, tsub, NS, T_SA],'GroupD_calibration_res.xls') 
    
    % Return values 
    [NS_val, idx] = max(NS) ; 
    Ksat_val = Ksat(idx);
    c_val = c(idx);
    z_val = z(idx);
    tsub_val = tsub(idx);

end