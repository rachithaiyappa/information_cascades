%personally setting inital conditions

%N unique links instead of N/2
%Marc's initial conditions of having non zero dosages for initially
%infected individuals

%corrected the bug of overwriting assigned doses to 0

clear all
tic
count = 1;
fileID = fopen('phi_star_data_corrections.txt','w');
for p = 0 : 0.1 : 1
    count = 1;
    phi_star = 999*ones(11,2);
for start = 0 : 0.1 : 1
    clearvars -except p start count fileID phi_star
    p 
    start 
    tic
    N = 100;     %number of individuals
    T = 1000;    %Total number of time steps
    mem_time_steps = 12; %indivdual remembers previous mem_time_steps

    %removed : r = 3 %to avoid confusion with matlab automatically alloting 0 to unknown values
    %susceptible : s = 1
    %infected i : = 2

%     p = 1; %probability of doze transfer
    rho = 1; % probability of becoming susceptible after recovery
    r = 1; %probability of recovery when dose level drops below threshold

    t = 1:1:T+1;      %obtaining the time steps
    focal(:,1) = (1:1:N)';   %storing all the individuals
    partner(:,1) = (1:1:N)';   %storing all the partners
    indi_state(:,1) = ones(N,1);

    %doses
    d_ini(:,1) = zeros(N,1);
    d(:,1) = zeros(N,1); %doses
    D(:,1) = zeros(N,1); %cumulative dose distribution

    d_star(:,1) = 3*ones(N,1); %threshold of each individual
%     split = floor(N/2); %dividing into three groups for i,s,r
    infec(:,1) = zeros(int32(start*N),1);
    % sus(:,1) = zeros(split,1);

    %randomly grouping into i, s or r
    infec(:,1) = randperm(N,int32(start*N))'; %infected
    indi_state(infec(:,1),1) = 2;   %infected
    % sus(:,1) = randi([1 N],split,1); %susceptible
    % indi_state(sus(:,1),1) = 1;   %susceptible
    
    I = find(indi_state(:,1)==2);   %indices of the initially infected individuals
    if size(I,1) ~= 0
        for i = 1:size(I,1)
            d_ini(I(i),1:1:mem_time_steps) = p; %settting non zero initial dosages for the infected individuals
        end
    else
        d_ini(:,1:1:mem_time_steps) = 0;
    end

    phi(1) = numel(find(indi_state(:,1)==2))/N;
    for k = 2:size(t,2)   %time loop
        if mod(k,500) == 0
            k
        end
    %     'time loop'
        for n = 1:N    %individual loop
    %         'individual loop 1'
            focal(n,k) = n; %focal individual
            partner(n,k) = randi(N); %the partner of the focal individual
            z1(focal(n,k),k) = rand; %random number to decide whether dosage is postive or zero
            %if ONE is infected and the OTHER is susceptible
            if (indi_state(focal(n,k),k-1) == 1 && indi_state(partner(n,k),k-1) == 2)
    %             'blah3'
                if p >= z1(focal(n,k),k)
    %                 k
    %                 'positive dose'
                    d(n,k) = lognrnd(1,0.658);     %dosage is log noramlly distributed
%                     d(focal(n,k),k) = 1;
%                     d(partner(n,k),k) = 0; %the other does not receive any dose
                else
    %                 'interaction but no dose'
                    d(focal(n,k),k) = 0;     %interaction but no dose
%                     d(partner(n,k),k) = 0;   %the other does not receive any dose
                end
            %elseif both are infected
            elseif (indi_state(focal(n,k),k-1) == 2 && indi_state(partner(n,k),k-1) == 2)
                if p >= z1(focal(n,k),k)
    %                 'positive dose'
                    d(n,k) = lognrnd(1,0.658);     %dosage is log noramlly distributed
%                     d(focal(n,k),k) = 1;
%                     d(partner(n,k),k) = 0; %STILL RECEIVES NO DOSE
                else
    %                 'interaction but no dose'
                    d(focal(n,k),k) = 0;     %interaction but no dose
%                     d(partner(n,k),k) = 0;   %the other does not receive any dose
                end   
            %else there is no interaction between the indivudal and the partner
            
            else
    %             'no interaction at all'
                d(focal(n,k),k) = 0; %no interaction at all
%                 d(partner(n,k),k) = 0;   %the other does not receive any dose
            end
        end


        %updating cumulative dose
        D(:,k) = 0;
    %     k
        m = k-(mem_time_steps);
        if m > 0
    %       'cumu dose later stages loop'
            D(:,k) = sum(d(:,k:-1:m),2);
        end
%         while m <= 0
%     %         'cumu dose initial stages loop'
%             m = m + 1;
%         end
%     %     'cumu dose calc'
%         D(:,k) = sum(d(:,k:-1:m),2);
        if m <= 0
    %         'cumu dose initial stages loop'
             while m <= 0
                m = m + 1;
             end
        mem_time_steps-(k-m);
            D(:,k) = sum(d(:,k:-1:m),2) + sum(d_ini(:,1:1:(mem_time_steps-(k-m))),2);
        end
        
        
        for n = 1:N
%             index = find(indi(:,k) == n);
%             index_prev = find(indi(:,k-1) == n);
    %         'updating indi_state loop'
            if (D(n,k) >= d_star(n,1)) && (indi_state(n,k-1) == 1)  %if threshold is met and is susceptible
    %             'above thresh becomes infected'
                indi_state(n,k) = 2; %susceptible becomes infected
            elseif (D(n,k) < d_star(n,1)) && ((indi_state(n,k-1)) == 2) %recovery : if drops bel threshold and is infected
                z2 = rand;
                if r >= z2 
    %                 'below thresh'
                    z3 = rand;
                    if rho >= z3
    %                     'below thresh becomes susceptible'
                        indi_state(n,k) = 1; %infected individual recovers but becomes susceptible
                    else
    %                     'below thresh becomes immune '
                        indi_state(n,k) = 3; %recovered becomes immune and thus removed
                    end
                else
    %                  'below thresh stays infected'
                    indi_state(n,k) = 2; %individual remains infected even though dose leve dropped below threshold
                end
            else
%                 'retains previous state'
                indi_state(n,k) = indi_state(n,k-1); %retains the previous state
            end
        end
%         phi(k) = numel(find(indi_state(:,k)==2))/N; %dynamics of phi for a particular value of p
    end

%     calculating avaerage phi_star from the last 100 time steps
    summ = 0;
    number = 0;
    for k = size(t,2):-1:size(t,2)-99
        if k == size(t,2) || k == size(t,2)-99
            'calculating average'
        end
        number = numel(find(indi_state(:,k)==2))/N;   %finding the number of infected in the given time step
        summ = summ + number; 
    end
    phi_star(count,1) = p
    phi_star(count,2) = summ/100
    
    toc
    
    fprintf(fileID,'%d\t%d\t%d\t%d\n',p,start,phi_star(count,2),toc);
    count = count + 1;
end
end
toc
