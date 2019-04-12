%personally setting inital conditions
clear all
tic
count = 1;
fileID = fopen('phi_star_data.txt','w');
for p = 0.2 : 0.01 : 0.3
    count = 1;
    phi_star = 999*ones(10,2);
for start = 0.1 : 0.1 : 1
    clearvars -except p start count fileID phi_star
%     p = 0.3
%     start = 1
    tic
    N = 1e05;     %number of individuals
    T = 1e03;    %Total number of time steps
    mem_time_steps = 12; %indivdual remembers previous mem_time_steps

    %removed : r = 3 %to avoid confusion with matlab automatically alloting 0 to unknown values
    %susceptible : s = 1
    %infected i : = 2

%     p = 1; %probability of doze transfer
    rho = 1; % probability of becoming susceptible after recovery
    r = 1; %probability of recovery when dose level drops below threshold

    t = 1:1:T+1;      %obtaining the time steps
    indi(:,1) = (1:1:N)';   %storing all the individuals
    indi_state(:,1) = ones(N,1);

    %doses
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

    phi(count,1) = numel(find(indi_state(:,1)==2))/N;
    for k = 2:size(t,2)   %time loop
        if mod(k,500) == 0
            k
        end
    %     'time loop'
        indi(:,k) = randperm(N)'; %assigning partners to each indvidual
        for n = 1:N/2     %individual loop
    %         'individual loop 1'
    
            z1(indi(n,k),k) = rand; %random number to decide whether dosage is postive or zero
            z1(indi((N/2)+n,k),k) =  z1(indi(n,k),k);
            %if ONE is infected and the OTHER is susceptible
            if (indi_state(indi(n,k),k-1) == 1 && indi_state(indi((N/2)+n,k),k-1) == 2)
    %             'blah3'
                if p >= z1(indi(n,k),k)
    %                 k
    %                 'positive dose'
                    d(indi(n,k),k) = lognrnd(1,0.658);     %dosage is log noramlly distributed
%                     d(indi(n,k),k) = 1;
                    d(indi((N/2)+n,k),k) = 0; %the other does not receive any dose
                else
    %                 'interaction but no dose'
                    d(indi(n,k),k) = 0;     %interaction but no dose
                    d(indi((N/2)+n,k),k) = 0;   %the other does not receive any dose
                end
            elseif (indi_state(indi(n,k),k-1) == 2 && indi_state(indi((N/2)+n,k),k-1) == 1)
                if p >= z1(indi((N/2)+n,k),k)
    %                 k
    %                 'positive dose'
                    d(indi((N/2)+n),k) = lognrnd(1,0.658);     %constant doses
%                     dindi((N/2)+n),k) = 1;
                    d(indi(n,k),k) = 0;   %the one does not receive any dose  
                else
    %                 'interaction but no dose'
                      d(indi((N/2)+n),k) = 0;     %interaction but no dose
                      d(indi(n,k),k) = 0;   %the one does not receive any dose  
                end
            %elseif both are infected
            elseif (indi_state(indi(n,k),k-1) == 2 && indi_state(indi((N/2)+n,k),k-1) == 2)
                if p >= z1(indi(n,k),k)
    %                 'positive dose'
                    d(indi(n,k),k) = lognrnd(1,0.658);     %dosage is log noramlly distributed
%                     d(indi(n,k),k) = 1;
                    d(indi((N/2)+n,k),k) = lognrnd(1,0.658); %dosage is log noramlly distributed
%                     d(indi((N/2)+n,k),k) = 1;
                else
    %                 'interaction but no dose'
                    d(indi(n,k),k) = 0;     %interaction but no dose
                    d(indi((N/2)+n,k),k) = 0;   %the other does not receive any dose
                end   
            %else there is no interaction between the indivudal and the partner
            else
    %             'no interaction at all'
                d(indi(n,k),k) = 0; %no interaction at all
                d(indi((N/2)+n,k),k) = 0;   %the other does not receive any dose
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
        while m <= 0
    %         'cumu dose initial stages loop'
            m = m + 1;
        end
    %     'cumu dose calc'
        D(:,k) = sum(d(:,k:-1:m),2);
        

        for n = 1:N
%             index = find(indi(:,k) == n);
%             index_prev = find(indi(:,k-1) == n);
    %         'updating indi_state loop'
            if (D(n,k) >= d_star(n,1)) && (indi_state(n,k-1) == 1)  %if threshold is met and is susceptible
    %             'above thresh becomes infected'
                indi_state(n,k) = 2; %susceptible becomes infected
            elseif (D(n,k) < d_star(n,1)) && ((indi_state(n,k-1)) == 2) %recovery : if drops bel threshold and is infected
                z2(n,k) = rand;
                if r >= z2(n,k) 
    %                 'below thresh'
                    z3(n,k) = rand;
                    if rho >= z3(n,k)
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
%         phi(count,k) = numel(find(indi_state(:,k)==2))/N; %dynamics of phi for a particular value of p
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
