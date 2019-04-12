clear all
phi_star = 999*ones(11,2);
count = 1;
for p = 0:0.1:1
    p
    tic
    N = 1e05;     %number of individuals
    T = 1e02;    %Total number of time steps
    mem_time_steps = 12; %indivdual remembers previous mem_time_steps

    %removed : r = 0
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
    split = floor(N/2); %dividing into three groups for i,s,r
    infec(:,1) = zeros(split,1);
    % sus(:,1) = zeros(split,1);

    %randomly grouping into i, s or r
    infec(:,1) = randperm(N,split)'; %infected
    indi_state(infec(:,1),1) = 2;   %infected
    % sus(:,1) = randi([1 N],split,1); %susceptible
    % indi_state(sus(:,1),1) = 1;   %susceptible
    phi(count,1) = numel(find(indi_state(:,1)==2))/N;

    for k = 2:size(t,2)   %time loop
        if mod(k,1e03) == 0
            k
        end
    %     'time loop'
        indi(:,k) = randperm(N)'; %assigning partners to each indvidual
        for n = 1:N/2     %individual loop
            if mod(n,1e05) == 0
                n
            end
    %         'individual loop 1'

            z1 = rand; %random number to decide whether dosage is postive or zero
            %if ONE is infected and the OTHER is susceptible
            if (indi_state(indi(n,k),k-1) == 1 && indi_state(indi((N/2)+n,k),k-1) == 2)
    %             'blah3'
                if p >= z1
    %                 k
    %                 'positive dose'
                    d(n,k) = lognrnd(1,0.774);     %dosage is log noramlly distributed
                else
    %                 'interaction but no dose'
                    d(n,k) = 0;     %interaction but no dose
                end
            %else there is no interaction between the indivudal and the partner
            else
    %             'no interaction at all'
                d(n,k) = 0; %no interaction at all
            end
            %if OTHER is infected and the ONE is susceptible
            if (indi_state(indi(n,k),k-1) == 2 && indi_state(indi((N/2)+n,k),k-1) == 1)
    %             'blah3'
                if p >= z1
    %                 k
    %                 'positive dose'
                    d(indi((N/2)+n),k) = 1;     %constant doses
                else
    %                 'interaction but no dose'
                      d(indi((N/2)+n),k) = 0;     %interaction but no dose
                end
    %             else there is no interaction between the indivudal and the partner
            else
    %             'no interaction at all'
                d(indi((N/2)+n),k) = 0; %no interaction at all
            end
        end


        %updating cumulative dose
        D(:,k) = 0;
    %     k
        m = k-(mem_time_steps-1);
        if m > 0
    %         'cumu dose later stages loop'
            D(:,k) = sum(d(:,k:-1:m),2);
        end
        while m <= 0
    %         'cumu dose initial stages loop'
            m = m + 1;
        end
    %     'cumu dose calc'
        D(:,k) = sum(d(:,k:-1:m),2);


        for n = 1:N
    %         'updating indi_state loop'
            if (D(indi(n,k),k) >= d_star(indi(n,k),1)) && (indi_state(indi(n,k),k-1) == 1)  %if threshold is met and is susceptible
    %             'above thresh becomes infected'
                indi_state(indi(n,k),k) = 2; %susceptible becomes infected
            elseif (D(indi(n,k),k) < d_star(indi(n,k),1)) && ((indi_state(indi(n,k),k-1)) == 2) %recovery : if drops bel threshold and is infected
                z2 = rand;
                if r >= z2 
    %                 'below thresh'
                    z3 = rand;
                    if rho >= z3
    %                     'below thresh becomes susceptible'
                        indi_state(indi(n,k),k) = 1; %infected individual recovers but becomes susceptible
                    else
    %                     'below thresh becomes immune '
                        indi_state(indi(n,k),k) = 0; %recovered becomes immune and thus removed
                    end
                else
    %                  'below thresh stays infected'
                    indi_state(indi(n,k),k) = 2; %individual remains infected even though dose leve dropped below threshold
                end
            else
    %             'retains previous state'
                indi_state(indi(n,k),k) = indi_state(indi(n,k),k-1); %retains the previous state
            end
        end
        phi(count,k) = numel(find(indi_state(:,k)==2))/N; %dynamics of phi for a particular value of p
    end

%     %calculating avaerage phi_star from the last 100 time steps
%     summ = 0;
%     number = 0;
%     for k = size(t,2):-1:size(t,2)-100
%         if k == size(t,2) || k == size(t,2)-100
%             'calculating average'
%         end
%         number = numel(find(indi_state(:,k)==2))/N;   %finding the number of infected in the given time step
%         summ = summ + number; 
%     end
%     phi_star(count,1) = p
%     phi_star(count,2) = summ/100
    toc
    count = count + 1;
end
