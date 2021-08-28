function Output = CoarseGraining(E,par)
%=================================================================================
%       The dynamics-preserving clustering method to analyze single-cell data      
%       Input: 
%       E: the symmetric weight matrix construted from data
%       par.
%       @Date 2018.4.12, original version from SAVI package
%=================================================================================

% parameters of algorithms:

% Informations of networks:

% =================================================================================

D=diag(sum(E,2)); % The diagonal matrix
P=D\E; % D^{-1} *E, the stochastic matrix

mu=diag(D); % size : n * 1
mu=mu./sum(mu); % the equilibrium distribution

num=length(mu); % The size of the matrix.

if par.reduce_large_scale
    [~,ydata] = pca(par.data',20);
    k = par.reduce_num_meta_cell;
    stream = RandStream('mlfg6331_64');  % Random number stream
    options = statset('UseParallel',1,'UseSubstreams',1,'Streams',stream);
    classes = kmeans(ydata.M,k,'Replicates',150,'Options',options);
    classes = classes';
    P = Comput_Phat( P, mu, classes, k );
    mu_hat = zeros(k,1);
    data_orig = par.data;
    N_gene = size(data_orig,2);
    data_reduced = zeros(k,N_gene);
    
    for i =1:1:k
        mu_hat(i) = sum(mu(classes == i));
        data_reduced(i,:) = mean(data_orig(classes==i,:));
    end
    mu = mu_hat;
    Output.data = data_reduced;
    Output.reduce_class = classes;
    data = data_reduced;
else
    data = par.data;
end


network_out.mu = mu;
network_out.P = P;
%===============================================================

%   The Body of the Algorithm

%===============================================================
Trials = par.trials;   % The trials of algorithm
E_best = inf; 
% Start :

if strcmp(par.initial,'pca')
        [~,ydata] = pca(data',10);
end

for trials = 1: Trials
 % Initialization :
    
    if strcmp(par.options,'fixed')
        k = par.K_cluster;
    else
        % the initial guess of cluster numbers
        K_max = par.K_max;
        K_min = par.K_min;

        % settings in the simulated anealing 
        T_max = par.T_max;
        T_min = par.T_min; 

        Rloop = par.Rloop;  % The iterative times at each temperature
        alpha = par.alpha;  % The cooling rate
        k = fix( rand()*(K_max-K_min)+K_min );
        T = T_max;
    end
    
    if strcmp(par.initial,'pca')
        classes = kmeans(ydata.M,k,'Replicates',150);
        classes = classes';
    elseif strcmp(par.initial,'tsne')
        Dist = squareform(pdist (data,par.choice_distance));
        ydata = tsne_d(Dist);
        classes = kmeans(ydata,k,'Replicates',150);
        classes = classes';
    elseif strcmp(par.initial,'other')
        ydata = par.init_score;
        classes = kmeans(ydata,k,'Replicates',150);
        classes = classes';
    end
    

    if strcmp(par.initial,'random')
        randvec = rand(1,num);
        classes = floor(k*randvec)+1;
        classes( randsample(num,k) ) = 1:k; % Force the k clusters are not empty.
    end
    

    % Comput P_hat :
    P_hat = Comput_Phat( P, mu, classes, k );

    % Comput J :
    J_current = Objective_J( P, mu, P_hat, classes   );

    %     P_hat_old = P_hat;
    %     classes_old = classes;
    %     J_old = J_current;

    while 1
        % Compute P_hat_new :
        P_hat_new = Comput_Phat( P, mu, classes, k );

        % Compute the classes_new :
        [E_bar classes_new k] = Distance_NodeClass(  P, mu, P_hat_new, classes );

        % Compute P_hat_new again:
        P_hat_new = Comput_Phat( P, mu, classes_new, k );

        % Compute J_new :
        J_new = Objective_J( P, mu, P_hat_new, classes_new )

        if J_new < J_current
            P_hat = P_hat_new;
            classes = classes_new;
            J_current = J_new;
        else
            break;
        end
        
    end

    % Compute energy E :
    E_current = Validity( P, mu, P_hat, classes );

    % The initialized as the old :
    E_old = E_current;

    % update the optimal state :
    
    if (E_current <E_best)
    P_hat_best =  P_hat;
    classes_best = classes;
    E_best = E_current;
    end
    
    
    %------------------------------------------------
    %     Main Loop of Simulated Annealing       %
    %------------------------------------------------
%%
if strcmp(par.options,'optimal')   
    iterNum = 0; % times of cooling
    flag=1;

    while 1
        T = T
        iterNum = iterNum+1;

        % Start of subloop at given temperature :
        for R = 1:Rloop
            choice = fix( rand()*2+1);
            % Two functions of our proposal :
            switch choice
                case 1
                    [classes, k] = Split_class( classes, P_hat , mu );
                case 2
                    [classes, k] = Delete_class( classes, P_hat );
                                  
                otherwise
                    'Default';
                    flag=0;
                    break;
            end

            % Comput P_hat :
            P_hat = Comput_Phat( P, mu, classes, k );

            % Comput J :
            J_current = Objective_J( P, mu, P_hat, classes );

            %--------------------------------------------------------------
            while 1
                % Compute P_hat_new :
                P_hat_new = Comput_Phat( P, mu, classes, k );

                % Compute the classes_new :
                [E_bar, classes_new k] = Distance_NodeClass( P, mu,  P_hat_new, classes );

                % Compute P_hat_new again:
                P_hat_new = Comput_Phat( P, mu, classes_new, k );

                % Compute J_new :
                J_new = Objective_J( P, mu,  P_hat_new, classes_new );
                if J_new < J_current
                    P_hat = P_hat_new;
                    classes = classes_new;
                    J_current = J_new;
                else
                    break;
                end

            end

            %--------------------------------------------------------------
            % Compute the energy E:
            E_current = Validity( P, mu, P_hat, classes );
            k = k

            %--------------------------------------------------------------
            % Accept or reject :
            r = rand();
            accept_P = 1/(1+exp(-(E_current-E_old)/T));
            %accept_P =  exp(-(E_current-E_old)/T);

            if E_current < E_old | r<accept_P
                P_hat_old = P_hat;
                classes_old = classes;
                E_old = E_current;
            end

            %--------------------------------------------------------------
            % If we a new E is more optimal, then replace the existed
            % optimal state :
            if E_current < E_best
                P_hat_best = P_hat_old;
                classes_best = classes_old;
                E_best = E_old;
            end

        end
        % End of subloop at some temperature

        %--------------------------------------------------------------
        % Cooling :
        T = alpha * T

        if T < T_min | flag==0
            break;
        end

        %--------------------------------------------------------------
        % Set the optimal state at the current temperature as the initial
        % state at the next temperature :
        P_hat = P_hat_best;
        classes = classes_best;
        E_old = E_best;

    end

    % choose the optinal one among the trials :
    %{
    if trials==1 || E_best < E_final
        E_final = E_best;
        P_hat_final = P_hat_best;
        classes_final = classes_best;
    end
    %}
end

    % choose the optinal one among the trials :
    %
    if trials==1 || E_best < E_final
        E_final = E_best;
        P_hat_final = P_hat_best;
        classes_final = classes_best;
    end
    %



end
E_best
%}
% ===============================
%% Make orders :

class = classes_best;
k = size(P_hat_best,1);

for i = 1:k
    Index = find(class==i);
    first(i) = Index(1);
    clear Index;
end

firstIndex = sort(first);

perm = zeros(k,1);
for i = 1:k
    j = find(first==firstIndex(i));
    perm(i) = j;
    classSA_hard(find(class==j)) = i;
end

class_out = classSA_hard;
P_hat = P_hat_best(perm, perm);

Output.class_out = class_out;
Output.P_hat = P_hat;
Output.network_out = network_out;
end