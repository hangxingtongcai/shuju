function [Gb_Fit,Gb_Sol,Conv_curve]=ICPO(Pop_size,Tmax,lb,ub,dim,fobj)

%%%%-------------------Definitions--------------------------%%
Conv_curve=zeros(1,Tmax);
ub=ub.*ones(1,dim);
lb=lb.*ones(1,dim);

%%-------------------Controlling parameters--------------------------%%
N=Pop_size; %% Initial population size
N_min=round(0.8*Pop_size); %% Minimum population size
T=2; %% Number of cycles
Tf=0.8; %% Tradeoff percentage

%%----------------- Stagnation detection parameters -------------------%%
Stagnation_Limit = 10;
epsilon = 1e-6;
stagnation_counter = 0;
previous_best = inf;

%%---------------Initialization (with ROBL)----------------------%%
X=initialization(Pop_size,dim,ub,lb);
for i=1:Pop_size
    J=rand(1,dim);
    X_opposite(i,:) = J.*(ub + lb - X(i,:));
    f1 = fobj(X(i,:));
    f2 = fobj(X_opposite(i,:));
    if f2 < f1
        X(i,:) = X_opposite(i,:);
    end
    fitness(i) = fobj(X(i,:));
end

[Gb_Fit,index]=min(fitness);
Gb_Sol=X(index,:);    
Xp=X;

t=0;
while t<=Tmax
    r2=rand;
    alpha = 0.2 * (0.5 * (1 + cos(t / Tmax * pi)));

    for i=1:Pop_size
        U1=rand(1,dim)>rand;
        if rand<rand
            if rand<rand
                y=(X(i,:)+X(randi(Pop_size),:))/2;
                X(i,:)=X(i,:)+(randn).*abs(2*rand*Gb_Sol-y);
            else
                y=(X(i,:)+X(randi(Pop_size),:))/2;
                X(i,:)=(U1).*X(i,:)+(1-U1).*(y+rand*(X(randi(Pop_size),:)-X(randi(Pop_size),:)));
            end
        else
            Yt=2*rand*(1-t/(Tmax))^(t/(Tmax));
            U2=rand(1,dim)<0.5*2-1;
            S=rand*U2;
            if rand<Tf
                St=exp(fitness(i)/(sum(fitness)+eps));
                S=S.*Yt.*St;
                X(i,:)= (1-U1).*X(i,:)+U1.*(X(randi(Pop_size),:)+St*(X(randi(Pop_size),:)-X(randi(Pop_size),:))-S); 
            else
                Mt=exp(fitness(i)/(sum(fitness)+eps));
                vt=X(i,:);
                Vtp=X(randi(Pop_size),:);
                Ft=rand(1,dim).*(Mt*(-vt+Vtp));
                S=S.*Yt.*Ft;
                X(i,:)= (Gb_Sol+(alpha*(1-r2)+r2)*(U2.*Gb_Sol-X(i,:)))-S; 
            end
        end

        for j=1:size(X,2)
            if  X(i,j)>ub(j)
                X(i,j)=lb(j)+rand*(ub(j)-lb(j));
            elseif  X(i,j)<lb(j)
                X(i,j)=lb(j)+rand*(ub(j)-lb(j));
            end
        end

        nF=fobj(X(i,:));
        if  fitness(i)<nF
            X(i,:)=Xp(i,:);
        else
            Xp(i,:)=X(i,:);
            fitness(i)=nF;
            if  fitness(i)<=Gb_Fit
                Gb_Sol=X(i,:);
                Gb_Fit=fitness(i);
            end
        end
    end

    % -------- Stagnation Detection and ROBL Injection --------
    if abs(Gb_Fit - previous_best) < epsilon
        stagnation_counter = stagnation_counter + 1;
    else
        stagnation_counter = 0;
    end
    previous_best = Gb_Fit;

    if stagnation_counter >= Stagnation_Limit
        inject_num = round(0.2 * Pop_size);
        idx_replace = randperm(Pop_size, inject_num);
        X_new = ROBL_injection(inject_num, dim, ub, lb, fobj);
        X(idx_replace, :) = X_new;
        stagnation_counter = 0;
    end

    Pop_size=fix(N_min+(N-N_min)*(1-(rem(t,Tmax/T)/Tmax/T)));
    t=t+1;
    if t>Tmax
        break
    end
    Conv_curve(t)=Gb_Fit;
end
end

% ---------------- Initialization ----------------
function Positions=initialization(SearchAgents_no,dim,ub,lb)
Boundary_no= size(ub,2);
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
else
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end
end

% --------------- ROBL Injection ------------------
function X_new = ROBL_injection(N_new, dim, ub, lb, fobj)
X_base = rand(N_new, dim) .* (ub - lb) + lb;
J = rand(N_new, dim);
X_opposite = J .* (ub + lb - X_base);
X_new = X_base;
for i = 1:N_new
    f1 = fobj(X_base(i,:));
    f2 = fobj(X_opposite(i,:));
    if f2 < f1
        X_new(i,:) = X_opposite(i,:);
    end
end
end
