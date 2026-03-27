%%  闲鱼：深度学习与智能算法
%%  唯一官方店铺：https://mbd.pub/o/autMPAr-aWWbm3BtZw==
%%  微信公众号：强盛机器学习，关注公众号获得更多免费代码！
function [bestcord,bestpos,BestCost] = ABC(nPop,MaxIt,VarMin,VarMax,nVar,fobj)

% MaxIt = 200;                  % 最大迭代次数
% nPop = 100;                   % 蜂群大小
VarSize = [1 nVar];     % 变量矩阵
nOnlooker = nPop;             % 侦察蜂个数
L = round(0.6*nVar*nPop);     % 探索极值限制参数
a = 1;                        % 加速度系数上限


%% Initialization
%% 初始化
% 置空蜜蜂矩阵
empty_bee.Position = [];
empty_bee.Cost = [];
% 初始化蜂群数组
pop = repmat(empty_bee, nPop, 1);
% 初始化最优解
BestSol.Cost = inf;
% 产生初始种群
for i = 1:nPop
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    pop(i).Cost = fobj(pop(i).Position);
    if pop(i).Cost <= BestSol.Cost
        BestSol = pop(i);
    end
end
% 丢解计数器
C = zeros(nPop, 1);
% 保存最优函数值的数组
BestCost = zeros(MaxIt, 1);

%% ABC迭代
for it = 1:MaxIt
    % 引领蜂
    for i = 1:nPop
        % 随机选择不等于i的k
        K = [1:i-1 i+1:nPop];
        k = K(randi([1 numel(K)]));
        % 定义加速度系数
        phi = a*unifrnd(-1, +1, VarSize);
        % 新的蜜蜂位置
        newbee.Position = pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        % 边界处理
        newbee.Position = max(newbee.Position, VarMin);
        newbee.Position = min(newbee.Position, VarMax);
        % 新的蜜蜂函数值
        newbee.Cost = fobj(newbee.Position);
        % 比较
        if newbee.Cost <= pop(i).Cost
            pop(i) = newbee;
        else
            C(i) = C(i)+1;
        end
    end
    % 计算适应度值和选择概率
    F = zeros(nPop, 1);
    MeanCost = mean([pop.Cost]);
    for i = 1:nPop
        % 将函数值转换为适应度
        if pop(i).Cost >= 0
            F(i) = 1/(1+pop(i).Cost);
        else
            F(i) = 1+abs(pop(i).Cost);
        end
    end
    P = F/sum(F);
    % 跟随蜂
    for m = 1:nOnlooker
        % 选择食物源
        i = RouletteWheelSelection(P);
        % 随机选择不等于i的k
        K = [1:i-1 i+1:nPop];
        k = K(randi([1 numel(K)]));
        % 定义加速度系数
        phi = a*unifrnd(-1, +1, VarSize);
        % 新的蜜蜂位置
        newbee.Position = pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        % 边界处理
        newbee.Position = max(newbee.Position, VarMin);
        newbee.Position = min(newbee.Position, VarMax);
        % 新的蜜蜂函数值
        newbee.Cost = fobj(newbee.Position);
        % 比较
        if newbee.Cost <= pop(i).Cost
            pop(i) = newbee;
        else
            C(i) = C(i) + 1;
        end
    end
    % 侦察蜂
    for i = 1:nPop
        if C(i) >= L    % 超出探索极值参数
            maxPos = max(pop(i).Position);
            minPos = min(pop(i).Position);
            for j = 1:numel(pop(i).Position)
                pop(i).Position(j) = minPos+rand*(maxPos-minPos);
            end
            pop(i).Cost = fobj(pop(i).Position);
            C(i) = 0;
        end
    end
    % 更新每轮最优解
    for i = 1:nPop
        if pop(i).Cost <= BestSol.Cost
            BestSol = pop(i);
        end
    end
    % 保存每轮最优解
    BestCost(it) = BestSol.Cost;

end

bestcord = BestSol.Cost;
bestpos = BestSol.Position;
end
function i=RouletteWheelSelection(P)

    r=rand;
    
    C=cumsum(P);
    %%%返回满足r<=C的第一个元素在C中的位置
    i=find(r<=C,1,'first');

end
