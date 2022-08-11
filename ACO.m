%% 蚁群算法求解TSP问题
% ACO 蚁周模型，规定:完成一次路径循环后，蚂蚁才释放信息素。deltaTau=Q/d(i,j)
%% 清空环境变量
clear
clc

%% 导入数据
load citys_data.mat

%% 计算城市间相互距离
D = squareform(pdist(citys));                                           % 第一处

%% 初始化参数
n = size(citys,1);                      %城市数量
m = 50;                                 % 蚂蚁数量
alpha = 1;                              % 信息素重要程度因子
beta = 5;                               % 启发函数重要程度因子
rho = 0.1;                              % 信息素挥发因子
Q = 10;                                % 信息素常系数
Eta = 1./D;                             % 启发因子函数
Tau = ones(n,n)*Q;                      % 信息素矩阵                     第二处
Table = zeros(m,n);                     % 路径记录表
iter = 1;                               % iterations迭代次数初值
iter_max = 200;                         % 最大迭代次数 
Route_best = zeros(iter_max,n);         % 各代最佳路径       
distance_best = zeros(iter_max,1);      % 前iter次中的最佳路径的长度
iter_distance_best = zeros(iter_max,1); % 各代最佳路径的长度
distance_mean = zeros(iter_max,1);      % 各代路径的平均长度

%% 迭代寻找最佳路径
tic
while iter <= iter_max
    % 随机产生各个蚂蚁的起点城市
      start = zeros(m,1);
      for i = 1:m
          start(i) = randperm(n,1);                                    % 第三处
      end
      Table(:,1) = start; 
      % 构建解空间
      citys_index = 1:n;
      
      % 逐个蚂蚁路径选择
      for i = 1:m
          % 逐个城市路径选择
         for j = 2:n
             tabu = Table(i,1:(j - 1));           % 已访问的城市集合(禁忌表)
             allow_index = ~ismember(citys_index,tabu);
             allow = citys_index(allow_index);  % 待访问的城市集合
%              allow = setdiff(citys_index,tabu);   %该函数效率太慢了
             P = allow;
             % 计算城市间转移概率
             for k = 1:length(allow)
                 P(k)   = Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;
             end
             P = P/sum(P);
             % 轮盘赌法选择下一个访问城市
             Pc = cumsum(P);     
            target_index = find(Pc >= rand);    %返回的是Pc中大于0-1的随机数的元素的位置
            target = allow(target_index(1));
            Table(i,j) = target;    %下一个目标城市
         end
      end
      
      % 计算各个蚂蚁的路径距离
      distance = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);    %各个蚂蚁路径
          for j = 1:(n - 1)
              distance(i) = distance(i) + D(Route(j),Route(j + 1));
          end
          distance(i) = distance(i) + D(Route(n),Route(1));  %加上回到原点的距离
      end
      
      % 计算最短路径和距离及各次迭代平均距离
      [min_distance,min_index] = min(distance);                        %第四处
      iter_distance_best(iter) = min_distance;
      distance_mean(iter) = mean(distance);
      if iter == 1
          distance_best(iter) = min_distance;
          Route_best(iter,:) = Table(min_index,:);
      else
          distance_best(iter) = min(distance_best(iter - 1),min_distance);
          if distance_best(iter) == min_distance
              Route_best(iter,:) = Table(min_index,:);
          else
              Route_best(iter,:) = Route_best((iter-1),:);
          end
      end
      disp(['第',num2str(iter),'迭代最短距离为:',num2str(iter_distance_best(iter)),...
          '循环最短距离为:',num2str(distance_best(iter))])
      
      % 更新信息素
      Delta_Tau = zeros(n,n);
      % 逐个蚂蚁计算
      for i = 1:m
          % 逐个城市计算
          for j = 1:(n - 1)
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/distance(i);
          end
          Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/distance(i);
      end
      Tau = (1-rho) * Tau + Delta_Tau;
      
      %信息素显示
      maxTau = max(Tau(:));
      minTau = min(Tau(:));
      
      tau_normalized = (Tau - minTau) ./ (maxTau - minTau);
      figure(3)
      for i = 1 : n-1
        for j = i + 1 : n
            x1 = citys(i, 1);
            y1 = citys(i, 2);
            
            x2 = citys(j, 1);
            y2 = citys(j, 2);         
            
            Tau_X = [x1, x2];
            Tau_Y = [y1, y2];
            
            plot(Tau_X, Tau_Y, 'color',[0, 0, (1-tau_normalized(i,j)), tau_normalized(i,j)], 'lineWidth',tau_normalized(i,j)+1)
        end
      end
      
      for i = 1 : n
          hold on
          X = citys( : , 1) ;
          Y = citys( : , 2) ;
          plot(X, Y,  'ok',  'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [1 .6 .6])
      end
      
      for i = 1:size(citys,1)
          text(citys(i,1),citys(i,2),['   ' num2str(i)]);
      end
      xlabel('城市位置横坐标')
      ylabel('城市位置纵坐标')
      title('路径上的信息素变化')
      box on
      drawnow
      
    % 迭代次数加1，清空路径记录表
    iter = iter + 1;
    Table = zeros(m,n);
end
toc

%% 结果显示
[Shortest_Length,index] = min(distance_best);
Shortest_Route = Route_best(index,:);
disp(['最短距离:' num2str(Shortest_Length)]);
disp(['最短路径:' num2str([Shortest_Route Shortest_Route(1)])]);

%% 绘图
figure(1)
best_X=[citys(Shortest_Route,1);citys(Shortest_Route(1),1)];
best_Y=[citys(Shortest_Route,2);citys(Shortest_Route(1),2)];
plot(best_X,best_Y,'o-', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [1 .6 .6]);
grid on
for i = 1:size(citys,1)
    text(citys(i,1),citys(i,2),['   ' num2str(i)]);
end
text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       起点');
text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       终点');
xlabel('城市位置横坐标')
ylabel('城市位置纵坐标')
title(['蚁群算法优化路径(最短距离:' num2str(Shortest_Length) ')'])

figure(2)
plot(1:iter_max,distance_best,'b',1:iter_max,distance_mean,'r:',1:iter_max,iter_distance_best,'g-.')
legend('最短距离','平均距离','各代最短距离')
xlabel('迭代次数')
ylabel('距离')
title('最短距离、平均距离以及各代最短距离对比')

