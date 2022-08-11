%% ��Ⱥ�㷨���TSP����
% ACO ����ģ�ͣ��涨:���һ��·��ѭ�������ϲ��ͷ���Ϣ�ء�deltaTau=Q/d(i,j)
%% ��ջ�������
clear
clc

%% ��������
load citys_data.mat

%% ������м��໥����
D = squareform(pdist(citys));                                           % ��һ��

%% ��ʼ������
n = size(citys,1);                      %��������
m = 50;                                 % ��������
alpha = 1;                              % ��Ϣ����Ҫ�̶�����
beta = 5;                               % ����������Ҫ�̶�����
rho = 0.1;                              % ��Ϣ�ػӷ�����
Q = 10;                                % ��Ϣ�س�ϵ��
Eta = 1./D;                             % �������Ӻ���
Tau = ones(n,n)*Q;                      % ��Ϣ�ؾ���                     �ڶ���
Table = zeros(m,n);                     % ·����¼��
iter = 1;                               % iterations����������ֵ
iter_max = 200;                         % ���������� 
Route_best = zeros(iter_max,n);         % �������·��       
distance_best = zeros(iter_max,1);      % ǰiter���е����·���ĳ���
iter_distance_best = zeros(iter_max,1); % �������·���ĳ���
distance_mean = zeros(iter_max,1);      % ����·����ƽ������

%% ����Ѱ�����·��
tic
while iter <= iter_max
    % ��������������ϵ�������
      start = zeros(m,1);
      for i = 1:m
          start(i) = randperm(n,1);                                    % ������
      end
      Table(:,1) = start; 
      % ������ռ�
      citys_index = 1:n;
      
      % �������·��ѡ��
      for i = 1:m
          % �������·��ѡ��
         for j = 2:n
             tabu = Table(i,1:(j - 1));           % �ѷ��ʵĳ��м���(���ɱ�)
             allow_index = ~ismember(citys_index,tabu);
             allow = citys_index(allow_index);  % �����ʵĳ��м���
%              allow = setdiff(citys_index,tabu);   %�ú���Ч��̫����
             P = allow;
             % ������м�ת�Ƹ���
             for k = 1:length(allow)
                 P(k)   = Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;
             end
             P = P/sum(P);
             % ���̶ķ�ѡ����һ�����ʳ���
             Pc = cumsum(P);     
            target_index = find(Pc >= rand);    %���ص���Pc�д���0-1���������Ԫ�ص�λ��
            target = allow(target_index(1));
            Table(i,j) = target;    %��һ��Ŀ�����
         end
      end
      
      % ����������ϵ�·������
      distance = zeros(m,1);
      for i = 1:m
          Route = Table(i,:);    %��������·��
          for j = 1:(n - 1)
              distance(i) = distance(i) + D(Route(j),Route(j + 1));
          end
          distance(i) = distance(i) + D(Route(n),Route(1));  %���ϻص�ԭ��ľ���
      end
      
      % �������·���;��뼰���ε���ƽ������
      [min_distance,min_index] = min(distance);                        %���Ĵ�
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
      disp(['��',num2str(iter),'������̾���Ϊ:',num2str(iter_distance_best(iter)),...
          'ѭ����̾���Ϊ:',num2str(distance_best(iter))])
      
      % ������Ϣ��
      Delta_Tau = zeros(n,n);
      % ������ϼ���
      for i = 1:m
          % ������м���
          for j = 1:(n - 1)
              Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1)) + Q/distance(i);
          end
          Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1)) + Q/distance(i);
      end
      Tau = (1-rho) * Tau + Delta_Tau;
      
      %��Ϣ����ʾ
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
      xlabel('����λ�ú�����')
      ylabel('����λ��������')
      title('·���ϵ���Ϣ�ر仯')
      box on
      drawnow
      
    % ����������1�����·����¼��
    iter = iter + 1;
    Table = zeros(m,n);
end
toc

%% �����ʾ
[Shortest_Length,index] = min(distance_best);
Shortest_Route = Route_best(index,:);
disp(['��̾���:' num2str(Shortest_Length)]);
disp(['���·��:' num2str([Shortest_Route Shortest_Route(1)])]);

%% ��ͼ
figure(1)
best_X=[citys(Shortest_Route,1);citys(Shortest_Route(1),1)];
best_Y=[citys(Shortest_Route,2);citys(Shortest_Route(1),2)];
plot(best_X,best_Y,'o-', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', [1 .6 .6]);
grid on
for i = 1:size(citys,1)
    text(citys(i,1),citys(i,2),['   ' num2str(i)]);
end
text(citys(Shortest_Route(1),1),citys(Shortest_Route(1),2),'       ���');
text(citys(Shortest_Route(end),1),citys(Shortest_Route(end),2),'       �յ�');
xlabel('����λ�ú�����')
ylabel('����λ��������')
title(['��Ⱥ�㷨�Ż�·��(��̾���:' num2str(Shortest_Length) ')'])

figure(2)
plot(1:iter_max,distance_best,'b',1:iter_max,distance_mean,'r:',1:iter_max,iter_distance_best,'g-.')
legend('��̾���','ƽ������','������̾���')
xlabel('��������')
ylabel('����')
title('��̾��롢ƽ�������Լ�������̾���Ա�')

