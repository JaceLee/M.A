clc;
clear;
close all;
 [data,~,raw] = xlsread('C:\Users\Jace\Desktop\实习\杂活\1010\逐时异常值处理\2_溶解氧项目_部分断面逐时数据整理_20190907_汇总.xlsx','汇总');
% [data,~,raw] = csvread('C:\Users\Jace\Desktop\实习\杂活\1010\逐时异常值处理\2_溶解氧项目_部分断面逐时数据整理_20190907_大墩.csv');
% 数据预处理
x1=data(:,1:8);x2=data(:,10:13);DO=data(:,9);                              %把自变量和因变量数据分别提取出来
index1=raw(1,3:10);index2=raw(1,12:15);                                     %把自变量和因变量标题分别提取出来
data1=[x1 x2];                                                             %合并所有自变量数据
TITLE=[index1 index2];                                                     %合并所有自变量标题
[C,IA,IC] = unique(raw(2:end,1),'stable');                                 %获取断面唯一值不排序，C为断面名字唯一值，IA为C在原数据中的位置（第一个出现的位置），IC为C在IA中的位置
for ii = 1:length(C)
    index_DM=find(IC==ii);                                                 
    Xdata_DM = data1(index_DM,:);
    Ydata_DM = DO(index_DM,:);
    [m,n] = size(Xdata_DM); %获取数据维度
    mic=cell(n,1);
    for i=1:n
        x=Xdata_DM(:,i)';                                                   %自变量行列倒置，因为mine函数的输入数据需为列数据
        [mx,nx]=find(isnan(x));    % 找出NaN数据位置
        x(nx)=[];
        y=Ydata_DM';
        y(nx)=[];                                                          %提取溶解氧序列
        [my,ny]=find(isnan(y));
        x(ny)=[];
        y(ny)=[];
        minestats = mine(x,y);                                            %计算最大信息系数
        mic{i}=minestats.mic;                                              %仅提取其中的mic数据
    end
        %根据最大信息系数进行排序，筛选出排名靠前的若干因子组成因子集XX1
        mic_variable=[TITLE' mic];
        [newtemp ind] = sortrows(mic_variable(:,2),-1);
        XX1 = mic_variable(ind,:); 
        mic_min=0.2;                                                         %假定指标保留率
        for k1=1:n
            if XX1{k1,2} < mic_min
                XX1_num=k1-1;
                break
            end
        end
        %计算XX2数据集的MIC值及特征值
        XX2=cell(XX1_num,m+1);
        XX2(:,1)=XX1(1:XX1_num,1);
        for j=1:XX1_num
            [bool,inx]=find(strcmp(raw,XX2(j,1)));
            XX2(j,2:end)=raw(index_DM+1,inx);                              %提取因子集XX1的数据
        end
        variablename2=XX2(:,1);
        MIC_matrix = zeros(XX1_num,XX1_num);
        data2=cell2mat(XX2(:,2:end));
        for kk=1:XX1_num
            for jj=1:XX1_num
                x2=data2(kk,:);
                y2=data2(jj,:);
                minestats = mine(x2,y2);                                   %计算最大信息系数
                MIC_matrix(kk,jj)=minestats.mic;                           %仅提取其中的mic数据    
            end
        end
        %% PCA主成分分析
        [SA,minp,maxp]=premnmx(data2);                                     %数据标准化变换
        [a,b]=size(data2);                                                 %得到的数据矩阵的行数和列数
        [V,D]=eig(MIC_matrix);                                             %计算CM的特征值和特征向量
        for j1=1:a
            DS(1,j1)=D(a+1-j1,a+1-j1);                                     %将特征值按降序排列到DS中
        end
        %计算贡献率
        for i1=1:a
            DS(2,i1)=DS(1,i1)/sum(DS(1,:));%单个贡献率
            DS(3,i1)=sum(DS(1,1:i1))/sum(DS(1,:));%累计贡献率
        end
        %假定主成分的信息保留率
        T=0.90;
        for k2=1:a
            if DS(3,k2) >= T
                com_num=k2;
                break;
            end
        end
        %提取主成分的特征向量
        for j2=1:com_num
            PV(j2,:)=V(a+1-j2,:);
        end
        %计算主成分得分
        [Y1]=postmnmx(SA,minp,maxp);
        new_score=PV*Y1;
        for i2=1:b
            total_score(1,i2)=sum(new_score(:,i2));
            total_score(2,i2)=i2;
        end
        %强主成分得分与总分放到同一个矩阵中
        result_report=[new_score;total_score];
        %按总分降序排列
        result_report=sortrows(result_report',-(com_num+1));
        
        % 输出结果
        resultl = cell(3,a+1);
        resultl(:,1) = {'特征值','贡献率','累积贡献率'};
        resultl(1:3,2:end) = num2cell(DS);
        resultl
        
        %取前主成分
        result2 = cell(com_num+1,a);
        result2(1,:) = variablename2';
        result2(2:end,:) = num2cell(PV);
        result2
        
        mic_DM{ii} = XX1;
        mic_DM_TQ{ii} = XX2;
        MIC_matrix_DM{ii} = MIC_matrix;
        GX_DM{ii} =resultl;
        DF_DM{ii} =result2;
        TS_DM{ii} =result_report;
        clear SA minp maxp DS PV Y1 new_score result_report resultl result2 total_score
end

for i =1:length(C)
    xlswrite('C:\Users\Jace\Desktop\实习\杂活\1010\逐时异常值处理\2_溶解氧项目_部分断面逐时数据整理_20190907_汇总.xlsx',mic_DM{1,i},char(C(i,1)));
end

% xlswrite('C:\Users\Lwjing\Desktop\MOD-SVM\PCA-MIC.xls',result2,'主成分');
% xlswrite('C:\Users\Lwjing\Desktop\MOD-SVM\PCA-MIC.xls',XX1,'MIC值');
% xlswrite('C:\Users\Lwjing\Desktop\MOD-SVM\PCA-MIC.xls',Y1','得分');
% xlswrite('C:\Users\Lwjing\Desktop\MOD-SVM\PCA-MIC.xls',resultl','贡献率');
% disp('主成分得分及排序（按第4列的总分进行降序排列，前4列为个各成分得分，第6列为指标编号）')
% result_report
% %% box plot forratings data
% % To get a quickimpression of the ratings data, make a box plot
% figure;
% boxplot(data2','orientation','horizontal','labels',variablename2);
% grid on;
% %% 通过看图可以看出前7个主成分可以表示出原始数据的90%.
% figure;
% percent_explained= 100*DS(2,:); %cumsum(latent)./sum(latent)
% pareto(percent_explained);
% xlabel('PrincipalComponent');
% ylabel('VarianceExplained (%)');
% %% 查看前2个主成分的主成分分析图
% biplot(PV(1:2,:)','Score',new_score(1:2,:)','VarLabels',variablename2)
% title('前两个主成分的主成分分析图')
% % 红点代表二个或三个观测值主成分值，如二维图，红点横坐标就是主成分1的值，纵坐标就是主成分2的值。但都是经过一定的量化
% % 蓝点和蓝线代表的是载荷值，坐标与红点意义相类似，所以线在横坐标的投影就是对主成分1的载荷（系数），在纵坐标的投影就是对主成分2的载荷。其值是或取相反数了
% %%  然后对其进行聚类分析，并与上次实验结果进行比较。
% % 提取前4个主成分得分
% scoredata = new_score(1:2,:)';
% % 点的大小为1.5
% figure()
% name = {'第一主成分得分','第二主成分得分'}
% gplotmatrix(scoredata,[],[],'r','o',1.5,'off',[],name,name)
% %矩阵的第 i 行、第 j 列中的子图是 X 的第 i 列相对于 X 的第 j 列的散点图。沿对角线方向是 X 的每一列的直方图
% title('主成分得分散点图')