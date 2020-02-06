clc;
clear;
close all;
 [data,~,raw] = xlsread('C:\Users\Jace\Desktop\ʵϰ\�ӻ�\1010\��ʱ�쳣ֵ����\2_�ܽ�����Ŀ_���ֶ�����ʱ��������_20190907_����.xlsx','����');
% [data,~,raw] = csvread('C:\Users\Jace\Desktop\ʵϰ\�ӻ�\1010\��ʱ�쳣ֵ����\2_�ܽ�����Ŀ_���ֶ�����ʱ��������_20190907_���.csv');
% ����Ԥ����
x1=data(:,1:8);x2=data(:,10:13);DO=data(:,9);                              %���Ա�������������ݷֱ���ȡ����
index1=raw(1,3:10);index2=raw(1,12:15);                                     %���Ա��������������ֱ���ȡ����
data1=[x1 x2];                                                             %�ϲ������Ա�������
TITLE=[index1 index2];                                                     %�ϲ������Ա�������
[C,IA,IC] = unique(raw(2:end,1),'stable');                                 %��ȡ����Ψһֵ������CΪ��������Ψһֵ��IAΪC��ԭ�����е�λ�ã���һ�����ֵ�λ�ã���ICΪC��IA�е�λ��
for ii = 1:length(C)
    index_DM=find(IC==ii);                                                 
    Xdata_DM = data1(index_DM,:);
    Ydata_DM = DO(index_DM,:);
    [m,n] = size(Xdata_DM); %��ȡ����ά��
    mic=cell(n,1);
    for i=1:n
        x=Xdata_DM(:,i)';                                                   %�Ա������е��ã���Ϊmine����������������Ϊ������
        [mx,nx]=find(isnan(x));    % �ҳ�NaN����λ��
        x(nx)=[];
        y=Ydata_DM';
        y(nx)=[];                                                          %��ȡ�ܽ�������
        [my,ny]=find(isnan(y));
        x(ny)=[];
        y(ny)=[];
        minestats = mine(x,y);                                            %���������Ϣϵ��
        mic{i}=minestats.mic;                                              %����ȡ���е�mic����
    end
        %���������Ϣϵ����������ɸѡ��������ǰ����������������Ӽ�XX1
        mic_variable=[TITLE' mic];
        [newtemp ind] = sortrows(mic_variable(:,2),-1);
        XX1 = mic_variable(ind,:); 
        mic_min=0.2;                                                         %�ٶ�ָ�걣����
        for k1=1:n
            if XX1{k1,2} < mic_min
                XX1_num=k1-1;
                break
            end
        end
        %����XX2���ݼ���MICֵ������ֵ
        XX2=cell(XX1_num,m+1);
        XX2(:,1)=XX1(1:XX1_num,1);
        for j=1:XX1_num
            [bool,inx]=find(strcmp(raw,XX2(j,1)));
            XX2(j,2:end)=raw(index_DM+1,inx);                              %��ȡ���Ӽ�XX1������
        end
        variablename2=XX2(:,1);
        MIC_matrix = zeros(XX1_num,XX1_num);
        data2=cell2mat(XX2(:,2:end));
        for kk=1:XX1_num
            for jj=1:XX1_num
                x2=data2(kk,:);
                y2=data2(jj,:);
                minestats = mine(x2,y2);                                   %���������Ϣϵ��
                MIC_matrix(kk,jj)=minestats.mic;                           %����ȡ���е�mic����    
            end
        end
        %% PCA���ɷַ���
        [SA,minp,maxp]=premnmx(data2);                                     %���ݱ�׼���任
        [a,b]=size(data2);                                                 %�õ������ݾ��������������
        [V,D]=eig(MIC_matrix);                                             %����CM������ֵ����������
        for j1=1:a
            DS(1,j1)=D(a+1-j1,a+1-j1);                                     %������ֵ���������е�DS��
        end
        %���㹱����
        for i1=1:a
            DS(2,i1)=DS(1,i1)/sum(DS(1,:));%����������
            DS(3,i1)=sum(DS(1,1:i1))/sum(DS(1,:));%�ۼƹ�����
        end
        %�ٶ����ɷֵ���Ϣ������
        T=0.90;
        for k2=1:a
            if DS(3,k2) >= T
                com_num=k2;
                break;
            end
        end
        %��ȡ���ɷֵ���������
        for j2=1:com_num
            PV(j2,:)=V(a+1-j2,:);
        end
        %�������ɷֵ÷�
        [Y1]=postmnmx(SA,minp,maxp);
        new_score=PV*Y1;
        for i2=1:b
            total_score(1,i2)=sum(new_score(:,i2));
            total_score(2,i2)=i2;
        end
        %ǿ���ɷֵ÷����ַܷŵ�ͬһ��������
        result_report=[new_score;total_score];
        %���ֽܷ�������
        result_report=sortrows(result_report',-(com_num+1));
        
        % ������
        resultl = cell(3,a+1);
        resultl(:,1) = {'����ֵ','������','�ۻ�������'};
        resultl(1:3,2:end) = num2cell(DS);
        resultl
        
        %ȡǰ���ɷ�
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
    xlswrite('C:\Users\Jace\Desktop\ʵϰ\�ӻ�\1010\��ʱ�쳣ֵ����\2_�ܽ�����Ŀ_���ֶ�����ʱ��������_20190907_����.xlsx',mic_DM{1,i},char(C(i,1)));
end

% xlswrite('C:\Users\Lwjing\Desktop\MOD-SVM\PCA-MIC.xls',result2,'���ɷ�');
% xlswrite('C:\Users\Lwjing\Desktop\MOD-SVM\PCA-MIC.xls',XX1,'MICֵ');
% xlswrite('C:\Users\Lwjing\Desktop\MOD-SVM\PCA-MIC.xls',Y1','�÷�');
% xlswrite('C:\Users\Lwjing\Desktop\MOD-SVM\PCA-MIC.xls',resultl','������');
% disp('���ɷֵ÷ּ����򣨰���4�е��ֽܷ��н������У�ǰ4��Ϊ�����ɷֵ÷֣���6��Ϊָ���ţ�')
% result_report
% %% box plot forratings data
% % To get a quickimpression of the ratings data, make a box plot
% figure;
% boxplot(data2','orientation','horizontal','labels',variablename2);
% grid on;
% %% ͨ����ͼ���Կ���ǰ7�����ɷֿ��Ա�ʾ��ԭʼ���ݵ�90%.
% figure;
% percent_explained= 100*DS(2,:); %cumsum(latent)./sum(latent)
% pareto(percent_explained);
% xlabel('PrincipalComponent');
% ylabel('VarianceExplained (%)');
% %% �鿴ǰ2�����ɷֵ����ɷַ���ͼ
% biplot(PV(1:2,:)','Score',new_score(1:2,:)','VarLabels',variablename2)
% title('ǰ�������ɷֵ����ɷַ���ͼ')
% % ����������������۲�ֵ���ɷ�ֵ�����άͼ����������������ɷ�1��ֵ��������������ɷ�2��ֵ�������Ǿ���һ��������
% % ��������ߴ�������غ�ֵ�������������������ƣ��������ں������ͶӰ���Ƕ����ɷ�1���غɣ�ϵ���������������ͶӰ���Ƕ����ɷ�2���غɡ���ֵ�ǻ�ȡ�෴����
% %%  Ȼ�������о�������������ϴ�ʵ�������бȽϡ�
% % ��ȡǰ4�����ɷֵ÷�
% scoredata = new_score(1:2,:)';
% % ��Ĵ�СΪ1.5
% figure()
% name = {'��һ���ɷֵ÷�','�ڶ����ɷֵ÷�'}
% gplotmatrix(scoredata,[],[],'r','o',1.5,'off',[],name,name)
% %����ĵ� i �С��� j ���е���ͼ�� X �ĵ� i ������� X �ĵ� j �е�ɢ��ͼ���ضԽ��߷����� X ��ÿһ�е�ֱ��ͼ
% title('���ɷֵ÷�ɢ��ͼ')