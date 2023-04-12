function ExperimentalDesign_SBS
clc
clear
close all
warning off
 
Need_Experiment   =   1;
Need_CalMetrics   =   0;
Need_Statistics     =   0;
 
%运行次数
Run_Times     =    11;
%测试问题
Problems = {...   
    'SharedBS6G_t1',...
    'SharedBS6G_t2',...
    'SharedBS6G_t3',...
    'SharedBS6G_t4',...
    'SharedBS6G_t5',...
    'SharedBS6G_t6',...
    'SharedBS6G_t7',...
    'SharedBS6G_t8',...
    'SharedBS6G_t9'...
%                'SharedBS6G_2000',...
%                 'SharedBS6G_s2',...
%                 'SharedBS6G_m1',...
%                 'SharedBS6G_m2',...
%                 'SharedBS6G_l1',...
%                 'SharedBS6G_l2',...
%                 'SharedBS6G_l3',...
%                 'SharedBS6G_l4',...
%                 'SharedBS6G_l5',...
%                 'SharedBS6G_l6',...
           };
%测试算法
Algorithms = {...

                'DPL_BSP',...
                'stMOBEA'

    };
 
  
%测试目标维数
% Obj_M        = [...
%     1
%     ];
%自变量空间维数
% Dec_D        = [...
%      30,...
%      50,...
%     100
%     ];
%种群规模
Size         = [...
 
%      20,...
%      20,...
%      20,...
%      20,...
%      40,...
%      40,...
%      40,...
%      40,...
%      40,...
%      40,...
50
    ];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pool = {};
if Need_Experiment
    for a=1:length(Algorithms)
        TheAlgorithm = eval(['@',Algorithms{a}]);
        for p = 1:1:length(Problems)
            TheProblem = eval(['@',Problems{p}]);
            for j=1:Run_Times
             

               filename = [Algorithms{a},'_',Problems{p},'_uM2_lM2_uD',num2str(Size(p)),'_lD',num2str(Size(p)),'_',num2str(j),'.mat'];


                if ~exist(filename)
                    Inputs = {'-algorithm',TheAlgorithm,'-problem',TheProblem,'-Nu',Size(p),'-Nl',Size(p),'-run',j,'-save',1};
                    Pool = cat(1,Pool,Inputs);
                end
            end
        end
    end
end

%并行开关
if ~isempty(Pool)
    Isparallel = 1;
    parallelnum = 22;
    if ~Isparallel
        for k=1:size(Pool,1)
            main(Pool{k,:},folder);
        end
    else
        parpool('local',parallelnum)
        parfor k=1:size(Pool,1)
%             Notcompelted = true; 
%             while Notcompelted
%                 try
                    main(Pool{k,:});
%                     Notcompelted =false;
%                 catch
%                     
%                 end
%             end
        end
        delete(gcp);
    end
end

% if Need_Experiment
%  
%     for a=1:length(Algorithms)
%         TheAlgorithm = eval(['@',Algorithms{a}]);
%         for p = 1:1:length(Problems)
%             TheProblem = eval(['@',Problems{p}]);
% %             for j = 1:1:length(Obj_M)
% %                 for i = 1:1:length(Dec_D)
% %                     M   = Obj_M(j);
% %                     D   = Dec_D(i);
% %                     N   = Size(i);                    
%                     %             Gen = Max_Gen(p,j);
%                       for j = 1:Run_Times/7
%                        parpool('local',7);
%                        parfor k = 1:7                                   
%                         main('-algorithm',TheAlgorithm,'-problem',TheProblem,'-Nu',50,'-Nl',50,'-run',(j-1)*7+k,'-save',1);
%                        end
%                          delete(gcp);
%                       end
% %                 end
% %             end
%         end
%     end
%  
% end
 
%% 计算指标
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Metrics = {...
%     'runtime'...
%     'GD'...
%     'HV'...
%     'IGD'...
%     'NHV'...    
%     'Spread'...
%     'Coverage'...
%     'DM'...
%     'PD'...
%     'Spacing'...
    };
 
 
if Need_CalMetrics 
    
    main
    close
 
    for m=1 : length(Metrics)
        TheMetric = Metrics{m};
    for j = 1 : 1:length(Obj_M)
        M = Obj_M(j);
        for a=1 : length(Algorithms)
            TheAlgorithm = Algorithms{a};
            for p = 1 : 1 : length(Problems)
                TheProblem = Problems{p};
                parfor r=1 : Run_Times
                    CalMetrics(TheAlgorithm,TheProblem,M,r,TheMetric);
                end
            end
        end
    end
    end
 
end
 
 
if Need_Statistics
 
%% 统计实验结果
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
main
close
% 说明实验的指标
Metrics={...
     'runtime'...
     'GD'...
     'HV'...
     'IGD'...
     'NIGD'...
     'IGD_epsilon'...
     'NHV'...    
%     'Spread'...
%     'Coverage'...
%     'DM'...
%     'PD'...
%     'Spacing'...
    };
 
%   需要统计的统计量
Statistics={...
    'mean'...       %均值
%     'var'...        %方差
%     'std'...        %标准差
%     'median'...     %中位数
%    'min'...
    };
 
if exist('Data','dir')==0
   mkdir('Data');
end
 
for s=1:length(Statistics)
    for m=1:length(Metrics)
        for j = 1:1:length(Obj_M)
            Matrix=[];
            for a=1:length(Algorithms)
                vector.([Metrics{m},'_M_',Statistics{s}])=[];
                for p = 1:1:length(Problems)
                    DataToUse.(Metrics{m})=[];
                    for r=1:Run_Times
                        Filename=['Data\',Algorithms{a},'\',Algorithms{a},'_',Problems{p},'_M',num2str(Obj_M(j)),'_D',num2str(Obj_D(p)),'_',num2str(r),'.mat'];
                        Data=load(Filename);
                        metric = Data.metric;             
                        %DataToUse.IGD = [ DataToUse.IGD,metric.IGD ];  
                        DataToUse.(Metrics{m}) =[ DataToUse.(Metrics{m}),metric.(Metrics{m})];          
                    end
                    %IGD_M_mean(p,a) = mean(IGD);   %????????eval?????????ü??
                    Statistics_func = eval(['@',Statistics{s}]);
                    Statistics_ind = Statistics_func(DataToUse.(Metrics{m}));
                    vector.([Metrics{m},'_M_',Statistics{s}]) = [vector.([Metrics{m},'_M_',Statistics{s}]);Statistics_ind];
                end
                Matrix = [Matrix, vector.([Metrics{m},'_M_',Statistics{s}])];
            end
        
            file_path=['Data\',Metrics{m},'.xlsx'];
            Sheet_Name=['M',num2str(Obj_M(j)),'_',Statistics{s}];
            xlswrite(file_path,{Sheet_Name}, Sheet_Name,'A1');
            xlswrite(file_path, Problems',   Sheet_Name,'A2');
            xlswrite(file_path, Algorithms,  Sheet_Name,'B1');
            xlswrite(file_path, Matrix,      Sheet_Name,'B2');
        end
    end
end
 
end
 
end


