clc
close all
clear
% Initialization
allreflag1 = [];
allreflag2 = [];
allreflag3 = [];
allcrflag1 = [];
allcrflag2 = [];
allcrflag3 = [];
PortPosition = [0 50];
MaxRobotNo = 3;
global XWind YWind RWind
XWind = 69.73;
YWind = 58.15;
RWind = 15;
load('x.mat') 
load('y.mat')
load('r.mat')
XVessel = x(1:6);
YVessel = y(1:6);
RVessel = r(1:6); 
reflag1 = [];
reflag2 = [];
reflag3 = [];
crflag1 = []; 
crflag2 = [];
crflag3 = [];
DistanceWeight = zeros(1,length(XVessel)); 
for i = 1:length(XWind)
    distance = sqrt((XVessel-XWind(i)).^2+(YVessel-YWind(i)).^2);
    DistanceWeight(find(distance<RWind(i))) = DistanceWeight(find(distance<RWind(i))) + distance(find(distance<RWind(i)))./RWind(i);
end
GenotypeLength = [4 4 2 2];
PopulationSize = 100; 
MaxGeneration = 100; 
MemorySize = 5;
ArtiSize = 10;
Pc = 0.9;
Pm = 0.1;
N = PopulationSize/2+10;
VarietyNumber = length(XVessel); 
FinalPop = [];
AllPop = [];
UCB_std = 0;
cu_re_1 = cell(1, sum(GenotypeLength)*VarietyNumber);
cu_re_0 = cell(1, sum(GenotypeLength)*VarietyNumber);
dis_factor = 0.9;
InitPop = randi([0 1],PopulationSize,VarietyNumber*sum(GenotypeLength)+VarietyNumber*4);
DecodePop = DecodeFunction(InitPop,XVessel,YVessel,RVessel,GenotypeLength,PopulationSize,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
pa = DecodePop(:,(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
IsDomi = IDAf(pa);
DomiIndex = find(IsDomi==1);
DomiPop = DecodePop(DomiIndex,:);
hillcoff = hillco(DomiPop(:,end-2:end)');
useindex = find(hillcoff>0);
DomiPop = DomiPop(useindex,:);
% whether iteration ends?
iter = 1;
while iter <= MaxGeneration
    newpop = DecodePop;
    % Actor-critic-inspired crossover
    crossPop = [];
    for i = 1:PopulationSize
        index1 = randi(PopulationSize);
        index2 = randi(PopulationSize);
        parent1 = newpop(index1,:);
        parent2 = newpop(index2,:);
        tempPa1 = parent1((sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
        tempGene1 = parent1(1:sum(GenotypeLength)*VarietyNumber);
        tempPa2 = parent2((sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
        tempGene2 = parent2(1:sum(GenotypeLength)*VarietyNumber);
        newIndi = [];
        for j = 1:sum(GenotypeLength)*VarietyNumber
            newGene1 = [tempGene1(1:j-1) 1-tempGene1(j) tempGene1(j+1:end)];
            decodeGene1 = DecodeFunction(newGene1,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
            decodePa1 = decodeGene1((sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
            rela1 =  IDAf([tempPa1; decodePa1]);
            stp1 = find(newpop(:, j) == tempGene1(j));
            stp2 = find(newpop(:, j) == 1 - tempGene1(j));
            stp1Pa = newpop(stp1, (sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
            stp2Pa = newpop(stp2, (sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
            rela =  IDAf([stp1Pa; stp2Pa]);
            if length(stp1) ~= 0
                cureward1 = sum(rela(1:length(stp1))) / length(stp1);
            else
                cureward1 = 0;
            end
            if length(stp2) ~= 0
                cureward2 = sum(rela((length(stp1)+1):PopulationSize)) / length(stp2);
            else
                cureward2 = 0;
            end
            if tempGene1(j) == 1
                cu_re_1{j} = [cu_re_1{j}, cureward1];
                cu_re_0{j} = [cu_re_0{j}, cureward2];
            else
                cu_re_0{j} = [cu_re_0{j}, cureward1];
                cu_re_1{j} = [cu_re_1{j}, cureward2];
            end
            newGene2 = [tempGene2(1:j-1) 1-tempGene2(j) tempGene2(j+1:end)];
            decodeGene2 = DecodeFunction(newGene2,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
            decodePa2 = decodeGene2((sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
            rela2 =  IDAf([tempPa2; decodePa2]);
            if (rela1(1) == rela2(1)) && (rela1(1) == 1) && rand < Pc 
                if tempGene1(j) == tempGene2(j)
                    newIndi = [newIndi tempGene1(j)];
                else 
                    cuallreward1 = 0;
                    cuallreward0 = 0;
                    for ir = 1 : (length(cu_re_1{j}) - 1)
                        cuallreward1 = cuallreward1 + cu_re_1{j}(ir) * dis_factor ^ (length(cu_re_1{j}) - ir);
                        cuallreward0 = cuallreward0 + cu_re_0{j}(ir) * dis_factor ^ (length(cu_re_0{j}) - ir);
                    end
                    if tempGene1(j) == 1
                        cureward1 = cu_re_1{j}(end);
                        cu_re1 = cuallreward1;
                        cureward2 = cu_re_0{j}(end);
                        cu_re2 = cuallreward0;
                    else
                        cureward1 = cu_re_0{j}(end);
                        cu_re1 = cuallreward0;
                        cureward2 = cu_re_1{j}(end);
                        cu_re2 = cuallreward1;
                    end
                    if cureward1 + cu_re1 / iter / MaxGeneration > cureward2 + cu_re2 / iter / MaxGeneration
                        newIndi = [newIndi tempGene1(j)];
                    else
                        newIndi = [newIndi tempGene2(j)];
                    end
                end
            elseif rela1(1) == 1 && rand < Pc
                newIndi = [newIndi tempGene1(j)];
            elseif rela2(1) == 1 && rand < Pc
                newIndi = [newIndi tempGene2(j)];
            else
                rela = IDAf([tempPa1; tempPa2]);
                if rela(2) == 0
                    newIndi = [newIndi tempGene1(j)];
                else
                    newIndi = [newIndi tempGene2(j)];
                end
            end
        end
        crossPop = [crossPop; newIndi];
    end
    DecodePop = DecodeFunction(crossPop,XVessel,YVessel,RVessel,GenotypeLength,PopulationSize,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
    % Upper Confidence Bound-based mutation
    mutatePop = [];
    for i = 1:PopulationSize
        parent = DecodePop(i,:);
        tempPa = parent((sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
        tempGene = parent(1:sum(GenotypeLength)*VarietyNumber);
        newIndi = [];
        for j = 1:sum(GenotypeLength)*VarietyNumber
            if rand < Pm
                newGene = [tempGene(1:j-1) 1-tempGene(j) tempGene(j+1:end)];
                decodeGene = DecodeFunction(newGene,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
                decodePa = decodeGene((sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
                rela =  IDAf([tempPa; decodePa]);
                x_j = rela(2);
                n = PopulationSize;
                nj = sum(newpop(:,j) == tempGene(j));
                UCB = x_j + sqrt(2*log(n)/nj);
                if UCB > UCB_std
                    newIndi = [newIndi 1-tempGene1(j)];
                    UCB_std = UCB;
                else
                    newIndi = [newIndi tempGene1(j)];
                end
            else
                newIndi = [newIndi tempGene1(j)];
            end
        end
        mutatePop = [mutatePop; newIndi];
    end
    DecodePop = DecodeFunction(crossPop,XVessel,YVessel,RVessel,GenotypeLength,PopulationSize,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
    % Mutation process regulation
    renumber1 = 0;
    renumber2 = 0;
    renumber3 = 0;
    repath = [];
    nonrepath = [];
    nonreno = [];
    crnumber1 = 0;
    crnumber2 = 0;
    crnumber3 = 0;
    crpath = [];
    noncrpath = [];
    noncrno = [];
    for itry = 1:PopulationSize
        [flag,positionindex] = trytest2(DecodePop(itry,sum(GenotypeLength)*VarietyNumber+1:end-3),GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition);
        renumber1 = renumber1+flag(1);
        if flag(1) ~= 0 
            temp1 = MultiMutationFunction1(positionindex,VarietyNumber,GenotypeLength,1,DecodePop(itry,1:sum(GenotypeLength)*length(XVessel))); 
            temp1(:,sum(GenotypeLength)*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3) = 0;
            decodetemp = DecodeFunction(temp1,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
            [flag,positionindex] = trytest2(decodetemp(1,sum(GenotypeLength)*VarietyNumber+1:end-3),GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition);
            renumber2 = renumber2+flag(1);
            repath = [repath;DecodePop(itry,:)];
        else
            nonrepath = [nonrepath;DecodePop(itry,:)];
            nonreno = [nonreno itry];
        end
        if flag(2) ~= 0
            crnumber1 = crnumber1+1;
            temp1 = MultiMutationFunction2(VarietyNumber,GenotypeLength,1,DecodePop(itry,1:sum(GenotypeLength)*length(XVessel)));
            temp1(:,sum(GenotypeLength)*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3) = 0;
            decodetemp = DecodeFunction(temp1,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
            [flag,positionindex] = trytest2(decodetemp(1,sum(GenotypeLength)*VarietyNumber+1:end-3),GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition);
            if flag(2) ~= 0
                crnumber2 = crnumber2+1;
            end
            crpath = [crpath;DecodePop(itry,:)];
        else
            noncrpath = [noncrpath;DecodePop(itry,:)];
            noncrno = [noncrno itry];
        end
    end
    newrepath = [];
    for ttry = 1:size(repath,1)
        redun = repath(ttry,:); 
        for t = 1:sum(GenotypeLength)*VarietyNumber
            gene = redun(t);
            mugene = 1-redun(t);
            pb1 = size(nonrepath,1)/PopulationSize;
            pab1 = length(find(nonrepath(:,t)==mugene))/size(nonrepath,1);
            pb2 = size(repath,1)/PopulationSize;
            pab2 = length(find(repath(:,t)==mugene))/size(repath,1);
            pb1a =  (pb1*pab1)/(pb1*pab1+pb2*pab2);
            if rand <= pb1a
                redun(t) = mugene;
            end
        end
        redunnew = DecodeFunction(redun,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
        [flagnew,positionindex] = trytest2(redunnew(1,sum(GenotypeLength)*VarietyNumber+1:end-3),GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition);
        renumber3 = renumber3+flagnew(1);
        newrepath = [newrepath;redunnew];
    end
    newcrpath = [];
    for ttry = 1:size(crpath,1)
        crdun = crpath(ttry,:);
        for t = 1:sum(GenotypeLength)*VarietyNumber
            gene = crdun(t);
            mugene = 1-crdun(t);
            pb1 = size(noncrpath,1)/PopulationSize;
            pab1 = length(find(noncrpath(:,t)==mugene))/size(noncrpath,1);
            pb2 = size(crpath,1)/PopulationSize;
            pab2 = length(find(crpath(:,t)==mugene))/size(crpath,1);
            pb1a =  (pb1*pab1)/(pb1*pab1+pb2*pab2);
            if rand <= pb1a
                crdun(t) = mugene;
            end
        end
        crdunnew = DecodeFunction(crdun,XVessel,YVessel,RVessel,GenotypeLength,1,VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
        [flagnew,positionindex] = trytest2(crdunnew(1,sum(GenotypeLength)*VarietyNumber+1:end-3),GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition);
        if flagnew(2) ~= 0
            crnumber3 = crnumber3+1;
        end
        newcrpath = [newcrpath;crdunnew];
    end
    DecodePop(setdiff([1:PopulationSize],nonreno),:) = newrepath;
    DecodePop(setdiff([1:PopulationSize],noncrno),:) = newcrpath;
    reflag1 = [reflag1 renumber1];
    reflag2 = [reflag2 renumber2];
    reflag3 = [reflag3 renumber3];
    crflag1 = [crflag1 crnumber1];
    crflag2 = [crflag2 crnumber2];
    crflag3 = [crflag3 crnumber3];
    % Evaluation
    pa = DecodePop(:,(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
    % Selection process regulation of memory cell
    IsDomi = IDAf(pa); 
    DomiIndex = find(IsDomi==1);
    DomiPop1=[DomiPop;DecodePop(DomiIndex,:)];
    pa = DomiPop1(:,(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
    IsDomi = IDAf(pa); 
    DomiIndex = find(IsDomi==1); 
    DomiPop2=DomiPop1(DomiIndex,:); 
    hillcoff = hillco(DomiPop2(:,end-2:end)');
    useindex = find(hillcoff>0);
    DomiPop = DomiPop2(useindex,:);
    tempsize = size(DomiPop);
    if tempsize(1) == 0
        stop
    end
    FinalPop = [FinalPop;DomiPop];
    AllPop = [AllPop;DecodePop];
    randorder = randperm(PopulationSize);
    DecodePop = [DomiPop(1:min(size(DomiPop,1),MemorySize),:);DecodePop(randorder(1:PopulationSize-min(size(DomiPop,1),MemorySize)),:)];
    iter = iter+1
end
DomiPop = DecodeFunction(DomiPop,XVessel,YVessel,RVessel,GenotypeLength,size(DomiPop,1),VarietyNumber,MaxRobotNo,DistanceWeight,PortPosition);
Epa = DomiPop(:,(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
load train
sound(y,Fs)
% Output solutions in memory cell
f1 = Epa(:,1);
f2 = Epa(:,2);
f3 = Epa(:,3);
figure(2)
axis([min(f1)-50 max(f1)+50 min(f2)-50 max(f2)+50 min(f3)-0.5 max(f3)+0.5])
hold on
h = plot3(f1,f2,f3,'o','MarkerFaceColor',[214/255, 39/255, 40/255],'MarkerEdgeColor',[214/255, 39/255, 40/255]); grid on;
hXLabel = xlabel('Objective Function 1')
hYLabel = ylabel('Objective Function 2');
hZLabel = zlabel('Objective Function 3');
title('Pareto-optimal Solution Space')
hold on
set(h,'LineWidth',1.5)
set(gca, 'Box', 'on', ...                                
    'LineWidth', 1,...                                       
    'TickDir', 'in', 'TickLength', [.015 .015], ...         
    'XMinorTick', 'off', 'YMinorTick', 'off', ...         
    'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1])    
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set(gca, 'FontSize', 13)
set([hXLabel, hYLabel, hZLabel], 'FontSize', 14)
VesselPro = 1;
VesselPose(1,:) = [-2*VesselPro 0.8*VesselPro];
VesselPose(2,:) = [-VesselPro -0.8*VesselPro];
VesselPose(3,:) = [VesselPro -0.8*VesselPro];
VesselPose(4,:) = [2*VesselPro 0.8*VesselPro];
VesselPose(5,:) = [-2*VesselPro 0.8*VesselPro];
VesselPose1(1,:) = [-VesselPro 5*VesselPro];
VesselPose1(2,:) = [-VesselPro 0.8*VesselPro];
VesselPose1(3,:) = [VesselPro 2.5 * VesselPro];
VesselPose1(4,:) = [-VesselPro 5*VesselPro];
for iplot = 1:size(Epa,1)
    figure
    grid on
    axis([0 100 0 100])
    hold on
    axis square;
    hold on
    for i = 1:length(XWind) 
        theta = 0:pi/50:2*pi;
        x0 = XWind(i)+RWind(i)*cos(theta);
        y0 = YWind(i)+RWind(i)*sin(theta);
        plot(x0,y0,':','Color',[214/255,39/255,40/255],'LineWidth',1.2);
        hold on
        x0 = XWind(i)+12*cos(theta);
        y0 = YWind(i)+12*sin(theta);
        plot(x0,y0,':','Color',[214/255,39/255,40/255],'LineWidth',1.2);
        hold on
        x0 = XWind(i)+9*cos(theta);
        y0 = YWind(i)+9*sin(theta);
        plot(x0,y0,':','Color',[214/255,39/255,40/255],'LineWidth',1.2);
        hold on
        x0 = XWind(i)+6*cos(theta);
        y0 = YWind(i)+6*sin(theta);
        plot(x0,y0,':','Color',[214/255,39/255,40/255],'LineWidth',1.2);
        hold on
        x0 = XWind(i)+3*cos(theta);
        y0 = YWind(i)+3*sin(theta);
        plot(x0,y0,':','Color',[214/255,39/255,40/255],'LineWidth',1.2);
        hold on
        x0 = XWind(i)+0.5*cos(theta);
        y0 = YWind(i)+0.5*sin(theta);
        plot(x0,y0,':','Color',[214/255,39/255,40/255],'LineWidth',1.2);
        hold on
    end
    for i = 1:length(XVessel) 
        fill([VesselPose(1,1)+XVessel(i),VesselPose(2,1)+XVessel(i),VesselPose(3,1)+XVessel(i),VesselPose(4,1)+XVessel(i),VesselPose(5,1)+XVessel(i)], [VesselPose(1,2)+YVessel(i),VesselPose(2,2)+YVessel(i),VesselPose(3,2)+YVessel(i),VesselPose(4,2)+YVessel(i),VesselPose(5,2)+YVessel(i)],[31/255,119/255,180/255])
        theta = 0:pi/50:2*pi;
        x0 = XVessel(i)+RVessel(i)*cos(theta);
        y0 = YVessel(i)+RVessel(i)*sin(theta);
        plot(x0,y0,'-.','Color',[31/255,119/255,180/255],'LineWidth',1.2);
        for ixp = 1:4
            plot([VesselPose(ixp,1)+XVessel(i) VesselPose(ixp+1,1)+XVessel(i)],[VesselPose(ixp,2)+YVessel(i) VesselPose(ixp+1,2)+YVessel(i)],'Color',[31/255,119/255,180/255])
            hold on
        end
        for ixp = 1:3
            plot([VesselPose1(ixp,1)+XVessel(i) VesselPose1(ixp+1,1)+XVessel(i)],[VesselPose1(ixp,2)+YVessel(i) VesselPose1(ixp+1,2)+YVessel(i)],'Color',[31/255,119/255,180/255],'LineWidth',1.2)
            hold on
        end
    end 
    h=plot(PortPosition(1),PortPosition(2),'s','Color','k','MarkerFaceColor','k','MarkerSize',5);
    hold on
    hXLabel=xlabel('x');
    hYLabel=ylabel('y');
    result = DomiPop(iplot,VarietyNumber*sum(GenotypeLength)+1:end-3);
    catetype=[];
    rou=[];
    cita=[];
    dijige=[];
    for i=1:length(result)
        if mod(i,4)==1
            rou=[rou result(i)];
        end
        if mod(i,4)==2
            cita=[cita result(i)];
        end
        if mod(i,4)==3
            catetype=[catetype result(i)];
        end
        if mod(i,4)==0
            dijige=[dijige result(i)];
        end
    end
    for i=1:length(XVessel)
        x111(i)=XVessel(i)+rou(i)*cos(cita(i));
        y111(i)=YVessel(i)+rou(i)*sin(cita(i));
    end
    colorarray = 'kkkkkkkkkkkkkkkkkkkkkkk';
    for i=1:MaxRobotNo
        cate=find(catetype==i);
        if isempty(cate)
            continue
        end
        [vv,v]=sort(dijige(cate));
        po=[PortPosition;[x111(cate(v))' y111(cate(v))'];PortPosition];
        for j = 1:size(po',2)-1
            arrow([po(j,1),po(j,2)],[po(j+1,1),po(j+1,2)],'Color','k','Type','Line','Width',2);
            hold on
        end
        plot(x111(cate(v))',y111(cate(v))','.','Color',[255/255,127/255,14/255],'MarkerSize',10);
    end
    set(h,'LineWidth',1.5)
    set(gca, 'Box', 'on', ...        
        'LineWidth', 1,...                               
        'XGrid', 'off', 'YGrid', 'off', ...                
        'TickDir', 'in', 'TickLength', [.015 .015], ...      
        'XMinorTick', 'off', 'YMinorTick', 'off', ...       
        'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1])    
    set(gca, 'FontName', 'Helvetica')
    set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
    set(gca, 'FontSize', 15)
    set([hXLabel, hYLabel], 'FontSize', 16)
end
Epa = DomiPop(:,(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
f1 = Epa(:,1);
f2 = Epa(:,2);
f3 = Epa(:,3);
plot3(f1,f2,f3,'ro','MarkerFaceColor','r'); grid on;
legend(['All solution count ',num2str(size(AllPop,1))],['Optimal solution count ',num2str(size(DomiPop,1))])
figure(21)
grid on
subplot(2,1,1)
hold on
plot([1:MaxGeneration],reflag1,'d-','Color',[0.1,0.1,0.1],'MarkerFaceColor',[0.1,0.1,0.1],'MarkerSize',3)
hold on
plot([1:MaxGeneration],reflag2,'d-','Color',[0.1,0.7,0.2],'MarkerFaceColor',[0.1,0.7,0.2],'MarkerSize',3)
hold on
plot([1:MaxGeneration],reflag3,'d-','Color',[0.8,0.1,0.6],'MarkerFaceColor',[0.8,0.1,0.6],'MarkerSize',3)
hold on
legend(['No mutation ',num2str(sum(reflag1))],['Random multi-point mutation ',num2str(sum(reflag2))],['Mutation process regulation ',num2str(sum(reflag3))])
xlabel('Number of Iterations')
ylabel('Number of Redundant Patrol')
subplot(2,1,2)
grid on
plot([1:MaxGeneration],crflag1,'d-','Color',[0.1,0.1,0.1],'MarkerFaceColor',[0.1,0.1,0.1],'MarkerSize',3)
hold on
plot([1:MaxGeneration],crflag2,'d-','Color',[0.1,0.7,0.2],'MarkerFaceColor',[0.1,0.7,0.2],'MarkerSize',3)
hold on
plot([1:MaxGeneration],crflag3,'d-','Color',[0.8,0.1,0.6],'MarkerFaceColor',[0.8,0.1,0.6],'MarkerSize',3)
legend(['No mutation ',num2str(sum(crflag1))],['Random multi-point mutation ',num2str(sum(crflag2))],['Mutation process regulation ',num2str(sum(crflag3))])
xlabel('Number of Iteration')
ylabel('Number of Overlappint Patrol')
figure(20)
subplot(2,1,1)
grid on
hold on
plot([1:MaxGeneration],reflag1,'d-','MarkerFaceColor','default','MarkerSize',3)
hold on
plot([1:MaxGeneration],reflag2,'rd-','MarkerFaceColor','default','MarkerSize',3)
hold on
plot([1:MaxGeneration],reflag3,'kd-','MarkerFaceColor','default','MarkerSize',3)
hold on
legend(['No mutation'],['Random multi-point mutation'],['Mutation process regulation'])
xlabel('Number of Iterations')
ylabel('Number of Redundant Patrol')
title('Redundant Patrol')
subplot(2,1,2)
grid on
hold on
plot([1:MaxGeneration],crflag1,'d-','MarkerSize',3)
hold on
plot([1:MaxGeneration],crflag2,'rd-','MarkerSize',3)
hold on
plot([1:MaxGeneration],crflag3,'kd-','MarkerSize',3)
legend(['No mutation'],['Random multi-point mutation'],['Mutation process regulation'])
xlabel('Number of Iteration')
ylabel('Number of Overlappint Patrol')
title('Overlapping Patrol')
allreflag1 = [allreflag1; reflag1];
allreflag2 = [allreflag2; reflag2];
allreflag3 = [allreflag3; reflag3];
allcrflag1 = [allcrflag1; crflag1];
allcrflag2 = [allcrflag2; crflag2];
allcrflag3 = [allcrflag3; crflag3];