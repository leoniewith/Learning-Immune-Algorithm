function [flag,positionindex] = trytest2(x,GenotypeLength,VarietyNumber,XVessel,YVessel,RVessel,MaxRobotNo,PortPosition)
global XWind YWind RWind
positionindex = [];
flag = [0 0];
xnew = reshape(x,length(GenotypeLength),VarietyNumber)';
rou = xnew(:,1);
cita = xnew(:,2);
for i=1:length(XVessel)
    x111(i)=XVessel(i)+rou(i)*cos(cita(i));
    y111(i)=YVessel(i)+rou(i)*sin(cita(i));
end
arrayvessel = [x111' y111' xnew(:,3:4)];
for i = 1:MaxRobotNo
    index = find(arrayvessel(:,3)==i);
    if isempty(index)
        coortemp{i} = [];
        continue
    end
    arraytemp = arrayvessel(index,:);
    [~,v] = sort(arraytemp(:,4));
    arraytemp = arraytemp(v,:);
    coortemp{i} = [PortPosition;arraytemp(:,1) arraytemp(:,2);PortPosition];
end
for i = 1:MaxRobotNo
    index = find(arrayvessel(:,3)==i);
    if isempty(index)
        isuseless(i) = 0;
        continue
    end
    co1 = [XVessel(index)' YVessel(index)'];
    dp = RVessel(index);
    d = [];
    for j = 1:size(co1,1)
        dtemp = [];
        for k = 1:MaxRobotNo
            if i ~= k
                kindex = find(arrayvessel(:,3)==k);
                if isempty(kindex)
                    dtemp = ones(1,length(index))*1000;
                else
                    co2 = coortemp{k};
                    dtemp = [dtemp p_poly_dist(co1(j,1),co1(j,2),co2(:,1),co2(:,2))];
                end
            end
        end
        d = [d min(dtemp)];
    end
    sindex = d < dp;  
    if sum(sindex)==size(co1,1)
        isuseless(i) = 1;
        positionindex = [positionindex;index];
    else
        isuseless(i) = 0;
    end
end
if sum(isuseless)>0 
    flag(1) = flag(1)+1;
end
if length(coortemp) < 2 
    flag(2) = flag(2)+0;
else
    for i = 1:length(coortemp)-1
        if isempty(coortemp{i})
            continue
        end
        for j = i+1:length(coortemp)
            if isempty(coortemp{j})
                continue
            end
            co1 = coortemp{i};
            co2 = coortemp{j};
            [x0,y0] = polyxpoly(co1(:,1),co1(:,2),co2(:,1),co2(:,2));
            for k = length(x0):-1:1
                if x0(k) == PortPosition(1) && y0(k) == PortPosition(2)
                    x0(k) = [];
                    y0(k) = [];
                end
            end
            if isempty(x0)
                flag(2) = flag(2)+0;
            else
                flag(2) = flag(2)+1;
            end
        end
    end
end