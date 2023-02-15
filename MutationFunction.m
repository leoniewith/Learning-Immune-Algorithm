function child=MutationFunction(MaxRobotNo,RVessel,VarietyNumber,GenotypeLength,Pm,parent)
CateNumber = [];
for j = 1:VarietyNumber
    for k = 1:length(GenotypeLength)
        if k == 1
            startpoint = 1;
            endpoint = startpoint + GenotypeLength(k) - 1;
        else
            startpoint = endpoint + 1;
            endpoint = startpoint + GenotypeLength(k) - 1;
        end
        varietycode = [];
        varietycode = parent((j-1)*(sum(GenotypeLength))+startpoint:(j-1)*(sum(GenotypeLength))+endpoint);
        if mod(k,4) == 1
            parent(VarietyNumber*sum(GenotypeLength)+length(GenotypeLength)*(j-1)+k)=(2.^(GenotypeLength(k)-1:-1:0)*...
                varietycode')'/(2.^GenotypeLength(k)-1)*RVessel(j);
        end
        if mod(k,4) == 2
            parent(VarietyNumber*sum(GenotypeLength)+length(GenotypeLength)*(j-1)+k)=(2.^(GenotypeLength(k)-1:-1:0)*...
                varietycode')'/(2.^GenotypeLength(k)-1)*2*pi;
        end
        if mod(k,4) == 3
            parent(VarietyNumber*sum(GenotypeLength)+length(GenotypeLength)*(j-1)+k)=min(floor((2.^(GenotypeLength(k)-1:-1:0)*...
                varietycode')'/(2.^GenotypeLength(k)-1)*MaxRobotNo)+1,MaxRobotNo);
            CateNumber = [CateNumber min(floor((2.^(GenotypeLength(k)-1:-1:0)*...
                varietycode')'/(2.^GenotypeLength(k)-1)*MaxRobotNo)+1,MaxRobotNo)];
        end
    end
end
CateType = [];
for s = 1:MaxRobotNo 
    vindex = find(CateNumber==s); 
    if isempty(vindex) 
        continue
    end
    if length(vindex)>0
        CateType = [CateType s];
    end
    RobotNoinCate = length(vindex);
    vesselvalue = [];
    for m = 1:RobotNoinCate
        startpoint = sum(GenotypeLength)-GenotypeLength(4)+1;
        endpoint = sum(GenotypeLength);
        varietycode = [];
        varietycode = parent((vindex(m)-1)*(sum(GenotypeLength))+startpoint:(vindex(m)-1)*(sum(GenotypeLength))+endpoint);
        vesselvalue = [vesselvalue (2.^(GenotypeLength(4)-1:-1:0)*varietycode')'/(2.^GenotypeLength(4)-1)];
    end
    [~,bindex] = sort(vesselvalue);
    vindex = vindex(bindex);
    for m = 1:RobotNoinCate
        parent(VarietyNumber*sum(GenotypeLength)+length(GenotypeLength)*(vindex(m)-1)+4) = m;
    end
end
mutations = [];
if(rand<Pm)
    mpoint = randi(VarietyNumber*sum(GenotypeLength),1,1);
    mutations = [mutations mpoint];
    if mpoint-(ceil(mpoint/sum(GenotypeLength))-1)*sum(GenotypeLength)>=1 && mpoint-(ceil(mpoint/sum(GenotypeLength))-1)*sum(GenotypeLength)<=(GenotypeLength(1)+GenotypeLength(2)) 
        Pm=0.1;
        whichvessel = ceil(mpoint/sum(GenotypeLength));
        whichindex = find(CateNumber==CateNumber(whichvessel)); 
        whichindex = setdiff(whichindex,whichvessel);
        for i = 1:length(whichindex)
            candipoint = (whichindex(i)-1)*sum(GenotypeLength)+1:(whichindex(i)-1)*sum(GenotypeLength)+GenotypeLength(1)+GenotypeLength(2);
            if rand<Pm
                temp = randi(length(candipoint));
                mutations = [mutations candipoint(temp)];
            end
        end   
    end
    if mpoint-(ceil(mpoint/sum(GenotypeLength))-1)*sum(GenotypeLength)>=GenotypeLength(1)+GenotypeLength(2)+1 && mpoint-(ceil(mpoint/sum(GenotypeLength))-1)*sum(GenotypeLength)<=(GenotypeLength(1)+GenotypeLength(2)+GenotypeLength(3))
        Pm=1;
        whichvessel = ceil(mpoint/sum(GenotypeLength));
        candipoint = [];
        for i = 1:VarietyNumber
            if i ~= whichvessel
                candipoint = [candipoint (i-1)*sum(GenotypeLength)+GenotypeLength(1)+GenotypeLength(2)+1:(i-1)*sum(GenotypeLength)+GenotypeLength(1)+GenotypeLength(2)+GenotypeLength(3)]; 
            else
                candipoint = [candipoint (i-1)*sum(GenotypeLength)+1:(i-1)*sum(GenotypeLength)+GenotypeLength(1)+GenotypeLength(2)]; 
            end
        end
        if rand<Pm
            temp = randi(length(candipoint));
            mutations = [mutations candipoint(temp)];
        end
    end
    if mpoint-(ceil(mpoint/sum(GenotypeLength))-1)*sum(GenotypeLength) >= GenotypeLength(1)+GenotypeLength(2)+GenotypeLength(3)+1 
        Pm=1;
        whichvessel = ceil(mpoint/sum(GenotypeLength));
        whichindex = find(CateNumber==CateNumber(whichvessel));
        whichindex = setdiff(whichindex,whichvessel);
        for i = 1:length(whichindex)
            candipoint = (whichindex(i)-1)*sum(GenotypeLength)+GenotypeLength(1)+GenotypeLength(2)+GenotypeLength(3)+1:whichindex(i)*sum(GenotypeLength);
            if rand<Pm
                temp = randi(length(candipoint));
                mutations = [mutations candipoint(temp)];
            end
        end  
    end
    child = parent;
    child(mutations) = 1-parent(mutations);
else
    child=parent;
end
child = child(1:VarietyNumber*sum(GenotypeLength));