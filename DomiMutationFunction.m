function child=DomiMutationFunction(VarietyNumber,GenotypeLength,Pm,parent)
if(rand<Pm)
    position = [];
    for i = 1:VarietyNumber
        position = [position (i-1)*sum(GenotypeLength)+3,(i-1)*sum(GenotypeLength)+4 (i-1)*sum(GenotypeLength)+7 (i-1)*sum(GenotypeLength)+8];
    end 
    index = randperm(length(position));
    mpoint = position(index(1));
    child = parent;
    child(mpoint) = 1-parent(mpoint);
else
    child=[];
end
s = setdiff(1:2*VarietyNumber,i);
t = randperm(length(s));
i = s(t(1));
if(rand<Pm)
    if floor((i-1)/6)==0
        position = [(mod(i,7)-1)*sum(GenotypeLength)+1,(mod(i,7)-1)*sum(GenotypeLength)+2,(mod(i,7)-1)*sum(GenotypeLength)+5,(mod(i,7)-1)*sum(GenotypeLength)+6]; 
    else
        position = [(mod(i-6,7)-1)*sum(GenotypeLength)+1,(mod(i-6,7)-1)*sum(GenotypeLength)+2,(mod(i-6,7)-1)*sum(GenotypeLength)+5,(mod(i-6,7)-1)*sum(GenotypeLength)+6]; 
    end
    index = randi(length(position),1,1);
    mpoint = position(index);
    child = parent;
    child(mpoint) = 1-parent(mpoint);
else
    child=parent;
end