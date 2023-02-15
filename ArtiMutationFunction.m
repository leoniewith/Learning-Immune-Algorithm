function child=ArtiMutationFunction(i,VarietyNumber,GenotypeLength,Pm,parent)
if(rand<Pm)
    if floor((i-1)/VarietyNumber)==0
        position = [(mod(i,VarietyNumber+1)-1)*sum(GenotypeLength)+GenotypeLength(1)-1,(mod(i,VarietyNumber+1)-1)*sum(GenotypeLength)+GenotypeLength(1),(mod(i,VarietyNumber+1)-1)*sum(GenotypeLength)+GenotypeLength(1)+GenotypeLength(2)-1,(mod(i,VarietyNumber+1)-1)*sum(GenotypeLength)+GenotypeLength(1)+GenotypeLength(2)]; 
    else
        position = [(mod(i-VarietyNumber,VarietyNumber+1)-1)*sum(GenotypeLength)+GenotypeLength(1)-1,(mod(i-VarietyNumber,VarietyNumber+1)-1)*sum(GenotypeLength)+GenotypeLength(1),(mod(i-VarietyNumber,VarietyNumber+1)-1)*sum(GenotypeLength)+GenotypeLength(1)+GenotypeLength(2)-1,(mod(i-VarietyNumber,VarietyNumber+1)-1)*sum(GenotypeLength)+GenotypeLength(1)+GenotypeLength(2)]; 
    end
    index = randi(length(position),1,1);
    mpoint = position(index);
    child = parent;
    child(mpoint) = 1-parent(mpoint);
else
    child=parent;
end
s = setdiff(1:2*VarietyNumber,i);
t = randperm(length(s));
i = s(t(1));
if(rand<Pm)
    if floor((i-1)/6)==0
        position = [(mod(i,VarietyNumber+1)-1)*sum(GenotypeLength)+1,(mod(i,VarietyNumber+1)-1)*sum(GenotypeLength)+2,(mod(i,VarietyNumber+1)-1)*sum(GenotypeLength)+GenotypeLength(1)+1,(mod(i,VarietyNumber+1)-1)*sum(GenotypeLength)+GenotypeLength(1)+2]; 
    else
        position = [(mod(i-VarietyNumber,VarietyNumber+1)-1)*sum(GenotypeLength)+1,(mod(i-VarietyNumber,VarietyNumber+1)-1)*sum(GenotypeLength)+2,(mod(i-VarietyNumber,VarietyNumber+1)-1)*sum(GenotypeLength)+GenotypeLength(1)+1,(mod(i-VarietyNumber,VarietyNumber+1)-1)*sum(GenotypeLength)+GenotypeLength(1)+2]; 
    end
    index = randi(length(position),1,1);
    mpoint = position(index);
    child = parent;
    child(mpoint) = 1-parent(mpoint);
else
    child=parent;
end