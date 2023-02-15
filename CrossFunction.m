function [child1,child2] = CrossFunction(VarietyNumber,GenotypeLength,parent1,parent2,Pc,iter,Gen)
Pc = (Gen-iter)/Gen*Pc;
if(rand<Pc)
    cpoint1 = max(randi(VarietyNumber*sum(GenotypeLength),1,1)-1,1);
    cpoint2 = max(randi(VarietyNumber*sum(GenotypeLength),1,1)-1,1);
    while cpoint1 == cpoint2
        cpoint2 = max(randi(VarietyNumber*sum(GenotypeLength),1,1)-1,1);
    end
    if cpoint1>cpoint2
        cp = cpoint2;
        cpoint2 = cpoint1;
        cpoint1 = cp;
    end
    child1=[parent1(1:cpoint1) parent2(cpoint1+1:cpoint2) parent1(cpoint2+1:VarietyNumber*sum(GenotypeLength))];
    child2=[parent2(1:cpoint1) parent1(cpoint1+1:cpoint2) parent2(cpoint2+1:VarietyNumber*sum(GenotypeLength))];
else
    child1=parent1;
    child2=parent2;
end