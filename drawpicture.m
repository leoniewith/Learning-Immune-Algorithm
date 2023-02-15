figure(12)
Epa = DomiPop(:,(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
XWind = 69.73;
YWind = 58.15;
RWind = 15;
PortPosition = [0 50];
Epa = DomiPop(:,(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+1:(sum(GenotypeLength)+length(GenotypeLength))*length(XVessel)+3);
f1 = Epa(:,1);
f2 = Epa(:,2);
f3 = Epa(:,3);
for iplot = 1:size(Epa,1)
    axis([0 100 0 100])
    hold on
    axis square;
    hold on
    for i = 1:length(XVessel) 
        theta = 0:pi/50:2*pi;
        x0 = XVessel(i)+RVessel(i)*cos(theta);
        y0 = YVessel(i)+RVessel(i)*sin(theta);
        plot(x0,y0,'r:',XVessel(i),YVessel(i),'rs','MarkerFaceColor','r');
    end
    for i = 1:length(XWind) 
        theta = 0:pi/50:2*pi;
        x0 = XWind(i)+RWind(i)*cos(theta);
        y0 = YWind(i)+RWind(i)*sin(theta);
        plot(x0,y0,'k:',XWind(i),YWind(i),'ko','MarkerFaceColor','k');
    end
    plot(PortPosition(1),PortPosition(2),'gs','MarkerFaceColor','g');
    hold on
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
    colorarray = 'rgcmbrgcmbrgcmbrgcmb';
    for i=1:MaxRobotNo
        cate=find(catetype==i);
        if isempty(cate)
            continue
        end
        [vv,v]=sort(dijige(cate));
        po=[PortPosition;[x111(cate(v))' y111(cate(v))'];PortPosition];
        for j = 1:size(po',2)-1
            arrow([po(j,1),po(j,2)],[po(j+1,1),po(j+1,2)],'Color',colorarray(i),'Type','line');
            hold on
        end
        plot(x111(cate(v))',y111(cate(v))','k.');
    end
end