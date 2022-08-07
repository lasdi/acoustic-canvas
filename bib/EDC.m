
function [edc,t60, edt] = EDC(RIR, fs, flag)

th = 60; 
decai=pow2db(RIR.^2)-max(pow2db(RIR.^2));

th_RT = -(th)-5; 

ind = find(decai>th_RT);
td = ind(end);

edc(td:-1:1)=10*log10(cumsum(RIR(td:-1:1).^2)/sum(RIR(1:td).^2));
edc = [edc(:); ones(length(decai)-length(edc),1)*edc(end)]; edc=edc(:);

ind1 = find(edc>-5);    ind1=ind1(end);
ind2 = find(edc>-35);   ind2=ind2(end);
ind3 = find(edc>-15);   ind3=ind3(end);

t = (0:length(RIR)-1)/fs; t=t(:);

p = polyfit(t(ind1:ind2),edc(ind1:ind2),1);
y   = p(1) * (t(ind1:ind2) - t(ind1)) + p(2);
b = y(1)-(p(1)* t(ind1));

pp = polyfit(t(ind1:ind3),edc(ind1:ind3),1);
yy   = pp(1) * (t(ind1:ind3) - t(ind1)) + pp(2);
bb = yy(1)-(pp(1)* t(ind1));

t60 = (-60- b)/p(1);
edt = (-60- bb)/pp(1);

    if flag 
        figure
        plot(t*1000,decai);
        hold on
        plot(t*1000,edc,'LineWidth',2);
        plot(t(ind1:ind2)*1000,p(1)*t(ind1:ind2)+p(2),'LineWidth',3);
        plot(t(ind1:ind2)*1000,pp(1)*t(ind1:ind2)+pp(2),'LineWidth',3);
        plot(t*1000,-5*ones(size(edc)),'y--');
        plot(t*1000,-35*ones(size(edc)),'y--');
        xlabel('tempo (ms)');
        ylabel('Energia normalizada (dB)');
        legend('h^2','EDC','Ajuste linear [RT30]','EDT',['Intervalo RT30: -5dB a -35dB'], 'Location', 'South');
        text(((ind2+15)/fs)*1000,-30,['T60 = ' num2str(t60*1000) 'ms'])
        text(((ind2+15)/fs)*1000,-10,['EDT = ' num2str(edt*1000) 'ms'])
        axis([0 t(end)*1000 -inf 0])
    end
end

















