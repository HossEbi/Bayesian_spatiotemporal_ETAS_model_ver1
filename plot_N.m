%   Ploting forecasted Number of events for the whole aftershock zone
%   comparison with the observed number of events with M>=m

%% Plot Errorbar

indx_m = find(vec_m==Ml_cat); %  vec_m(indx_m) == ML

gcf=figure;

errorx = 0.011;
errory = 0;

errorL = N_robust_p50(indx_m)-N_robust_p16(indx_m);
errorU = N_robust_p84(indx_m)-N_robust_p50(indx_m);
errorbar(1,N_robust_p50(indx_m),errorL,errorU,'--b','linewidth',2,'CapSize',28);
hold on
errorL = N_robust_p50(indx_m)-N_robust_p02(indx_m);
errorU = N_robust_p98(indx_m)-N_robust_p50(indx_m);
e = errorbar(1,N_robust_p50(indx_m),errorL,errorU,'-ks','CapSize',30,'MarkerEdgeColor','k','MarkerFaceColor',[0.85,0.85,0.85],'MarkerSize',35,'linewidth',3);
hold on
plot(1,N_observed(indx_m),'*','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12,'linewidth',2.5);

text(1,N_robust_p50(indx_m),num2str(round(N_robust_p50(indx_m))),'FontName','Times New Roman','FontSize',18,'FontWeight','normal','FontAngle','Normal','Color','k','HorizontalAlignment','center')
text(1+errorx,N_robust_p16(indx_m),num2str(round(N_robust_p16(indx_m))),'FontName','Times New Roman','FontSize',18,'FontWeight','normal','FontAngle','Normal','Color','b')
text(1+errorx,N_robust_p84(indx_m),num2str(round(N_robust_p84(indx_m))),'FontName','Times New Roman','FontSize',18,'FontWeight','normal','FontAngle','Normal','Color','b')
text(1+errorx,N_robust_p02(indx_m)-errory,num2str(round(N_robust_p02(indx_m))),'FontName','Times New Roman','FontSize',18,'FontWeight','normal','FontAngle','Normal','Color','k')
text(1+errorx,N_robust_p98(indx_m)+errory,num2str(round(N_robust_p98(indx_m))),'FontName','Times New Roman','FontSize',18,'FontWeight','normal','FontAngle','Normal','Color','k')

text(1-errorx,N_observed(indx_m),num2str(round(N_observed(indx_m))),'FontName','Times New Roman','FontSize',18,'FontWeight','normal','FontAngle','Normal','Color','r')

limy = ylim;

xlim([0.98 1.02])
ylim([0 limy(2)+N_robust_p02(indx_m)])
% ylabel(['\bfN (m \geq ',num2str(vec_m(indx_m),'%2.1f'),')'],'fontsize',18)
set(gca,'XTick',[],'YTick',[]);
% set(gca,'visible','off')

limx = xlim;
limy = ylim;
text(limx(1)+errorx/2,limy(2)/2,['\bfN (m \geq ',num2str(vec_m(indx_m),'%2.1f'),')'],'FontName','Times New Roman','FontSize',18,'FontWeight','normal','FontAngle','Normal','Color','k','HorizontalAlignment','center','Rotation',90)

set(gcf,'Position',[680,625,210,350])

saveas(gcf,[output_Dir,'/','ErrorBar.tiff']) 





