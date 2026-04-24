
function b=addLinTrend(x,y,c)

[b,~,stats]=glmfit(x,y,'normal');
x0 = linspace(min(x),max(x),256)';
[yhat,dlo,dhi]=glmval(b,x0,'identity',stats);
hold on
shadedErrorBar(x0,yhat,[dlo dhi]',c);
hold off