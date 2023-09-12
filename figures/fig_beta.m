clear, clc, close all
addpath ../src_figures
%%% Plotting defaults %%%
[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

%% convergence rate for different beta
ref = 3.183099555546359;
% test for 
beta = [1 2 3 4 5];
ilevel = 3:8;
% dimension = 2;
value = zeros(size(beta,2), size(ilevel,2));

value(1,:) = [2.948949116402245   3.187983360746620   3.188038515370521   3.183899095187627   3.183145582488622   3.183109660188394];
value(2,:) = [2.945697165806303   3.225479793839332   3.198047431189706   3.185679524722279   3.183324528947399   3.183102410912782];
value(3,:) = [2.987030335715596   3.278918299953618   3.222112405537001   3.191811872733767   3.184141886962451   3.183190855947344];
value(4,:) = [3.033606402672489   3.347410032525615   3.256085215324941   3.202486278047244   3.185770267935921   3.183366318089950];
value(5,:) = [3.116421579742593   3.423194946244803   3.297732235197555   3.219526457341364   3.188623349456059   3.183647803479647];
err = abs(value - ref)./ref;

h1 = figure(1);
for i=1:5
    hold on; plot(ilevel(1:end), err(i,1:end),'DisplayName',['\beta = ' num2str(i)],...
        'Color',colors{i}, 'Marker' ,markers{i},'markersize',ms,'MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i},'LineWidth',lw);
end
set(gca, 'YScale', 'log');
leg1 = legend('$\beta = 1$','$\beta = 2$','$\beta = 3$','$\beta = 4$','$\beta = 5$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',20);

xla = xlabel('interpolation level $L_I$');
set(xla,'Interpreter','latex');
set(xla,'FontSize',20);
ylabel('relative error', 'FontSize',20);

box on;
ax = gca;  ax.FontSize = 20; 
set(h1, 'Position', [100 100 700 650]);


