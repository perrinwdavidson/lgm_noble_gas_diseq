model_data_lon = readmatrix("data/sims/cmip5/sic/sic_lon.csv");
model_data_lat = readmatrix("data/sims/cmip5/sic/sic_lat.csv");
model_data_int = readmatrix("data/sims/cmip5/sic/sic_int.csv");

idx_cut = 4;

list_factory = fieldnames(get(groot, "factory"));
index_interpreter = find(contains(list_factory, "Interpreter"));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)}, "factory", "default");
    set(groot, default_name, "latex");
end

int_plt = figure;
tiledlayout(1, 2);
nexttile(); 
scatter(model_lon(:), model_lat(:), 40, model_data_mean(:), "filled");
colorbar();
colormap(winter(9));
title("Raw", "FontSize", 16);
set(gca, "Box", "on", "LineWidth", 0.5, "xlim", [0, 360], "ylim", [-80, 80]);
nexttile(); 
scatter(model_data_lon(:), model_data_lat(:), 40, model_data_int(:), "filled"); 
colorbar();
colormap(winter(9));
title("Interpolated", "FontSize", 16)
set(gca, "Box", "on", "LineWidth", 0.5, "xlim", [0, 360], "ylim", [-80, 80]);
set(gcf, "Position", [10, 10, 1000, 250]);
exportgraphics(int_plt, "plots/ice/pic/int_comp.png", "Resolution", 300);

int_plt_close = figure;
tiledlayout(1, 2);
nexttile(); 
scatter(model_lon(:), model_lat(:), 50, model_data_mean(:), "filled");
colorbar();
colormap(winter(9));
title("Raw", "FontSize", 16);
set(gca, "Box", "on", "LineWidth", 0.5, "xlim", [150, 250], "ylim", [-80, -60]);
nexttile(); 
scatter(model_data_lon(:), model_data_lat(:), 100, model_data_int(:), "filled"); 
colorbar();
colormap(winter(9));
title("Interpolated", "FontSize", 16)
set(gca, "Box", "on", "LineWidth", 0.5, "xlim", [150, 250], "ylim", [-80, -60]);
set(gcf, "Position", [10, 10, 1000, 100]);
exportgraphics(int_plt_close, "plots/ice/pic/int_comp_close.png", "Resolution", 300)