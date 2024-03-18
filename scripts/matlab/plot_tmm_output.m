%%=========================================================================
%   plot_tmm_output
%%-------------------------------------------------------------------------
%   purpose: to plot tmm simulation output
%   author: perrin w. davidson
%   contact: perrinwdavidson@gmail.com
%   date: 01.08.22
%%=========================================================================
%%  configure
%   environment ::
clear; 
close;
clc;

%   plotting ::
set(0, 'defaultAxesFontName', 'CMU Sans Serif', 'defaultAxesFontSize', 12); 

%%  load tmm data
%   get filenames ::
fnames = dir('data/sims/tmm');
fnames = fnames(~[fnames.isdir], :);

%   load ::
for iVar = 1 : 1 : length(fnames)

	if ~strcmp(fnames(iVar).name, '.DS_Store')
	
		load(fnames(iVar).name); 

	end

end

%%  load cmip5 wind factor data
load(fullfile('data', 'sims', 'cmip5', 'wind_factor', 'wind_factor'), 'wind_factor'); 
cmip5_wind_factor = wind_factor; 
clear('wind_factor'); 

%%  single gas anamolies
%   start tiled layout ::
figure;
tplot = tiledlayout(2, 6, 'tileSpacing', 'compact'); 

%   set filenames ::
wind_factors = 0.5 : 0.1 : 1.5; 
NUMFACTS = length(wind_factors);

%   set conversion factor ::
percent = 100; 

%   set pic data ::
pic_dat = pic_output_wind_factor_1x;

%   loop through all values and plot ::
for iFact = 1 : 1 : NUMFACTS

	%    get filename ::
	wind_factor = erase(num2str(wind_factors(iFact)), '.');
	if strcmp(wind_factor, '1.0')

		wind_factor = '1';

	end
	filename = append('lgm_output_wind_factor_', wind_factor, 'x'); 

	%    set data ::
	eval(append('lgm_dat = ', filename)); 

	%    calculate sensitivity ::
	%%%  lgm ::
	lgm_kr_data_plot = lgm_dat.kr_anom_profile(:, 1); 
	lgm_kr_z = lgm_dat.kr_anom_profile(:, 2); 
	lgm_xe_data_plot = lgm_dat.xe_anom_profile(:, 1); 
	lgm_xe_z = lgm_dat.xe_anom_profile(:, 2); 

	%%%  pic ::
	pic_kr_data_plot = pic_dat.kr_anom_profile(:, 1); 
	pic_kr_z = pic_dat.kr_anom_profile(:, 2); 
	pic_xe_data_plot = pic_dat.xe_anom_profile(:, 1); 
	pic_xe_z = pic_dat.xe_anom_profile(:, 2); 

	%    plot ::
	nexttile(); 
	hold('on');
	plot([0 0], [0 5000], '--k'); 
	plot(pic_kr_data_plot .* 100, pic_kr_z, 'r', 'linewidth', 1.5);
	plot(lgm_kr_data_plot .* 100, lgm_kr_z, '--r', 'linewidth', 1.5);
	plot(pic_xe_data_plot .* 100, pic_xe_z, 'b', 'linewidth', 1.5);
	plot(lgm_xe_data_plot .* 100, lgm_xe_z, '--b', 'linewidth', 1.5);
	hold('off');
	title(append(num2str(wind_factors(iFact)), 'x Wind Factor'), 'fontName', 'CMU Sans Serif', 'fontSize', 14)
	set(gca, 'ydir', 'reverse', 'ylim', [0 5000], 'xlim', [-11, 5]);
	if (iFact == 1) || (iFact == 7)

		set(gca, 'ydir', 'reverse', 'ylim', [0 5000], 'xlim', [-11, 5]);

	else

		set(gca, 'ytick', [], 'ydir', 'reverse', 'ylim', [0 5000], 'xlim', [-11, 5]);
	end

end

%   make look nice ::
set(gcf, 'position', [0 0 2400 1200]); 
xlabel(tplot, 'Solubility anamoly, \Delta (%)', 'fontName', 'CMU Sans Serif', 'fontSize', 16); 
ylabel(tplot, 'Depth (m)', 'fontName', 'CMU Sans Serif', 'fontSize', 16); 
leg = legend('Equilibrium', 'PIC Kr', 'LGM Kr', 'PIC Xe', 'LGM Xe', 'fontName', 'CMU Sans Serif'); 
leg.Layout.Tile = 'east'; 

%   export graphics ::
exportgraphics(tplot, fullfile('plots', 'lgm_wind_perturbation', 'lgm_wind_perturbation_single_profiles.png'), 'resolution', 300); 

%%  plot depth profiles
%   start tiled layout ::
figure;
tplot = tiledlayout(2, 6, 'tileSpacing', 'compact'); 

%   set filenames ::
wind_factors = 0.5 : 0.1 : 1.5; 
NUMFACTS = length(wind_factors); 

%   set pic data ::
pic_dat = pic_output_wind_factor_1x;

%   loop through all values and plot ::
for iFact = 1 : 1 : NUMFACTS

	%    get filename ::
	wind_factor = erase(num2str(wind_factors(iFact)), '.');
	if strcmp(wind_factor, '1.0')

		wind_factor = '1';

	end
	filename = append('lgm_output_wind_factor_', wind_factor, 'x'); 

	%    set data ::
	eval(append('lgm_dat = ', filename)); 

	%    calculate sensitivity ::
	sensitivity(iFact, :) = [pic_dat.kr_n2_anom_global - lgm_dat.kr_n2_anom_global, ...
	             	         pic_dat.xe_n2_anom_global - lgm_dat.xe_n2_anom_global]; 

	%    plot ::
	nexttile(); 
	hold('on');
	plot([0 0], [0 5000], '--k'); 
	plot(pic_dat.kr_n2_anom_profile(:, 1) .* 100, pic_dat.kr_n2_anom_profile(:, 2), 'r', 'linewidth', 1.5);
	plot(lgm_dat.kr_n2_anom_profile(:, 1) .* 100, lgm_dat.kr_n2_anom_profile(:, 2), '--r', 'linewidth', 1.5);
	plot(pic_dat.xe_n2_anom_profile(:, 1) .* 100, pic_dat.xe_n2_anom_profile(:, 2), 'b', 'linewidth', 1.5);
	plot(lgm_dat.xe_n2_anom_profile(:, 1) .* 100, lgm_dat.xe_n2_anom_profile(:, 2), '--b', 'linewidth', 1.5);
	hold('off');
	title(append(num2str(wind_factors(iFact)), 'x Wind Factor'), 'fontName', 'CMU Sans Serif', 'fontSize', 14)
	if (iFact == 1) || (iFact == 7)

		set(gca, 'ydir', 'reverse', 'ylim', [0 5000], 'xlim', [-5, 0.5]);

	else

		set(gca, 'ytick', [], 'ydir', 'reverse', 'ylim', [0 5000], 'xlim', [-5, 0.5]);
	end

end

%   make look nice ::
set(gcf, 'position', [0 0 2400 1200]); 
xlabel(tplot, 'Solubility anamoly, \Delta (%)', 'fontName', 'CMU Sans Serif', 'fontSize', 16); 
ylabel(tplot, 'Depth (m)', 'fontName', 'CMU Sans Serif', 'fontSize', 16); 
leg = legend('Equilibrium', 'PIC Kr/N_2', 'LGM Kr/N_2', 'PIC Xe/N_2', 'LGM Xe/N_2', 'fontName', 'CMU Sans Serif'); 
leg.Layout.Tile = 'east'; 

%   export graphics ::
exportgraphics(tplot, fullfile('plots', 'lgm_wind_perturbation', 'lgm_wind_perturbation_ratio_profiles.png'), 'resolution', 300); 

%%  delta delta
%   change wind factors to percent change ::
wind_factors = (wind_factors - 1) * 100; 

%   fit function basis independent matrix ::
E = [ones(NUMFACTS, 1), wind_factors', wind_factors.^2'];

%   kr:n2 fit ::
x_kr_tilde = inv(E' * E) * E' * sensitivity(:, 1); 
n_kr_tilde = sensitivity(:, 1) - (E * x_kr_tilde); 
P_kr = sqrt(var(n_kr_tilde) .* inv(E' * E)); 
x_min_kr = -x_kr_tilde(2) / (2 .* x_kr_tilde(3));

%   xe:n2 fit ::
x_xe_tilde = inv(E' * E) * E' * sensitivity(:, 2); 
n_xe_tilde = sensitivity(:, 2) - (E * x_xe_tilde); 
P_xe = sqrt(var(n_xe_tilde) .* inv(E' * E)); 
x_min_xe = -x_xe_tilde(2) / (2 .* x_xe_tilde(3));

%   set plotting values for fit ::
x_mod_fit = linspace(0.4, 1.6, 1000)
x_mod_fit = (x_mod_fit - 1) * 100; 

%   make function ::
fit_fun = @(x, coeffs) coeffs(1) + (coeffs(2) .* x) + (coeffs(3) .* (x .^ 2)); 
err_fun = @(x, err) sqrt((err(1) ^ 2) + ((err(2) .* x) .^ 2) + ((err(3) .* (x .^ 2)) .^ 2));

%   start plot ::
figure;
tplot = tiledlayout(3, 1, 'tileSpacing', 'compact'); 

%   set plotting options ::
plot_err = 0;  % 0 (no error bar ploting sqrt{C_rr}), 1 (plotting of regression error)

%  kr:n2 ::
%%% fill between ::
if plot_err

	x_err_plot = [x_mod_fit, fliplr(x_mod_fit)];
	error_plot = [fit_fun(x_mod_fit, x_kr_tilde) - err_fun(x_mod_fit, diag(P_kr)), ...
	      	      fliplr(fit_fun(x_mod_fit, x_kr_tilde) + err_fun(x_mod_fit, diag(P_kr)))];

end

%%% plot ::
nexttile([2 1]);
hold('on');
plot(x_mod_fit, -fit_fun(x_mod_fit, x_kr_tilde).*100, '--r', 'lineWidth', 2); 
if plot_err

	fill(x_err_plot, error_plot, 'k', 'faceAlpha', 0.3, 'edgeColor', 'none'); 
	plot(x_mod_fit, fit_fun(x_mod_fit, x_kr_tilde) + err_fun(x_mod_fit, diag(P_kr)), '--k', 'lineWidth', 1); 
	plot(x_mod_fit, fit_fun(x_mod_fit, x_kr_tilde) - err_fun(x_mod_fit, diag(P_kr)), '--k', 'lineWidth', 1);

end
scatter(wind_factors, -sensitivity(:, 1).*100, 300, 'd', 'filled', 'r')

%  xe:n2 ::
%%% fill between ::
if plot_err
	
	x_err_plot = [x_mod_fit, fliplr(x_mod_fit)];
	error_plot = [fit_fun(x_mod_fit, x_xe_tilde) - err_fun(x_mod_fit, diag(P_xe)), ...
	      	      fliplr(fit_fun(x_mod_fit, x_xe_tilde) + err_fun(x_mod_fit, diag(P_xe)))];
end

%%% plot ::
plot(x_mod_fit, -fit_fun(x_mod_fit, x_xe_tilde).*100, '--b', 'lineWidth', 2); 
if plot_err

	fill(x_err_plot, error_plot, 'k', 'faceAlpha', 0.3, 'edgeColor', 'none'); 
	plot(x_mod_fit, fit_fun(x_mod_fit, x_xe_tilde) + err_fun(x_mod_fit, diag(P_xe)), '--k', 'lineWidth', 1); 
	plot(x_mod_fit, fit_fun(x_mod_fit, x_xe_tilde) - err_fun(x_mod_fit, diag(P_xe)), '--k', 'lineWidth', 1);

end
scatter(wind_factors, -sensitivity(:, 2).*100, 300, 'd', 'filled', 'b')
hold('off');
set(gca, 'xlim', ([0.5, 1.5]-1)*100, 'box', 'on', 'xTickLabels', {}); 

%   make look nice ::
ylabel('\Delta_{LGM} - \Delta_{PIC} (%)', 'fontName', 'CMU Sans Serif', 'fontSize', 16); 
legend('', 'Kr/N_2', '', 'Xe/N_2', 'fontName', 'CMU Sans Serif', 'fontSize', 16); 

%   plot box plot ::
nexttile(); 
boxplot(([cmip5_wind_factor.south(1:end-1); cmip5_wind_factor.north(1:end-1)]'-1).*100, 'orientation', 'horizontal'); 
ylabel('CMIP5 Data Products'); 

%   make look nice ::
xlabel(tplot, 'Percent Wind Change in High Latitudes (> 50\circ Lat)', 'fontName', 'CMU Sans Serif', 'fontSize', 16); 
set(gca, 'yTickLabels', {'South', 'North'}, 'xLim', [-50, 50]); 
set(gcf, 'position', [0 0 600 900]); 

%   export graphics ::
exportgraphics(tplot, fullfile('plots', 'lgm_wind_perturbation', 'lgm_wind_perturbation_sensitivity.png'), 'resolution', 300); 

%%  end program
