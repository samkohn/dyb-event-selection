%pylab
mpl.rc_file('plot_generators/matplotlibrc')
import sqlite3
conn = sqlite3.Connection('/global/cscratch1/sd/skohn/dyb30/fitter_validation_poisson.db')
cursor = conn.cursor()
cursor.execute('SELECT ChiSquare_Fit FROM fitter_validation_results WHERE Category = "ToyMC 26 all systematics and stat flucts, rel-escale only pulls 4/11/2021" AND IsAvgNear = 1 AND IsRateOnly = 1')
results = array(cursor.fetchall()).flatten()
hist(results, bins=bins, density=True, histtype='step')
bins = np.linspace(0, 12, 120)
hist(results, bins=bins, density=True, histtype='step')
from scipy.stats import chi2
xvals = bins
yvals = chi2.pdf(xvals, df=3)
plot(xvals, yvals)
cursor.execute('SELECT FitSinSqT13 FROM fitter_validation_results WHERE Category = "ToyMC 26 all systematics and stat flucts, all pulls 4/11/2021" AND IsAvgNear = 1 AND IsRateOnly = 1 ORDER BY `Index`')
results_all = np.array(cursor.fetchall()).flatten()
results_all.shape
cursor.execute('SELECT FitSinSqT13 FROM fitter_validation_results WHERE Category = "ToyMC 26 all systematics and stat flucts, all pulls 1/12/2021" AND IsAvgNear = 1 AND IsRateOnly = 1 ORDER BY `Index`')
results_all = np.array(cursor.fetchall()).flatten()
results_all.shape
cursor.execute('SELECT FitSinSqT13 FROM fitter_validation_results WHERE Category = "ToyMC 26 all systematics and stat flucts, no pulls 1/12/2021" AND IsAvgNear = 1 AND IsRateOnly = 1 ORDER BY `Index`')
results_none = np.array(cursor.fetchall()).flatten()
plot(results_none, results_all, '.')
w160-1710
2160-1710
cursor.execute('SELECT FitSinSqT13 FROM fitter_validation_results WHERE Category = "ToyMC 26 all systematics and stat flucts, no pulls 1/12/2021" AND IsAvgNear = 1 AND IsRateOnly = 1 ORDER BY `Index` LIMIT 190')
results_none = np.array(cursor.fetchall()).flatten()
plot(results_none, results_all, '.')
plot(results_none, results_all[:190], '.')
hist(results_all[:190]/results_none - 1)
cursor.execute('SELECT FitSinSqT13 FROM fitter_validation_results WHERE Category = "ToyMC 26 all systematics and stat flucts, no pulls 1/12/2021" AND IsAvgNear = 1 AND IsRateOnly = 1 ORDER BY `Index`')
results_none = np.array(cursor.fetchall()).flatten()
hist(results_all, bins=bins, histtype='step', linecolor='k')
hist(results_all, bins=bins, histtype='step', color='k')
cursor.execute('SELECT ChiSquare_Fit FROM fitter_validation_results WHERE Category = "ToyMC 26 all systematics and stat flucts, no pulls 1/12/2021" AND IsAvgNear = 1 AND IsRateOnly = 1')
results_none = np.array(cursor.fetchall()).flatten()
cursor.execute('SELECT ChiSquare_Fit FROM fitter_validation_results WHERE Category = "ToyMC 26 all systematics and stat flucts, all pulls 1/12/2021" AND IsAvgNear = 1 AND IsRateOnly = 1')
results_all = np.array(cursor.fetchall()).flatten()
results.none.shape
results_none.shape
cursor.execute('SELECT ChiSquare_Fit FROM fitter_validation_results WHERE Category = "ToyMC 26 all systematics and stat flucts, all pulls 1/12/2021" AND IsAvgNear = 1 AND IsRateOnly = 1 AND MOD(`Index`, 240) < 190')
cursor.execute('SELECT PI() FROM fitter_validation_results LIMIT 1')
cursor.execute('SELECT 10 % 7 FROM fitter_validation_results LIMIT 1')
cursor.fetchall()
cursor.execute('SELECT ChiSquare_Fit FROM fitter_validation_results WHERE Category = "ToyMC 26 all systematics and stat flucts, all pulls 1/12/2021" AND IsAvgNear = 1 AND IsRateOnly = 1 AND (`Index` % 240) < 190')
cursor.execute('SELECT ChiSquare_Fit FROM fitter_validation_results WHERE Category = "ToyMC 26 all systematics and stat flucts, no pulls 1/12/2021" AND IsAvgNear = 1 AND IsRateOnly = 1 AND (`Index` % 240) < 190')
results_none = np.array(cursor.fetchall()).flatten()
results_none.shape
cursor.execute('SELECT ChiSquare_Fit FROM fitter_validation_results WHERE Category = "ToyMC 26 all systematics and stat flucts, all pulls 1/12/2021" AND IsAvgNear = 1 AND IsRateOnly = 1')
results_all = np.array(cursor.fetchall()).flatten()
results_all.shape
hist(results_all, bins=bins, histtype='step', color='k')
hist(results_all, bins=bins, histtype='step', color='k')
rcParams['lines.linewidth']
hist(results_all, bins=bins, histtype='step', color='k', linewidth=2)
rcParams['hist.linewidth']
mpl.rc_file('plot_generators/matplotlibrc')
hist(results_all, bins=bins, histtype='step', color='k')
bins = np.linspace(0, 12, 61, endpoint=True)
hist(results_all, bins=bins, histtype='step', color='k')
hist(results_none, bins=bins, histtype='step', color='k', linestyle=':')
hist(results_none, bins=bins, histtype='stepfilled', color='k', facecolor='k', hatch='.')
hist(results_all, bins=bins, histtype='step', color='k')
hist(results_none, bins=bins, histtype='stepfilled', color='k', facecolor='k', hatch='*')
hist(results_all, bins=bins, histtype='step', color='k')
hist(results_none, bins=bins, histtype='stepfilled', color='k', hatch='*')
hist(results_none, bins=bins, histtype='stepfilled', color='k', hatch='/', linewidth=1)
hist(results_none, bins=bins, histtype='stepfilled', hatch='/', linewidth=1)
hist(results_none, bins=bins, histtype='stepfilled', facecolor=None, edgecolor='k', hatch='/', linewidth=1)
hist(results_none, bins=bins, histtype='stepfilled', facecolor=(255, 255, 255, 0), edgecolor='k', hatch='/', linewidth=1)
hist(results_none, bins=bins, histtype='stepfilled', facecolor='w', edgecolor='k', hatch='.', linewidth=1)
hist(results_none, bins=bins, histtype='stepfilled', facecolor='w', edgecolor='k', hatch='/')
hist(results_all, bins=bins, histtype='step', color='k')
hist(results_all, bins=bins, histtype='step')
hist(results_none, bins=bins, histtype='stepfilled', facecolor='w', hatch='/')
hist(results_none, bins=bins, histtype='stepfilled', facecolor='w', edgecolor='auto', hatch='/')
hist(results_none, bins=bins, histtype='stepfilled',hatch='/')
hist(results_none, bins=bins, histtype='step', hatch='/')
hist(results_all, bins=bins, histtype='step')
hist(results_none, bins=bins, histtype='step')
from scipy.stats import chi2
xvals = np.linspace(0, 12, 121, endpoint=True)
yvals = chi2.pdf(xvals, df=3) * results_none.shape * np.diff(results_none)[0]
plot(xvals, yvals, 'k')
results_none.shape
np.diff(results_none)[0]
hist(results_all, bins=bins, histtype='step')
hist(results_none, bins=bins, histtype='step')
yvals = chi2.pdf(xvals, df=3) * results_none.shape * np.diff(bins)[0]
plot(xvals, yvals, 'k')
obj_all = _87[2][0]
obj_none = _88[2][0]
obj_curve = _90[0]
legend([obj_curve, obj_all, obj_none], [r'$\chi^2$ with NDF=3', 'All pull parameters enabled', 'No pull parameters enabled'], frameon=False)
grid(False)
grid(axis='x')
grid(True)
ax2 = gca().twinx()
ax2.set_ylim([0, gca().get_ylim()[1]/results_none.shape/np.diff(bins)[0]])
gca().get_ylim()
gcf()
_101.get_axes()
_101.get_axes()[0]
_101.get_axes()[0].get_ylim()
ax2.set_ylim([0, 150.15/results_none.shape/np.diff(bins)[0]])
ax2.set_ylim([0, 150.15/len(results_none)/np.diff(bins)[0]])
ax2.grid(False)
ax1 = _103
ax2.ylabel(r'$P(\chi^2)$')
ax2.set_ylabel(r'$P(\chi^2)$')
ax1.set_ylabel(r'$N(\chi^2)$')
ax1.set_ylabel(r'$N(\chi^2)$ with 1710 trials')
ax1.set_xlabel(r'Fit $\chi^2_{min}$')
grid(False)
ax1.grid(False)
figure()
obj_all_line = plot([1, 2, 3], [4, 5, 6])
obj_none_line = plot([1, 2, 3], [4, 5, 6])
_101
legend([obj_curve, obj_all_line, obj_none_line], [r'$\chi^2$ with NDF=3', 'All pull parameters enabled', 'No pull parameters enabled'], frameon=False)
obj_all_line
legend([obj_curve, obj_all_line[0], obj_none_line[0]], [r'$\chi^2$ with NDF=3', 'All pull parameters enabled', 'No pull parameters enabled'], frameon=False)
legend([obj_curve, obj_all_line[0], obj_none_line[0]], [r'$\chi^2$ with NDF=3', 'All pull parameters enabled', 'No pull parameters enabled'], frameon=False)
_94.remove()
savefig('validation_chi2_dist.pdf', bbox_inches='tight')
%history
%history?
%history -f plot_generators/plot_validation_chi2_dist.py
