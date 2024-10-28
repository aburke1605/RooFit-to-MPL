import ROOT as R
import numpy as np
from scipy.stats import norm

from pathlib import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
plt.style.use(f"{Path(__file__).parent}/RF.mplstyle")
# this is optional:
plt.style.use("seaborn-v0_8-muted")

import pickle

from PDF import PDF

class Model:

    def __init__(self, pdfs: list[PDF], fracs: list[float], norms: list[float]|None=None, names: list[str]|None=None) -> None:
        """
        Assumption: pdfs, fracs and norms are ordered correctly
        """

        assert len(pdfs) == len(fracs), \
            "Must provide N fractions for N PDFs"

        self._pdfs = pdfs
        self._fracs = fracs
        
        if norms is not None:
            assert len(pdfs) == len(norms), \
                "Must provide N normalisations for N PDFs"
        self._norms = norms

        self._names = []
        if names is not None:
            self._names = names

        self._has_data = False # place holder

    def evaluate(self, x: float|np.ndarray) -> float|np.ndarray:
        assert self._norms!=None, \
            Exception("Normalisations have not been set")

        y = 0. if type(x) == float else np.zeros_like(x)

        # loop over PDFs summing each component
        for i, pdf in enumerate(self._pdfs):
            frac = self._fracs[i]
            norm = self._norms[i]

            y += pdf.evaluate(x, frac * norm)
        
        return y


    def import_data(self, data: R.TH1) -> None:
        """
        Only support importing as TH1 for now
        """

        self._nbins = data.GetNbinsX()
        self._bin_edges = np.ones(self._nbins+1, dtype=float)
        self._bin_widths = np.ones(self._nbins, dtype=float)

        # convert root histogram to numpy array
        self._counts = np.ndarray(self._nbins)
        self._errors = np.ndarray(self._nbins)
        for i in range(self._nbins):
            self._bin_edges[i] = data.GetBinLowEdge(i+1)
            self._bin_widths[i] = data.GetBinWidth(i+1)
            self._counts[i] = data.GetBinContent(i+1)
            self._errors[i] = data.GetBinError(i+1)

        self._bin_edges[-1] = data.GetBinLowEdge(self._nbins) + data.GetBinWidth(self._nbins)
        self._bin_centres = np.array([(self._bin_edges[i] + self._bin_edges[i+1]) / 2 for i in range(self._nbins)])

        self._min = self._bin_edges[0]
        self._max = self._bin_edges[-1]

        # if no normalisations have been passed,
        # i.e. it's not an extended likelihood fit,
        # then just normalise to the data passed
        # as per regular fits
        if self._norms == None:
            self._norms = []
            for i, pdf in enumerate(self._pdfs):
                self._norms.append(pdf.get_normalisation(self._nbins, np.sum(self._counts), self._bin_edges[0], self._bin_edges[-1]))

        self._has_data = True


    def plot(self, config: dict={}) -> float:

        #######
        # 
        # TODO
        #
        # * only add chi2 if requested, since most fits are NLL
        # * plot without pulls too
        # * ability to set y range
        #######

        ###############
        ## configration
        allowed_keys = ["x_label", "y_label", "x_unit", "y_unit", "domain", "legend_title", "nparams", "y_lim", "logy", "no_chi2", "no_pulls", "no_data", "no_pdfs", "signals", "backgrounds", "left_legend", "pull_dist", "save_to", "pickle_it"]
        for key in config.keys():
            assert key in allowed_keys, "Unrecognised configuration parameter"
        x_label = config.get("x_label", None)
        y_label = config.get("y_label", "Candidates")
        x_unit = config.get("x_unit", "")
        y_unit = config.get("y_unit", None)
        domain = config.get("domain", None)
        legend_title = config.get("legend_title", "")
        nparams = config.get("nparams", 0)
        y_lim = config.get("y_lim", None)
        logy = config.get("logy", False)
        no_chi2 = config.get("no_chi2", False)
        no_pulls = config.get("no_pulls", False)
        no_data = config.get("no_data", False)
        no_pdfs = config.get("no_pdfs", False)
        signals = config.get("signals", None)
        backgrounds = config.get("backgrounds", None)
        left_legend = config.get("left_legend", False)
        pull_dist = config.get("pull_dist", False)
        save_to = config.get("save_to", None)
        pickle_it = config.get("pickle_it", False)
        if no_data or no_pdfs:
            no_chi2 = True
            no_pulls = True


        #####################
        ## prepare the canvas
        fig, fit_ax = plt.subplots()
        if not no_pulls:
            divider = make_axes_locatable(fit_ax)
            pull_ax = divider.append_axes("bottom", size="25%", pad=0.35, sharex=fit_ax)
            plt.setp(fit_ax.get_xticklabels(), visible=False)
        if logy:
            fit_ax.semilogy()


        #####################
        ## plot the fit first

        # if no data has been imported, make some
        # fake histogram so everything else works
        if not self._has_data:
            assert type(domain)==tuple, \
                Exception("Pass a plotting domain tuple if not importing any data")
            self._min = domain[0]
            self._max = domain[1]
            empty_TH = R.TH1F("", "", 1, self._min, self._max)
            self.import_data(empty_TH)
            no_chi2 = True
            no_pulls = True
            no_data = True

        # plot data points
        if not no_data:
            assert self._has_data, \
                Exception("No imported data. Pass \"no_data\" in config to plot without")
            fit_ax.errorbar(self._bin_centres, self._counts, xerr=self._bin_widths/2, yerr=self._errors, markersize=5, elinewidth=1, capsize=1, capthick=1, fmt="k.", label="Data", zorder=100)

        # plot PDFs
        if not no_pdfs:
            assert self._norms!=None, \
                Exception("Normalisations have not been set")
            lin = np.linspace(self._min, self._max, 1_000)
            if len(self._pdfs) > 1:
                # plot total
                fit_ax.plot(lin, self.evaluate(lin), label="Fit")

                prop_cycle = plt.rcParams["axes.prop_cycle"]
                colors = prop_cycle.by_key()["color"]

                # plot components
                if signals is not None and backgrounds is not None:
                    for pdfs, label, zorder in zip([signals, backgrounds], ["Signal", "Background"], [1, 0]):
                        for pdf in pdfs:
                            assert pdf in range(len(self._pdfs)), f"Index {pdf} out of range"
                        fit_ax.plot(lin, np.sum([self._pdfs[i].evaluate(lin, self._fracs[i] * self._norms[i]) for i in pdfs], axis=0), linestyle="dashed", label=label, zorder=zorder)
                elif signals is None and backgrounds is not None:
                    for pdf in backgrounds:
                        assert pdf in range(len(self._pdfs)), f"Index {pdf} out of range"
                    fit_ax.fill_between(lin, 0, np.sum([self._pdfs[i].evaluate(lin, self._fracs[i] * self._norms[i]) for i in backgrounds], axis=0), color=colors[2], label="Background")
                else:
                    for i, pdf in enumerate(self._pdfs):
                        label = self._names[i] if i < len(self._names) else f"PDF {i+1 - len(self._names)}"
                        fit_ax.plot(lin, pdf.evaluate(lin, self._fracs[i] * self._norms[i]), linestyle="dashed", label=label)
            else:
                label = self._names[0] if len(self._names) != 0 else "PDF"
                fit_ax.plot(lin, self._pdfs[0].evaluate(lin, self._fracs[0] * self._norms[0]), label=label)


        fit_ax.set_xlim(self._min, self._max)
        if logy:
            fit_ax.set_ylim(0.1+0.5*min(self._counts), 2*max(self._counts))
        if y_lim is not None:
            assert len(y_lim)==2, "Range should contain 2 values (min and max)"
            fit_ax.set_ylim(y_lim[0], y_lim[1])
        fit_ax.set_ylabel(rf"{y_label} / ({self._bin_widths[0]:.2f}{' ' + x_unit if x_unit!='' else ''})" if y_unit is None and (self._bin_widths==self._bin_widths[0]).all() else rf"{y_label} [{y_unit}]")


        ###########################
        ## calculate chi2 and pulls

        chi2 = 0.
        good_pulls = np.zeros_like(self._counts) # plotted in black
        bad_pulls = np.zeros_like(self._counts) # plotted in red
        for i in range(self._nbins):
            if self._counts[i] == 0:
                continue
            f = 0.
            n_steps = 500
            bin_range = np.linspace(self._bin_edges[i], self._bin_edges[i+1], n_steps, endpoint=False)
            for x in bin_range:
                f += self.evaluate(x) / n_steps
            good_pulls[i] = (self._counts[i]-f)/self._errors[i]
            if abs(good_pulls[i]) >= 3:
                bad_pulls[i] = good_pulls[i]
            chi2 += np.power(good_pulls[i], 2)

        # add chi2 to plot
        if not no_chi2 and nparams is not None:
            result_string = (rf"$\chi^{{2}} / $nDOF $ = {chi2:.2f} / {self._nbins - nparams}$")
            fit_ax.errorbar(np.NaN, np.NaN, color="none", label=result_string, zorder=101)

        # add legend to plot
        # fit_ax.legend(title=legend_title, loc="upper left" if left_legend else "upper right", bbox_to_anchor=(0,1) if left_legend else (1, 1), ncol=1)
        fit_ax.legend(title=legend_title, loc="upper left" if left_legend else "upper right", bbox_to_anchor=(0,1) if left_legend else (1, 1), ncol=1, fontsize=20, title_fontsize=20)


        ######################
        ## plot the pulls next
        if no_pulls:
            fit_ax.set_xlabel(f"{x_label}{'' if x_unit=='' else ' ['+x_unit+']'}")
        else:
            pull_ax.hlines(0, self._min, self._max, linewidth=0.1, color="k")
            pull_ax.bar(self._bin_centres, good_pulls, width=self._bin_widths, align="center", color="k")
            pull_ax.bar(self._bin_centres, bad_pulls, width=self._bin_widths, align="center", color="r")

            pull_ax.set_xlim(self._min, self._max)
            pull_ax.set_ylim(-5, 5)
            pull_ax.set_ylabel(rf"Pull [$\sigma$]")
            if x_label is not None:
                pull_ax.set_xlabel(f"{x_label}{'' if x_unit=='' else ' ['+x_unit+']'}")

            if pull_dist:
                pulls_fig, pulls_ax = plt.subplots()
                nbins = 25
                pulls_hist, edges = np.histogram(good_pulls, range=(-5,5), bins=nbins, density=True)
                old_max = max(pulls_hist)
                pulls_hist, edges = np.histogram(good_pulls, range=(-5,5), bins=nbins)
                new_max = max(pulls_hist)
                centres = [(edges[i] + edges[i+1]) / 2 for i in range(nbins)]
                width = (edges[1] - edges[0]) / 2

                pulls_ax.plot(np.linspace(-5,5,100), new_max/old_max*norm.pdf(np.linspace(-5,5,100)), label="Normal distribution")
                pulls_ax.errorbar(centres, pulls_hist, xerr=width, yerr=np.sqrt(pulls_hist), fmt="k.", label="Fit pulls")

                pulls_ax.legend(loc="upper right", bbox_to_anchor=(1, 1), ncol=1, fontsize=20, title_fontsize=20)
                pulls_ax.set_xlabel(r"Pull [$\sigma$]")
                pulls_ax.set_ylabel(rf"Entries / ({width:.2f}$\sigma$)")

                if save_to is not None:
                    file_name, extension = save_to.split(".")
                    if pickle_it == True:
                        pickle.dump(pulls_fig, open(f"{file_name}_pulls.fig.pickle", "wb"))
                    else:
                        pulls_fig.savefig(f"{file_name}_pulls.{extension}", bbox_inches="tight")


        ################
        ## finally, save
        if save_to is not None:
            if pickle_it == True:
                file_name, extension = save_to.split(".")
                pickle.dump(fig, open(f"{file_name}.fig.pickle", "wb"))
            else:
                fig.savefig(save_to, bbox_inches="tight")
        else:
            plt.show()

        return chi2 / (self._nbins - nparams)