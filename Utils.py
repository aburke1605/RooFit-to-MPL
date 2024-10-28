import ROOT as R

from PDF import (
    GaussianPDF,
    JohnsonPDF,
    ArgusPDF,
)
from Model import Model

def Example(ran: tuple=(140, 155), nbins: int=100):
    """
    Example usage:
        test the addition of 2 PDFs with a fraction between them and a separate PDF for bkg

        first generate and fit some data with RooFit
        then pass the fitted values and binned data to mpl

        should compare the root canvas plot and the mpl figure plot,
        and also the chi2s
    """

    ####################
    ### first in RooFit:

    x = R.RooRealVar("x", "x", *ran)

    ## define fitting PDFs
    #     Gaussian:
    Gmu = R.RooRealVar("Gmu", "Gmu", 145.45, *ran)
    Gsigma = R.RooRealVar("Gsigma", "Gsigma", 0.17, 0.01, 2)
    G = R.RooGaussian("G", "G", x, Gmu, Gsigma)
    #     Johnson:
    Jmu = R.RooRealVar("Jmu", "Jmu", 145.39, *ran)
    Jsigma = R.RooRealVar("Jsigma", "Jsigma", 0.25, 0.01, 2)
    Jgamma = R.RooRealVar("Jgamma", "Jgamma", -0.23, -5, 0)
    Jdelta = R.RooRealVar("Jdelta", "Jdelta", 0.78, 0, 10)
    J = R.RooJohnson("J", "J", x, Jmu, Jsigma, Jgamma, Jdelta)
    #     Total signal:
    frac = R.RooRealVar("frac", "frac", 0.32, 0, 1)
    signal = R.RooAddPdf("signal", "signal", R.RooArgList(G,J), R.RooArgList(frac))
    #     Argus:
    c = R.RooRealVar("c", "c", -3.59, -10, 0)
    p = R.RooRealVar("p", "p", 0.92, 0, 10)
    x0 = R.RooRealVar("x0", "x0", 139.57061)
    flipped_x = R.RooFormulaVar("flipped_x", "@1-(@0-@1)", R.RooArgList(x, x0))
    background = R.RooArgusBG("A", "A", flipped_x, x0, c, p)

    ## generate some pseudo-data from PDFs
    n_sig = 100_000
    n_bkg = 25_000
    roo_data = G.generate(R.RooArgSet(x), NumEvents=int(frac.getValV()*n_sig))
    roo_data.append(J.generate(R.RooArgSet(x), NumEvents=int((1-frac.getValV())*n_sig)))
    roo_data.append(background.generate(R.RooArgSet(x), NumEvents=n_bkg))

    ## now do a fit to pseudo-data
    # alter initial guesses slightly
    Gmu.setVal(145.4)
    Gsigma.setVal(0.2)
    Jmu.setVal(145.4)
    Jsigma.setVal(0.3)
    Jgamma.setVal(-0.3)
    Jdelta.setVal(1.0)
    c.setVal(-3.4)
    p.setVal(1)
    # define yields for extended maximum likelihood fit
    Nsig = R.RooRealVar("Nsig", "Nsig", 50_000, 10_000, 500_000)
    Nbkg = R.RooRealVar("Nbkg", "Nbkg", 20_000, 10_000, 50_000)
    model = R.RooAddPdf("model", "model", R.RooArgList(signal, background), R.RooArgList(Nsig, Nbkg))
    model.fitTo(roo_data, Extended=True, Save=True, PrintLevel=-1)

    ## plot the pseudo-data and fit results
    frame = x.frame()
    # plot data points first to correctly normalise PDFs after
    roo_data.plotOn(frame, Binning=(nbins, *ran))
    # plot total PDF
    model.plotOn(frame)
    # plot PDF componente
    model.plotOn(frame, R.RooFit.Components("G"), R.RooFit.LineStyle(2), R.RooFit.LineColor(R.kSpring))
    model.plotOn(frame, R.RooFit.Components("J"), R.RooFit.LineStyle(2), R.RooFit.LineColor(R.kRed))
    model.plotOn(frame, R.RooFit.Components("A"), R.RooFit.LineStyle(2), R.RooFit.LineColor(R.kMagenta+2))
    # plot data points again so they're on top
    roo_data.plotOn(frame, Binning=(nbins, *ran))
    canvas = R.TCanvas("canvas", "canvas", 1125, 800)
    pad = R.TPad("pad", "pad", 0, 0.0, 1.0, 1)
    pad.Draw()
    pad.SetLogy()
    pad.cd()
    frame.Draw()
    canvas.SaveAs("./root_canvas_plot.pdf")

    ## get the chi2 from RooFit
    nparams = model.getParameters(roo_data).selectByAttrib("Constant", False).getSize()
    TH_data = roo_data.createHistogram("x", Binning=(nbins, *ran))
    roo_data_hist = R.RooDataHist("roo_data_hist", "roo_data_hist", R.RooArgList(x), R.RooFit.Import(TH_data))
    chi2 = model.createChi2(roo_data_hist, R.RooFit.Extended(True))
    print(f"{chi2.getValV():.2f}/{nbins-nparams}")



    ##########################
    ### now the important bit,
    ### use matplotlib.pyplot:

    G = GaussianPDF(Gmu.getValV(), Gsigma.getValV())
    Gnorm = G.get_normalisation(nbins, Nsig.getValV(), *ran)
    J = JohnsonPDF(Jmu.getValV(), Jsigma.getValV(), Jgamma.getValV(), Jdelta.getValV())
    Jnorm = J.get_normalisation(nbins, Nsig.getValV(), *ran)
    A = ArgusPDF(c.getValV(), p.getValV())
    Anorm = A.get_normalisation(nbins, Nbkg.getValV(), *ran)

    model = Model([G, J, A], [frac.getValV(), 1.-frac.getValV(), 1.], [Gnorm, Jnorm, Anorm], ["Gaussian", "Johnson", "Argus"])
    model.import_data(TH_data)
    model.plot({
        "x_label": r"$\Delta m$",
        "x_unit": r"MeV$/c^{2}$",
        "legend_title": "LHCb $5.6\,\mathrm{fb}^{-1}$",
        "nparams": nparams,
        "logy": True,
        "save_to": "./mpl_figure_plot.pdf",
    })
