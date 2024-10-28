# RooFit-to-MPL

to get the package, do \
`$ git clone git@github.com:aburke1605/RooFit-to-MPL.git`

add the package to your python interpreter by doing \
`$ export PYTHONPATH=$PYTHONPATH:path/to/RooFit-to-MPL`

allow for intellisense recognition in VS code by doing the following in your base directory \
`$ echo "PYTHONPATH=$PYTHONPATH:path/to/RooFit-to-MPL" > .env`



to run the example do
```shell
$ python -i
>>> from Utils import Example
>>> Example()
```
which generates the example plots `"mpl_figure_plot.pdf"` and `"root_canvas_plot.pdf"`:
<img src="plots/root_canvas_plot.png" alt="ROOT plot" width="49%"> <img src="plots/mpl_figure_plot.png" alt="matplotlib plot" width="49%">

a more extensive use of the plotter can be found at <a href="https://gitlab.cern.ch/lhcb-charm/kpi_dt/-/blob/master/analysis/fit_delta_mass.py?ref_type=heads">this link</a>
(requires a CERN account)