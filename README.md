# StructuralCoexistence

This is a Shiny app that allows the user to interactively explore the structural approach to coexistence framework described by Saavedra et al. (2017). For the app, I have simply extracted analysis and plotting code from a [Dryad repository](https://doi.org/10.5061/dryad.v9f5s) that accompanies the paper, making modifications for interactivity.

To use the app, clone this repo locally or try it out online at: https://ecodynamics.shinyapps.io/StructuralCoexistence/

Any errors are most likely mine—please submit an Issue or contact me directly if you find one.

### Citation:

Saavedra, S., Rohr, R. P., Bascompte, J., Godoy, O., Kraft, N. J. B. and Levine, J. M. (2017), A structural approach for understanding multispecies coexistence. Ecol Monogr, 87: 470–486. [doi:10.1002/ecm.1263](http://dx.doi.org/10.1002/ecm.1263)

### To do:
- explainer text on main panel
- legend for plot symbols
- better labels on cone plot
- row/column headers for the competition coefficient matrix
- graphically indicate which species pairs are feasible
- figure out subscript in RenderTable (expressions?)
- tab for 4 species (better way to show 3d plots?)
- remove +/- on sidebarPanel (ugh Javascript in shinyIncubator)
- speed up server-side
- re-do coordinate system to make ternary plot equilateral (will require re-write of most functions, possibly additional package dependency to make the Cartesian-to-ternary conversions easier)
