(TeX-add-style-hook
 "VOF_adv"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "amsmath"
    "tikz")
   (LaTeX-add-labels
    "eq:cut-1"
    "eq:cut1"
    "eq:wy"
    "eq:CIAM"
    "min"
    "taylor"))
 :latex)

