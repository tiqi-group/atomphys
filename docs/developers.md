{% include-markdown "../CONTRIBUTING.md" %}


This package is based on the package written by Matt Grau 'atomphys'. In order to change few data-structures and simplify overall package structure we (W.Adamczyk and C.Mordini) decided to rebuild it. This allowed us to add qutip helper functions, integrate other databases, and add visualising tools. States, Atoms and Transitions are now also not using registries, but graphs instead. P.Leindecker contributed valuable advice of how to build the package such that its API is well integrateable to web-development.

Wojciech Adamczyk <wadamczyk@phys.ethz.ch>,
Dr Carmelo Mordini <cmordini@phys.ehtz.ch>, Philipp Leindecker <pleindecker@phys.ethz.ch>,
Prof Matt Grau <mgrau@odu.edu>, Prof Jonathan Home <jhome@phys.ethz.ch>