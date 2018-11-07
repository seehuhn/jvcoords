Developers' Notes
=================

Style
-----

Several popular style guides for R code can be found on the internet,
but sadly, there is no agreement on how functions and variables should
be named.  After some consideration, I decided to follow the majority
choice of packages on CRAN, as reported by Rasmus Bååth's article `The
State of Naming Conventions in R`_.  Thus, The convention I use now
is:

  - ``lowerCamelCase`` for function names, and
  - ``period.separated`` for variable names and parameter names.

Except for this convention, I try to follow Google's `R Style Guide`_.

.. _`The State of Naming Conventions in R`: https://journal.r-project.org/archive/2012-2/RJournal_2012-2_Baaaath.pdf
.. _`R Style Guide`: https://google.github.io/styleguide/Rguide.xml

Release Checklist
-----------------

  - update ``Version`` and ``Date`` in the DESCRIPTION file
  - check ``git status'' and remove any unnecessary files
  - check that the current version of R is installed
  - check that the output of ``library(help=jvcoords)`` looks reasonable
  - build source package
  - check that the contents of the source package look reasonable
  - run ``R CMD check --as-cran``
