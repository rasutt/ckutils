# Following example from https://r-pkgs.org/whole-game.html

# Load devtools package for writing packages
# library(devtools)

# Set devtools to be attached automatically
# use_devtools()

devtools::dev_sitrep()

use_build_ignore("dev_notes.R")

# keyboard shortcuts, press Alt + Shift + K
# Command palette Ctrl + Shift + P

# Make new function file
use_r("sim_pop")
use_r("plot_exp_pop")
rename_files("sim_pop", "sim_pop_study")

# Load all package functions - simulates the process of building, installing,
# and attaching (using library()) the package. Ctrl + Shift + L
load_all()

# Try new functions in readme/tests

# Check the package works. Ctrl + Shift + E. Do this often!
check()

# Insert roxygen skeleton for function documentation. Inside function code Ctrl
# + Alt + Shift + R

# Put the dev code into the examples section

# Convert roxygen comments into R documentation, and update namespace. Ctrl +
# Shift + D
document()

# Check again
check()

# Check the help files
?plot_exp_pop
?sim_pop_study

# Install package. Ctrl + Shift + B
install()

# Load package!
library(ckutils)

# Make test files
usethis::use_testthat(3)
use_test("sim_pop_study")
use_test("plot_exp_pop")
use_test("find_exp_ns_kps")

# Adapt the example code as tests
# expect_equal includes numerical tolerance, expect_identical does not

# use_data()

# Run tests.
test_active_file() # Ctrl + t
test() # Ctrl + Shift + T

# Add packages to imports section in description file (required by CRAN). Need
# to use package::function_name for non-base packages
use_package("stats")
use_package("graphics")
use_package("utils")

# Check some more
check()

# Put examples in readme (https://r-pkgs.org/man.html#child-documents avoids
# repeating examples so much. I still think writing vignettes seems like a
# natural way to develop a package, cause you bring the functions together. But
# a website seems the main thing?)

# Knit readme with latest version of package
build_readme()

# Make website
# use_pkgdown_github_pages()
build_site()

# Git commit. ctrl + Alt + M
