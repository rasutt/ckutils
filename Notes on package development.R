# Following example from https://r-pkgs.org/whole-game.html

# Load devtools package for writing packages
library(devtools)

# Load all package functions - simulates the process of building, installing,
# and attaching (using library()) the package. Ctrl + Shift + L
load_all()

# Set random seed for testing
set.seed(1)

# Need a global variable setting births to be stochastic
stch.bths = T

# Simulate one population and study
pop_stud = SimPopStud(
  phi = 0.9, lambda = 1.05, N.init = 20, hist.len = 20, srvy.yrs = 20, k = 1,
  f.year = 20, p = 0.5, L = 10, imaf = 0.5, clvng.p = 0, tmp.emgn = 0,
  alpha = 5, clvng.ints = F
)

# Look at it
head(pop_stud, 1)
names(attributes(pop_stud))

# Check the package works. Ctrl + Shift + E
check()

# Convert roxygen comments into R documentation, and update namespace. Ctrl +
# Shift + D
document()

# Check the help file
?SimPopStud

# Install package. Ctrl + Shift + B
install()

# Load package!
library(ckutils)

# Make test file
use_test("simpop.R")

# Run tests. Ctrl + Shift + T
test()

# Add stats package to imports section in description file (required by CRAN).
# Need to use stats::function_name too I think, and in general for non-base
# packages
use_package("stats")
