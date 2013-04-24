# This runs the full-throated R-lang complete package system check,
# which also checks for code issues, missing dependencies, and other
# important issues. Will take much longer to run than the "unit tests only" script.
# Note that this will also run the unit tests though...
R CMD check biom
