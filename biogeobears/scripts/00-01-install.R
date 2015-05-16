#######################################################
# Installing BioGeoBEARS
#######################################################
# Uncomment this command to get everything
# I request that you use the "0-cloud" R repository at "http://cran.rstudio.com" as it is
# the only one that keeps download statistics. However, this is up to you.
#######################################################
#
# I recommend you install optimx first, with the dependencies=TRUE command, as
# optimx has many dependencies that people sometimes forget to install.

# Install optimx
install.packages("optimx", dependencies=TRUE, repos="http://cran.rstudio.com")

# Install BioGeoBEARS from CRAN 0-cloud:
install.packages("BioGeoBEARS", dependencies=TRUE, repos="http://cran.rstudio.com")

#######################################################
#
# IMPORTANT NOTE: ONCE BIOGEOBEARS IS INSTALLED,
# USE THE SOURCE() COMMANDS IN THE EXAMPLE SCRIPT AT
#
# http://phylo.wikidot.com/biogeobears#script
#
# TO GET BUG FIXES AND UPDATES.
#
# YOU WILL HAVE TO RUN THESE SOURCE COMMANDS
# AFTER EACH TIME YOU RUN library(BioGeoBEARS).
# (you can either source them online, or download the .R files
#  and source them from your hard drive)
#
#######################################################
