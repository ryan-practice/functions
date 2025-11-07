#https://www.youtube.com/watch?v=kQ5QkN4Kx4Q

# $ git commit -m "first commit"
# [master (root-commit) caa4c98] first commit
# Committer: Sangston <rms2jg@virginia.edu>
#   Your name and email address were configured automatically based
# on your username and hostname. Please check that they are accurate.
# You can suppress this message by setting them explicitly. Run the
# following command and follow the instructions in your editor to edit
# your configuration file:
#
#   git config --global --edit
#
# After doing this, you may fix the identity used for this commit with:
#
#   git commit --amend --reset-author
#
# 7 files changed, 78 insertions(+)
# create mode 100644 .Rbuildignore
# create mode 100644 .gitignore
# create mode 100644 DESCRIPTION
# create mode 100644 NAMESPACE
# create mode 100644 R/hello.R
# create mode 100644 functions.Rproj
# create mode 100644 man/hello.Rd

# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}
