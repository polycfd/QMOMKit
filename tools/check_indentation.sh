#!/bin/sh
## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2014 - 2018 by the deal.II authors
##
## This file was adapted for QMOMKit from Lethe:
##   https://github.com/chaos-polymtl/lethe/blob/master/contrib/utilities/check_indentation.sh
## which in turn adapted it from the deal.II library:
##   https://github.com/dealii/dealii/blob/master/contrib/utilities/check_indentation.sh
##
## Changes with respect to the Lethe version:
##   - The indent script is invoked from QMOMKit's tools/ directory instead
##     of contrib/utilities/.
##
## Detailed license information governing the original source code can be
## found in LICENSE.md and CONTRIBUTING.md at the top level directory of
## deal.II.
##
## -----------------------------------------------------------------------------

#
# This is a script that is used by the continuous integration servers
# to make sure that the currently checked out version of a git repository
# satisfies our indentation standards.
#
# It does so by running the 'indent-all' script (located in the current
# directory), calling 'git diff' to show what differences exist between
# the correctly indented code and what is in the git index (which is
# typically what is in the last commit), and then running a command
# that either returns success or failure, depending on whether or not
# there are differences. The continuous integration services return
# a failure code for a pull request if this script returns a failure.
#
# It requires clang-format 22.1, which './tools/download_clang_format'
# installs into './tools/programs'.
#

echo "Running indentation test"

./tools/indent-all || exit 1

git diff
git diff-files --quiet
