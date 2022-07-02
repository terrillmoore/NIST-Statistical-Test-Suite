#!/bin/bash
#
# Module: setup.sh
#
# Function:
#	Create additional directories and prepare for building after checkout.
#
# Copyright and license info:
#	See accompanying LICENSE.md file.
#
# Author:
#	Terry Moore, January 2019
#
# Description:
#	See the --help output.
#

PNAME=$(basename "$0")

function _help {
cat <<EOF
$PNAME is used after checking out a clean repository, to make sure that
all the required directories are present. sts/assess and sts/makefile
do not have the ability to create missing directories, and git doesn't
record empty directories. I'd like to keep the NIST distribution pristene
(at least for now), so you need to run this script after a clean checkout,
before doing your first build. It won't hurt to re-run it.

For testing, $PNAME --clean will remove all the subdirectories that this
script creates.
EOF
}

if [ X"$1" == "X--help" ]; then
	_help
	exit 0
fi

function _error {
	echo "$PNAME:" "$@" 1>&2
	exit 1
}

ROOTDIR=$(dirname "$0")/sts

if [ ! -d "$ROOTDIR" ]; then
	_error "STS distribution directory not found:" "$ROOTDIR/sts"
fi

function _prefix {
	local i PREFIX
	PREFIX="$1" ; shift
	if [ ! -d "$PREFIX" ]; then
		_error "Not a directory: $PREFIX"
	fi
	for i in "$@" ; do echo "$PREFIX/$i" ; done
}

# use an array for SUBDIRS, in case ROOTDIR expansion contains spaces.
typeset -a SUBDIRS
SUBDIRS+=("$(_prefix "$ROOTDIR" obj)")
for i in AlgorithmTesting BBS CCG G-SHA1 LCG MODEXP MS QCG1 QCG2 XOR; do
	SUBDIRS+=("$(_prefix "$ROOTDIR/experiments" $i)")
done

if [[ X"$1" == "X--clean" ]]; then
	echo "Cleaning up directory tree."
	#
	# expand one word per array entry, in case ROOTDIR expansion contains spaces
	#
	rm -rf "${SUBDIRS[@]}"
	exit $?
fi

if [[ $# -ne 0 ]]; then
	_error "Unrecognized arguments; use --help to get help"
fi

echo "Setting up directories in $ROOTDIR/experiments."

# expand one word per entry, in case ROOTDIR expansion contains spaces
for i in "${SUBDIRS[@]}" ; do
	if [ ! -d "$i" ]; then
		mkdir -p "$i" || _error "Can't create dir: $i"
        echo "Created $i."
    else
        echo "$i already exists."
	fi
done

echo "Creating the subdirectories via $ROOTDIR/experiments/create-dir-script."
	(
	cd "$ROOTDIR/experiments" || _error "Can't cd: $ROOTDIR/experiments"
	# this is a hack; we use the first directory created by create-dir-script
	# as the flag to skip create-dir-script.
	if [ ! -d AlgorithmTesting/Frequency ]; then
		if ! ./create-dir-script; then
			echo "$PNAME: some directory creations failed; this probably isn't a problem, but please check!"
		else
			echo "create-dir-script succeeded."
		fi
	else
		echo "Skipping create-dir-script. It appears that it was already run?"
	fi
	) || exit $?

echo "Directories are set up. Change directory to $ROOTDIR and say 'make'!"

exit 0
