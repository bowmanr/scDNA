#!/bin/bash

set -e

REPO_URL="https://github.com/bowmanr/scDNA.git"
BRANCH="dev_cleanup"

echo "Checking Git..."
if ! command -v git >/dev/null 2>&1; then
    echo "Git is not installed"
    echo "Please install GIT from https://git-scm.com/downloads"
    exit 1
fi

if ! command -v git-lfs >/dev/null 2>&1; then
    
  echo "Git LFS not found."
  
  if command -v brew >/dev/null 2>&1; then
    echo "installing Git LFS with homebrew..."
    brew install git-lfs
  else
    echo ""
    echo "homebrew not found"
    echo "https://brew.sh"
    exit 1
  fi
fi

git lfs install

TMPDIR=$(mktemp -d)
echo "Using temporary directory:"
echo "$TMPDIR"

echo "Cloning repository..."
git clone --branch "$BRANCH" "$REPO_URL" "$TMPDIR/package"

cd "$TMPDIR/package"

echo "Downloading LFS files..."
git lfs pull

echo "Installing package..."

Rscript -e '
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install(".", upgrade = "never",dependencies = TRUE)
'

echo "Cleaning up..."
cd /
rm -rf "$TMPDIR"

echo "Installation complete."