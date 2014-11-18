#!/usr/bin/env bash -xe

git submodule update --init --recursive
# Remove the template remote
git remote remove origin
# Configure IPYNB output filtering
git config --local filter.dropoutput_ipynb.clean scripts/utils/ipynb_output_filter.py
git config --local filter.dropoutput_ipynb.smudge cat
# Link README to project notes, instead of template notes.
unlink README.md
ln -s NOTE.md README.md
# Avoid accidentally re-running this initialization script.
chmod -x "$0"
# Remove all of the commits after the first, leaving files intact,
# add files created/changed during initialization,
# and amend the first commit with everything else.
# TL;DR Squash everything to a single first commit.
git reset --soft $(git rev-list --max-parents=0 HEAD)
git add -A
git commit --amend -em "Clean project.  Let's get started!"