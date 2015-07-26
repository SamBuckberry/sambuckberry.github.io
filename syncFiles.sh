#! /bin/bash

# rsync files to update site
rsync -a _site/about/index.html index.html
rsync -a ~/Desktop/Dropbox/PhD/R_Projects/RUN-WGCNA/readme.md _posts/2015-07-25-runWGCNA.md

# Get the site ready to serve
jekyll serve

