# Place holder for maintaining the results directory.

# For making movies from a list of .png files then use the command
mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800 -nosound -o eigs.avi
