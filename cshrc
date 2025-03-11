#ident  "@(#).cshrc     ver 1.0     Aug 20, 1996"
# Default user .cshrc file.
#
# This file is executed each time a shell is started.
# This includes the execution of shell scripts.


#####
# Source the site-wide syscshrc file.
# The syscshrc file defines some needed aliases (setup amd unsetup)
# and environment variables (PATH and MANPATH).  This line
# should not be deleted.  You do, however, have a choice of
# syscshrc files.  Uncomment the one that you prefer.
#####
#source /site/env/syscshrc    # Searches /usr/local/bin first.
#source /site/env/syscshrc.alt   # Searches /usr/local/bin last.

module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles

#####
# Set up the shell environment.  You may comment/uncomment
# the following entries to meet your needs.
#####
# Number of commands to save in history list.
set history=500
 
# Number of commands to save in ~/.history upon logout.
set savehist=500

# Notify user of completed jobs right away, instead of waiting
# for the next prompt.
set notify

# Don't redirect output to an existing file.
# CAD NOTE!  This must be commented out for proper ME10 functionality!!
set noclobber

# Set the file creation mode mask (default permissions for newly created files).
umask 022

#####
# Define your aliases.
#####
alias       h       history
alias       d       dirs
alias       pd      pushd
alias       pd2     pushd +2
alias       po      popd
alias       m       more
alias       ls      'ls -F'
alias	    GROOVY myClara/plugins/clas12/bin/run-groovy
alias 	    matlab /apps/matlab/matlab-R2012a/bin/matlab
alias       mathematica /apps/mathematica/10.4/bin/MathKernel
alias       hipodump /u/home/thayward/myClara/plugins/clas12/bin/hipo-utils -dump
alias       hipoutils /u/home/thayward/myClara/plugins/clas12/bin/hipo-utils
alias	    swif /site/bin/swif2



#####
# User specific additions should be added below.
#####


#### Environment set up; last updated July 17th, 2024
module load python/3.12.4
module load root/6.32.02
module load groovy/4.0.3


setenv PATH "/home/thayward/.local/bin:$PATH"