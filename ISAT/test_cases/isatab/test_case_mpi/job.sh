#$ -V                      # Inherit the submission environment
#$ -cwd                    # Start job in  submission directory
#$ -N test_mpi             # Job Name
#$ -j y                    # combine stderr & stdout into stdout    
#$ -pe 8way 24             # Requests 8 cores/node, 2 nodes total
#$ -q development          # Queue name
#$ -l h_rt=01:00:00        # Run time (hh:mm:ss) - 1.0 hours
ibrun ./test_mpi
