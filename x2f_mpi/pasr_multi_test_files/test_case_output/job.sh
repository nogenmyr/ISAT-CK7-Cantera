#$ -V                      # Inherit the submission environment
#$ -cwd                    # Start job in  submission directory
#$ -N pasr_multi           # Job Name
#$ -pe 16way 32            # Requests 16 cores/node, 16 cores total
#$ -q development          # Queue name
#$ -l h_rt=01:30:00        # Run time (hh:mm:ss) - 1.5 hours
ibrun ./pasr_multi
