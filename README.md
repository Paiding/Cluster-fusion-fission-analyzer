# Cluster-fusion-fission-analyzer

## Tutorial  
1. Prepare a trajectory, put the .gro or .xyz ... files in /conf.  
2. Run clustering(), read the files and cluster data will be placed in /Index.  
3. Set the value of max cluster size from last step,set the size of box, set the number of frames, run splitTraceBack() and rate coefficients alpha and beta will be placed in /rate.  
4. Run splitAverageRate(), compute the average rate matrix.  
5. Edit the file name in NeoVisual.py, and draw the pcolor figure.  
