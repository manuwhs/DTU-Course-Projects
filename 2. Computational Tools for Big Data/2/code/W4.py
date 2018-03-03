import pickle_lib as pklib
import CDBSCAN as CDBSCAN

problem_indx = 4
files_list = ["./data/data_10points_10dims.dat",
              "./data/data_100points_100dims.dat",
              "./data/data_1000points_1000dims.dat",
              "./data/data_10000points_10000dims.dat",
              "./data/data_100000points_100000dims.dat"]

eps_list = [0.4, 0.3, 0.15, 0.15, 0.15]
data = pklib.load_pickle (files_list[problem_indx],verbose = 1);

myDBSCAN = CDBSCAN.CDBSCAN(eps=eps_list[problem_indx], MinPts=2, reuse_computed = 0)
myDBSCAN.fit(data)

#print myDBSCAN.Already_Calculated_Regions[1]
print "Number of clusters " + str(myDBSCAN.K + 1)
myDBSCAN.print_clusters_sizes()

# 28470, 28121, 9, 19, 5, 13, 28412, 8, 2, 17, 8, 15, 16, 16, 16, 14, 14, 6, 13, 7, 21, 15, 15, 21, 10, 13, 10, 8, 12, 25, 13, 21