from mrjob.job import MRJob
from mrjob.step import MRStep

# This Job determines the number of triangle sin a graph using the algorithm
# Algortihm 3 MR-NodeIterator++(V,E) described in the document 
#"Counting Triangles and the Curse of the Last Reducer" 
# by Siddharth Suri Sergei and Vassilvitskii

class MRNodesNT(MRJob):
    
    ######################################
    ## Job to create the adjacency matrix
    ######################################
    def mapper_AGeneration(self, _, line):
        # This mapper assigns per every vertex, a 1 to every nodein the vextex.
        nodes_vertex = line.split()
        if (len(nodes_vertex) == 2):
            yield nodes_vertex[0], (nodes_vertex[1],1)
            yield nodes_vertex[1], (nodes_vertex[0],1)

    def combiner_AGeneration(self, key, values):
        # Combines the adjacencies of the nodes.
        yield key, list(values)

    def reducer_AGeneration(self, key, values):
        # Combines the adjacencies of the nodes again, all of them.
        final_list = []
        for value in values:
            final_list.extend(list(value))
        final_dict = dict(final_list)
        
        row_elem = [key, final_dict]
        # Each row contains on its data, the number of the row and the dictionary
        # of values in sparse way. [Nrow, Dict of cols]
        
        yield "A", row_elem  # So that the next reducer get
        
    ######################################
    ## Generate the 2-length paths
    ######################################
    # We need to have access to all rows of the adjacency matrix to create
    # the tasks so we implement a new returcer to join them
    def reducer_inter(self, key, values):
        val_list = list(values)
        Matrix_A = dict(val_list) # We make another dictionary with the cols
        A_keys = Matrix_A.keys()
        
        # Output all the pairs !! And their neighbours
        for k in A_keys:
            neighs = Matrix_A[k]
            neighs_keys = neighs.keys()
            Nneigh = len(neighs_keys)
            for i in range(Nneigh):
                for j in range(i,Nneigh):
                    nk1 = neighs_keys[i]
                    nk2 = neighs_keys[j]
                    if (nk1 != nk2):   
                        yield str(k) + "_" + str(nk1)+"_"+str(nk2), [nk1, Matrix_A[nk2]]

    ##################################################
    ## Check if you can form a triangle with every pair
    ##################################################
    def mapper_2Path(self, key, values):
        # This mapper outputs pairs of neighbours for each node.
        # So that latter the mapping can 
        lista = list(values)
#        print lista
        if (lista[0] in lista[1].keys()):
            yield None,1
            
    def reducer_2Path(self, key, values):
        yield "NT", sum(values)/3
    
    def steps(self):
        return [

            MRStep(mapper=self.mapper_AGeneration,
                   combiner=self.combiner_AGeneration,
                   reducer=self.reducer_AGeneration),
        
            MRStep(reducer=self.reducer_inter),

            MRStep(mapper=self.mapper_2Path,
                   reducer=self.reducer_2Path),
        ]
        
if __name__ == '__main__':
    MRNodesNT.run()
