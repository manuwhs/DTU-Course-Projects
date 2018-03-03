from mrjob.job import MRJob
from mrjob.step import MRStep

# This Job determines the number of triangle sin a graph using
# the trace of the A^3 (A = Adjacency matrix)

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
        # The combiner does need to do anything, the nodes will have all
        # its nodes joined already
        yield key, list(values)

    def reducer_AGeneration(self, key, values):
        # Same here
        final_list = []
        for value in values:
            final_list.extend(list(value))
        final_dict = dict(final_list)
        
        row_elem = [key, final_dict]
        # Each row contains the number of the row and the dictionary
        # of values in sparse way. [Nrow, Dict of cols]
        
        yield "Matrices", row_elem  # So that the next reducer get
        
    ######################################
    ## Matrix multiplication 
    ## We are gonna multiply the matrix A*A in sparse way,
    ######################################
    def reducer_MatrixMult(self, key, values):
        # This reducer output a par of (row_1, col_2) to multiply
        # by the combiner so that it outputs one element of the
        # output matrix
        # SINCE THE MATRIX IS SYMETRIC
        # We also want to obtain the number of nodes using the keys! 
        # For knowing how many multiplications

        val_list = list(values)
        Matrix_A = dict(val_list) # We make another dictionary with the cols
        # First we yield the A matrix to perform A (A^2) later
        yield "Matrices", val_list
        
        keys_all = Matrix_A.keys()
        nt = len(keys_all)
        c = 1   # Just for having a sense of how much time it will take.
        for i in keys_all:
            if (c % 50 == 0):
                print str(c)+" / "+ str(nt) + " getting"
            c+= 1
            for j in keys_all:
                key_aux = i + "_" + j
                yield key_aux, [Matrix_A[i], Matrix_A[j]]
 
    ######################################
    ## MAPPERS are gonna do multiplications
    ######################################
    def mapper_MatrixMult(self, key, values):
        # This mapper expects the 2 vectors to multiply.
        # The key contains the position
        # Multiply the vectors and output the value only if it is nonzero
        # The output key will be the row only.
        if (key == "Matrices"):
            yield key, list(values)  # The first time we do not do [0]  TODO
            
        else:
            val_list = list(values)
            v1 = val_list[0]
            v2 = val_list[1]
            #### Sparse multiplication of the vectors
            # Get the intersection of elements
            InterKeys = list(set(v1.keys()) & set(v2.keys()))
            
            if (len(InterKeys) > 0):  # If they have same keys in common
                total = 0
                for key_i in InterKeys:
                    total += v1[key_i] * v2[key_i]

                i = key.split("_")[0]
                j = key.split("_")[1]
                
                yield i, (j,total)  # Yield the result

    def combiner_AGeneration_2(self, key, values):
        # Just combine the cells to create the output matrix
        if (key == "Matrices"):
            yield key,dict(list(values)[0])
        else:
            yield key, list(values)
            
    def reducer_AGeneration_2(self, key, values):
        # Same here, more combinations to create the matrix
        if (key == "Matrices"):
            yield key,list(values)[0]
        else:
            # Same here
            final_list = []
            for value in values:
                final_list.extend(list(value))
            final_dict = dict(final_list)
            row_elem = [key, final_dict]
            yield "AA", row_elem  # So that the next reducer get

    ######################################
    ## JOIN MATRIX AA
    ######################################
    def reducer_AGeneration_3(self, key, values):
        if (key == "Matrices"):
            Matrix_A = list(values)[0] 
            yield "Matrices",Matrix_A
            
        elif(key == "AA"):
            val_list = list(values)
            Matrix_AA = dict(val_list) # Make another dictionary with the cols
            yield "Matrices",Matrix_AA
            
    #################################################################
    ## MULTIPLY AGAIN
    #################################################################

    def reducer_MatrixMult2(self, key, values):
        final_list = list(values)
        Matrix_A = final_list[0]
        Matrix_AA = final_list[1]
        # We are only interested in the diagonal elements so... 
        # we only calculate those. We do not need to calculate the whole matrix
        print Matrix_AA.keys()
        for i in Matrix_AA.keys():
                key_aux = i + "_" + i
                yield key_aux, [Matrix_AA[i], Matrix_A[i]]    
                
    ###################################
    #### Finally get the trace of elements
    ###################################
    def reducer_trace_rows(self, key, values):
        # Get the diagonal
        lista = list(values)[0]
        dictionary = dict(lista)
        
        if key in dictionary:  # Get the trace if it is non-zero
            yield None, dictionary[key]
        else:
            yield None, 0
            
    def reducer_trace_rows_2(self, key, values):
        # Sum t
        Trace = sum(values)
        print "Trace " + str(Trace)
        N_triangles = Trace/6
        
        print "Number of triangles " + str(N_triangles)
        yield "NT", N_triangles
      
    def steps(self):
        return [
            # Get the Adjacency Matrix, the key is the row, and it is
            # eliminated in the reduction so that  we can handle all the rows
            # in the next reduction
            MRStep(mapper=self.mapper_AGeneration,
                   combiner=self.combiner_AGeneration,
                   reducer=self.reducer_AGeneration),
             
#             Generate the vector multiplications
            MRStep(reducer=self.reducer_MatrixMult),
#            
            # Multiply the matrices and obtain the result,
            # While storing the initial matrix
            MRStep(mapper=self.mapper_MatrixMult,
                   combiner=self.combiner_AGeneration_2,
                   reducer=self.reducer_AGeneration_2),
#            
            # We join the AA matrix and the A matrix
            MRStep(reducer=self.reducer_AGeneration_3),
            MRStep(reducer=self.reducer_MatrixMult2),
#
            # We multiply them again
            MRStep(mapper=self.mapper_MatrixMult,
                   combiner=self.combiner_AGeneration_2,
                   reducer = self.reducer_trace_rows),
                   
            # Add all the traces and divide by 6
            MRStep(reducer=self.reducer_trace_rows_2),
        ]
        
if __name__ == '__main__':
    MRNodesNT.run()
