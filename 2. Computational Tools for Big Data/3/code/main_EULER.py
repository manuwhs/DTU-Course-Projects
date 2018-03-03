from mrjob.job import MRJob
from mrjob.step import MRStep
# This Job is meant to determine if a graph has an Euler tour 
# This happends if all vertices have even degree, or only 2 vertex
# have an even degree

class MRWordVC(MRJob):

    def mapper_VC(self, _, line):
        # This mapper assigns per every vertex, a 1 to every node in the vertex.
        nodes_vertex = line.split()
        if (len(nodes_vertex) == 2):
            yield nodes_vertex[0], 1
            yield nodes_vertex[1], 1
        
    def combiner_VC(self, key, values):
        # Each combiner will get a set of (key,value) pairs 
        # And it will join the ones it can.
        yield key, sum(values)

    def reducer_VC(self, key, values):
        # The reducer will do the same joining with all the joining
        # from the combiner. And it will also say if the number of 
        # vertex is odd or even for every node
        yield None, sum(values)%2

    def reducer_oddVertex(self, _, vertex_count_nodes):
        # This reducer gets all the counted words and yields
        # then in order of frequency
        suma_odd = sum(vertex_count_nodes)
        yield None, "Number of odd vertex nodes " + str(suma_odd)
        
        if (suma_odd == 0 or suma_odd == 2):
            yield 1, "There is an Euler Path"
        else:
            yield 0, "There is no Euler Path"
            
    def steps(self):
        return [
            MRStep(mapper=self.mapper_VC,
                   combiner=self.combiner_VC,
                   reducer=self.reducer_VC),
            MRStep(reducer=self.reducer_oddVertex)
        ]
        
if __name__ == '__main__':
    MRWordVC.run()
