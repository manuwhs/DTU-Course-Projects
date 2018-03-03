from mrjob.job import MRJob
from mrjob.step import MRStep

class MRWordWC(MRJob):

    def mapper_WC(self, _, line):
        # Each mapper will create the (key,value) pairs of words
        # for every line they get.
        words = line.split()
        for word in words:
            yield word.lower(), 1
            
    def combiner_WC(self, key, values):
        # Each combiner will get a set of (key,value) pairs 
        # And it will join the ones it can.
        yield key, sum(values)

    def reducer_WC(self, key, values):
        # The reducer will do the same joining with all the joining
        # from the combiner
    
        # Instead of yielding key-value pairs with yield key, sum(values)
        # We just yield values, with the words and their count.
        # This way, the next reducer will get all values at once
#        print list(key)
#        print list(values)
        yield None, (sum(values), key)

    def reducer_orderedWords(self, _, word_count_pairs):
        # This reducer gets all the counted words and yields
        # then in order of frequency
        def getKey(item):
            return item[0]
            
        ordered = sorted(word_count_pairs, key= getKey, reverse=True)
        yield None, ordered

    def steps(self):
        return [
            MRStep(mapper=self.mapper_WC,
                   combiner=self.combiner_WC,
                   reducer=self.reducer_WC),
            MRStep(reducer=self.reducer_orderedWords)
        ]
        
if __name__ == '__main__':
    MRWordWC.run()
