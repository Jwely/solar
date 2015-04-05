

class chunk_bundle():
    """
    Creates a chunk bundle object.

    it can be used to pass smaller pieces of information to complex functions
    and reduce memory consumption in those functions. It can drastically
    decrease performance, since it forces writing data to hard disk, so do not
    use this unless you would otherwise run out of memory.
    """

    def __init__(self, name, numpy_array, metadata, num_chunks, temp_dir):

        self.name = name            # name of object, used for saving chunks
        self.data = numpy_array     # input raster numpy array
        self.meta = metadata        # metadata object
        self.subs = num_chunks      # number of chunks to split raster into
        self.work = temp_dir        # temporary workspace to save chunks

        return

    def split(self, num_chunks): pass

    def pull(self, chunk_id):
        """ reads a specific chunk and returns the info"""
        pass

    def push(self, chunk_id):
        """ replaces a chunk with the input chunk data"""
        pass

    def save(self): pass


class chunk():
    """ creates an individual chunk object. chunk_bundles are made of these"""

    def __init__(self, chunk_id, numpy_array = None, filepath = None)


        self.name = chunk_id
        
        if numpy_array:
            self.data = numpy_array

        if filepath:
            self.path = filepath

            
    
