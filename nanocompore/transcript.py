class Transcript():
    def __init__(self, id, ref_id, ref_seq):
        self._id = id
        self._name = ref_id
        self._ref_seq = ref_seq
        self._length = len(self._ref_seq)


    @property
    def id(self):
        return self._id


    @property
    def name(self):
        return self._name


    @property
    def length(self):
        return self._length


    @property
    def seq(self):
        return self._ref_seq

