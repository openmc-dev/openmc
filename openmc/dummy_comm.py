import sys


class DummyCommunicator:
    rank = 0
    size = 1

    def allgather(self, sendobj):
        return [sendobj]

    def allreduce(self, sendobj, op=None):
        return sendobj

    def barrier(self):
        pass

    def bcast(self, obj, root=0):
        return obj

    def gather(self, sendobj, root=0):
        return [sendobj]

    def py2f(self):
        return 0

    def reduce(self, sendobj, op=None, root=0):
        return sendobj

    def scatter(self, sendobj, root=0):
        return sendobj[0]

    def Abort(self, exit_code_or_msg):
        sys.exit(exit_code_or_msg)
