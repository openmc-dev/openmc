class MsrContinuous:

    def __init__(self, model):
        self.model = model
        self.geometry = model.geometry
        self.removal_terms = []

    def add_transfer(self, elements, Lambda, dest_mat):
        transfer = OrderedDict()
        transfer['element'] = elements
        transfer['lambda'] = Lambda
        transfer['dest_mat'] = dest_mat
        return transfer

    def add_removal_term(self, mat, transfers):
        removal_term = OrderedDict()
        removal_term['mat'] = mat
        for transfer in transfers:
            removal_term['transfer'].append(transfer)
        self.removal_terms.append(removal_term)
