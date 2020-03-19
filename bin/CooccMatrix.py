import numpy as np
import scipy.sparse as sps

class CooccMatrix(object):
    
    def __init__(self):
        
        self.win = None
        
    def get_window(self, tokens):
        for x,i in enumerate(tokens):
            start = np.max([0,x-(self.win+1)])
            end = np.min([len(tokens), x+(self.win+1)])
            yield((tokens[x], np.array(tokens[start:end])))
           
    def set_item_map(self, item_map):
        self.item_map = item_map
        

    def build_coocc_matrix(self, data, win=None, threshold = None):
        self.win = win
        
        item_set = set(self.item_map.keys())
        n_items = len(self.item_map)
            
        coocc_matrix = sps.dok_matrix((n_items, n_items))
        
        for info in data:
            subVec = np.zeros(n_items)
            indices = []
            for item in info:
                ind = self.item_map.get(item)
                indices.append(ind)
            
            if self.win:
                for center_word,context_words in self.get_window(indices):
                    for cw in context_words:
                        coocc_matrix[center_word,cw] += 1
                    
            else:
                locs = np.meshgrid(indices, indices)
                coocc_matrix[locs[0], locs[1]] += 1
        
        self.coocc_matrix = coocc_matrix
    
    def get_coocc_dict(self,threshold=0):
        outDict = dict()
        for k,v in self.coocc_matrix.items():
            if k[0] not in outDict:
                outDict[k[0]] = dict()
            if k[1] not in outDict:
                outDict[k[1]] = dict()
            if v >= threshold:
                outDict[k[0]][k[1]] = v
                outDict[k[1]][k[0]] = v
                
        return(outDict) 
     
