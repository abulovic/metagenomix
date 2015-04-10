import os, sys
from collections import defaultdict

import metagenomix

NO_ENTRY_NAME = 'no_db_entry_name'
NO_ENTRY_RANK = 'no_db_entry_rank'

ranks = { 'superkingdom'        : 0,
          'kingdom'             : 1,
          'subkingdom'          : 2,
          'superphylum'         : 3,
          'phylum'              : 4,
          'subphylum'           : 5,
          'superclass'          : 6,
          'class'               : 7,
          'subclass'            : 8,
          'infraclass'          : 9,
          'superorder'          : 10,
          'order'               : 11,
          'suborder'            : 12,
          'infraorder'          : 13,
          'parvorder'           : 14,
          'superfamily'         : 15,
          'family'              : 16,
          'subfamily'           : 17,
          'tribe'               : 18,
          'subtribe'            : 19,
          'genus'               : 20,
          'subgenus'            : 21,
          'species group'       : 22,
          'species subgroup'    : 23,
          'species'             : 24,
          'subspecies'          : 25,
          'varietas'            : 26,
          'forma'               : 27,
          'no rank'             : 28 }

human           = 9606
mouse           = 10090
rats            = 10114
rodents         = 9989
primates        = 9443
animalia        = 33208
green_plants    = 33090
eukaryota       = 2759

archea          = 2157
bacteria        = 2
viruses         = 10239
fungi           = 4751

euglenozoa      = 33682
alveolata       = 33630
amoebozoa       = 554915
fornicata       = 207245
parabasalia     = 5719
heterolobosea   = 5752

viroids         = 12884
stramenopiles   = 33634

blastocladiomycota = 451459 #(ne)
chytridiomycota = 4761      #(ne)
cryptomycota    = 1031332   #(da)
dikarya         = 451864    #(ne)
entomophthoromycota = 1264859  #(da)
glomeromycota   = 214504    #(ne)
microsporidia   = 6029      #(DA)
neocallimastigomycota = 451455 #(da)

other           = 28384
unclassified    = 12908
artificial      = 81077


class TaxTree ():
    ''' Loads the NCBI taxonomy tree, creates both
        parent-child and child-parent relations,
        enables parent-child relationship testing and
        finding the least common ancestor.
    '''

    def __init__ (self, parent2child_fname=None, tax_nodes_fname=None):
        ''' Locates the ncbi taxonomy file and sets the important
            taxonomy assignments (such as animalia, bacteria ecc)

            :param parent2child_fname location of the ncbi taxonomy tree file
            :param tax_nodes_fname location of the file containing taxid,
            organism name and organism rank for each taxid in the tree.
        '''

        if not parent2child_fname:
            parent2child_fname = os.path.join(metagenomix.__path__[0], 'ncbi_tax_tree')
        self.load(parent2child_fname)
        if not tax_nodes_fname:
            tax_nodes_fname = os.path.join(metagenomix.__path__[0], 'taxid2namerank')
        self.load_taxonomy_data(tax_nodes_fname)

        #--------- RELEVANT TAXONOMY ASSIGNMENTS ----------#
        self._h_set_relevant_taxonomy_assignments()
        self._h_map_taxids_to_relevant_tax_nodes()

    def load (self, parent2child_fname):
        self.parent_nodes   = self._h_get_tax_nodes(parent2child_fname)
        self.child_nodes    = self._h_populate_child_nodes()

    def load_taxonomy_data(self, tax_nodes_fname):
        '''
        Uses data access object to find organism name and
        rank of each of the tax IDs.
        For each tax ID creates a node of type TaxNode
        After invoking this method, there is nodes parameter
        of type dict(key=tax_id:int, value=node:TaxNode)
        '''
        self.nodes = {}
        total = len(self.parent_nodes)
        current = 0
        tax_nodes_file = open(tax_nodes_fname, 'r')
        readline = tax_nodes_file.readline
        while (True):
            line = readline()
            if not line: break
            (taxid, org_name, rank) = line.strip().split('|')
            node = TaxNode(org_name, rank)
            self.nodes[int(taxid)] = node
        tax_nodes_file.close()

    def get_org_name(self, taxid):
        if taxid not in self.nodes:
            return NO_ENTRY_NAME
        else:
            return self.nodes[taxid].organism_name

    def get_org_rank(self, taxid):
        if taxid not in self.nodes:
            return NO_ENTRY_RANK
        else:
            return self.nodes[taxid].rank

    def is_child (self, child_taxid, parent_taxid):
        ''' Test if child_taxid is child node of parent_taxid
            Node is not the child of itself
        '''
        # check boundary conditions
        if child_taxid == parent_taxid:
            return False
        if parent_taxid == self.root:
            return True

        tmp_parent_taxid = child_taxid
        while True:
            if not self.parent_nodes.has_key(tmp_parent_taxid):
                return False
            tmp_parent_taxid = self.parent_nodes[tmp_parent_taxid]
            if tmp_parent_taxid == self.root:
                return False
            if tmp_parent_taxid == parent_taxid:
                return True

    def find_lca (self, taxid_list):
        ''' Finds the lowest common ancestor of a list of nodes

            Args:
                taxid_list ([int]): List of tax_ids
            Returns:
                (int): tax_id of LCA
        '''
        # Check if all nodes exist (and sum up blast scores)
        for taxid in taxid_list:
            if taxid != self.root and not self.parent_nodes.has_key(taxid):
                try:
                    raise Exception ("Key error, no element with id " + str(taxid))
                except Exception, e:
                    pass

        # Filter out invalid tax_ids - those without parents
        taxid_list = filter(lambda tax_id: tax_id == self.root or self.parent_nodes.has_key(tax_id) , taxid_list)

        # Check if list is empty
        if len(taxid_list) == 0:
            try:
                raise Exception ("taxid_list is empty, cannot find LCA!")
            except Exception:
                sys.stderr.write("{0}\n".format(e))
                return 1 # Assign to root

        # each of the visited nodes remembers how many
        # child nodes traversed it
        self.num_visited = defaultdict(int)

        current_taxids  = taxid_list
        num_of_nodes    = len(current_taxids)

        # now find the lowest common ancestor
        while (True):

            parent_nodes = []
            for taxid in current_taxids:
                # root node must not add itself to parent list
                if   taxid != self.root:    parent_taxid = self.parent_nodes[taxid]
                else:                       parent_taxid = None
                # if parent exists, append him to parent list
                # duplicates ensure that every traversal will count
                if parent_taxid:            parent_nodes.append(parent_taxid)

                # Check for LCA
                self.num_visited[taxid] += 1
                if self.num_visited[taxid] == num_of_nodes:
                    self.lca_root = taxid
                    return taxid

            # refresh current nodes
            current_taxids = parent_nodes


    def get_relevant_taxid (self, tax_id):
        return self.tax2relevantTax.get(tax_id, -1)

    def get_lineage(self,tax_id):
        lineage = []
        while (True):
            if tax_id == self.root:
                break
            lineage.append(tax_id)
            tax_id = self.parent_nodes[tax_id]
        return reversed(lineage)

    def get_parent_with_rank(self, tax_id, rank):
        if tax_id not in self.nodes:
	    return -1
        parent = 0
        while (True):
            if tax_id == self.root:
                return 0
            if self.nodes[tax_id].rank == rank:
                return tax_id
            tax_id = self.parent_nodes[tax_id]

    def _h_get_tax_nodes        (self, parent2child_fname):
        '''Loads the taxonomy nodes in a dictionary
           mapping the child to parent node.
        '''
        # file format per line: child_taxid parent_taxid
        with open(parent2child_fname) as fd:
            d = dict(self._h_from_parent_child_str (line) for line in fd)
        return d

    def _h_from_parent_child_str (self, line):
        '''Loads two integers (taxids) from a line
        '''
        key, sep, value = line.strip().partition(" ")
        if key == value: self.root = int(key)
        return int(key), int(value)

    def _h_populate_child_nodes (self):
        ''' Populates child nodes from parent to child
            mapping dictionary
        '''
        child_nodes = defaultdict(list)
        for (child, parent) in self.parent_nodes.iteritems():
            child_nodes[parent].append(child)
        return child_nodes

    def _h_set_relevant_taxonomy_assignments (self):
        ''' Sets some of the more important taxonomy
            assignments which can help in checking which kingdom
            an organism belongs to.
        '''
        self.potential_hosts = [human,
                                mouse,
                                rats,
                                rodents,
                                primates,
                                animalia,
                                green_plants]

        self.microbes =        [archea,
                                bacteria,
                                viruses,
                                fungi,
                                euglenozoa,
                                alveolata,
                                amoebozoa,
                                fornicata,
                                parabasalia,
                                heterolobosea,
                                viroids,
                                stramenopiles,
                                cryptomycota,
                                entomophthoromycota,
                                microsporidia,
                                neocallimastigomycota]


    def _h_map_taxids_to_relevant_tax_nodes(self):
        host_nodes = list(self.potential_hosts)
        microbe_nodes = list(self.microbes)

        self.tax2relevantTax = {}
        for microbe_node in self.microbes:
            microbe_children = self.get_all_children(microbe_node)
            for child in microbe_children:
                self.tax2relevantTax[child] = microbe_node

        for host_node in self.potential_hosts:
            host_children = self.get_all_children(host_node)
            for child in host_children:
                self.tax2relevantTax[child] = host_node

        tagged_nodes    = self.tax2relevantTax.keys()
        all_nodes       = self.parent_nodes.keys()
        untagged_nodes  = set(all_nodes).difference(tagged_nodes)
        for node in untagged_nodes:
            self.tax2relevantTax[node] = -1

    def get_all_children(self, tax_id):
        if not self.child_nodes.has_key(tax_id):
            return []
        one_step_children = self.child_nodes[tax_id]
        all_children = []
        while (True):
            if not one_step_children:
                break
            new_one_step_children = []
            all_children.extend(one_step_children)
            for child in one_step_children:
                if self.child_nodes.has_key(child):
                    new_one_step_children.extend(self.child_nodes[child])
            one_step_children = new_one_step_children
        return all_children



class TaxNode (object):
    '''
    Taxonomy nodes hold information on relevant
    taxonomy data, which are:
    * organism name (scientific)
    * taxonomy rank
    * score (arbitrary)
    '''

    def __init__(self, organism_name, rank, score=0.):
        self.organism_name = organism_name
        self.rank = rank
        self.score = score


def main():
    tt = TaxTree()
    print tt.nodes[9606].organism_name

if __name__ == '__main__':
    main()