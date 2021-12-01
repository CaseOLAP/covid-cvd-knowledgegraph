import copy

class CleanseData():
    '''
    This class contains methods for cleansing the data extracted from DrugBank containing info. about cardiovascular drugs. 
    There are 3 main parts of the data cleansing process:
    1. Remove drugs that do not interact with any entity or those that only interact with entities without a uniprot ID from list of CV drugs.
    This is implemented via the removeDrugs method.
    2. Remove proteins with a null uniprot ID. This is implemented via the removeNullUIDEntity method.
    3. Merge entities with the same UniProt ID. More specifically, this is done by replacing entity that is classified as 'protein group' on
    the corresponding DrugBank web page (the ParseWeb class is implemented to do web scraping to get  this info) with members of that protein
    group. The association between the protein group entities and their drug(s) are established via non-specific publication findings. In other
    words, the publication in question finds that the drug is related to the protein group, but not any of its members specifically. The 
    UniProt ID for the protein group is the same as its first member, ordered alpha-numerically. 
    
    For example, in our data, the entities 'Beta adrenergic receptor' and 'Beta-1 adrenergic receptor' both have the same UniProt ID (P08588).
    'Beta adrenergic receptor' is the protein group, while 'Beta-1 adrenergic receptor' is one of its members. The former is assigned the same
    UniProt ID as the latter rather than its other members (Beta-2 adrenergic receptor and Beta-3 adrenergic receptor) b/c Beta-1 adrenergic 
    receptor is the first member when all of them are ordered alpha-numerically. The associations between 'Beta adrenergic receptor' and its
    drugs are established based on non-specific publication findings. These findings show that the drugs are related to this protein group, but
    not any of its members specifically.
    '''

    def __init__(self,DATA):
        '''
        @param DATA is the data extracted from DrugBank containing info. about cardiovascular drugs. This data will be used to initialize
        the object. 
        '''
        self.data = DATA

    def getDrugIndexList(self):
        '''
        @return list(set(index_list)) is a list of index of drugs in the data that do not interact with with any entity (i.e. the value of 
        'targets', 'carriers', 'enzymes', and 'transporters' is an empty list). These will be removed.
        '''
        index_list = []
        for i,drug in enumerate(self.data):
            if drug['targets'] == [] and \
            drug['enzymes'] == [] and \
            drug['carriers'] == [] and \
            drug['transporters'] == []:
                index_list.append(i)
        return list(set(index_list))

    def getNullUIDIndexList(self):
        '''
        @return indexl is a list of index of drugs in the data that contain at least one entity with a null UniProt ID.
        '''
        indexl = []
        for i,drug in enumerate(self.data):
            for entity in ['targets','enzymes','carriers','transporters']:
                for j in drug[entity]:
                    if j['uniprot_id']=="Null":
                        indexl.append(i)
        indexl = list(set(indexl))
        return indexl

    def getNullUIDIndexInfo(self):
        '''
        This method calls the getNullUIDIndexList method.
        @return nullid_index_info is a list of dictionaries, each containing an index in the list returned by the getNullUIDIndexList method 
        as the key and a list of the names of ALL the entities (not just the entity with the null UniProt ID) for that particular index as
        the values.
        '''
        indexl = self.getNullUIDIndexList()
        nullid_index_info = []
        for i in indexl:
            names = []
            index_dict = {}
            for entity in ['targets','enzymes','carriers','transporters']:
                for j in self.data[i][entity]:    
                    names.append(j['name'])
            index_dict.update({i:names})
            nullid_index_info.append(index_dict)
        return nullid_index_info

    def getDupUIDs(self):
        '''
        @return dup_uids is a dictionary, with the keys being Uniprot IDs that correspond to multiple entities (i.e. duplicate UniProt IDs) from 
        the data. The values of the dictionary are DrugBank IDs of entities with the duplicated UniProt ID. 
        '''
        uid_list = []
        ent_dict = {}
        for drug in self.data:
            for entity in ['targets','enzymes','transporters','carriers']:
                for ent in drug[entity]:
                    uid = ent['uniprot_id']
                    if not uid in uid_list:
                        uid_list.append(uid)
                        ent_dict[uid] = [ent['drugbank_id']]
                    else:
                        if not ent['drugbank_id'] in ent_dict[uid]:
                            ent_dict[uid].append(ent['drugbank_id'])
        dup_uids = {}
        for item in list(ent_dict.items()):
            if len(item[1])>1:
                dup_uids[item[0]]=item[1]
        return dup_uids

    def removeNullUIDEntity(self):
        '''
        This method calls the getNullUIDIndexList method.
        @return self.data is the data after removing entities with a null UniProt ID.
        '''
        indexl = self.getNullUIDIndexList()
        for index in indexl:
            for entity in ['targets','enzymes','carriers','transporters']:
                for i,j in enumerate(self.data[index][entity]):
                    if j['uniprot_id']=='Null':
                        del self.data[index][entity][i]
        return self.data

    def removeDrugs(self):
        '''
        This method calls the getNullUIDIndexInfo and getDrugIndexList methods.
        @return self.data is the data after removing drugs (elements) in the data that interact exclusively with an entity that has a null 
        UniProt ID.
        '''
        nullid_index_info = self.getNullUIDIndexInfo()
        index_list = self.getDrugIndexList()
        for i in nullid_index_info:
            if len(list(i.values())[0])==1:
                index_list.append((list(i.keys())[0]))
        index_list.sort()
        for index in reversed(index_list):
            del self.data[index]
        return self.data

    def mergeDuplicateUIDs(self,ent_dbid,ent_list):
        '''
        @ent_dbid is the list containing the DrugBank IDs of the protein groups that have a duplicate UniProt ID
        @ent_list is the list containing dictionaries with each correspinding to protein groups that have a duplicate UniProt ID
        @return self.data is the data after replacing the entities with a duplicate UniProt ID and are classified as a 'protein group' with
        its members. Each member inherits the actions of the protein group. Because it is not known what the action of the drug is on each 
        specific member, instead of 'actions', this is represented by a new key called 'actions_of_group'. The name of the protein group
        from which each member inherits from is indicated by a new key called 'group_name'.
        '''
        for drug in self.data:
            for entity in ['targets','enzymes','carriers','transporters']:
                drug_ent = [i.get('drugbank_id') for i in drug[entity]]
                ent2del = [i for i in drug[entity] if i.get('drugbank_id') in ent_dbid]
                for ent in ent2del:
                    dbid = ent['drugbank_id']
                    actions = ent['actions']
                    group_name = ent['name']
                    ent_list_ele = [i for i in ent_list if i.get('drugbank_id')==dbid][0]
                    members = [i for i in ent_list_ele['members'] if i.get('drugbank_id') not in drug_ent]
                    for member in members:
                        tmp = copy.deepcopy(member)
                        tmp.update({'group_name':group_name,\
                            'actions_of_group':actions})
                        drug[entity].append(tmp)
                drug[entity] = [i for i in drug[entity] if i.get('drugbank_id') not in ent_dbid]
        return self.data