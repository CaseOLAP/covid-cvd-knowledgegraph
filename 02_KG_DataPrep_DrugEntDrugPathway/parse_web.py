from bs4 import BeautifulSoup
import requests

class ParseWeb():
    '''
    This class contains methods for web scraping. B/c the DrugBank API is inaccessible and finding whether an entity is a protein 
    member or a protein group is not possible from the DrugBank XML file, another approach is to implement web scraping to extract this info.
    from the HTML of the web pages corresponding to the DrugBank ID of the entities in question. This information is needed to merge entities 
    that have the same UniProt ID. The goal is to replace the entities that have a duplicate UniProt ID and are classified as a 'protein group' 
    with its members. The actions of the drugs on the protein group entities will be inherited by their members.
    '''

    def __init__(self,url,db_id):
        '''
        @param url is the URL of the web page on which web scraping will be performed
        @param db_id is the DrugBank ID of the entity whose info. will be extracted
        '''
        self.url = url
        self.db_id = db_id
        
    def getHTML(self):
        '''
        @return html is the HTML object containing the HTML of the web page
        '''
        entity_web = requests.get(self.url)
        html = BeautifulSoup(entity_web.text, 'html.parser')
        return html
    
    def getKind(self,html):
        '''
        @param html is the HTML object that is to be parsed
        @return kind is one of two strings: 'protein group' or 'protein'
        '''
        if html:
            body = html.find('body')
            details_tb = body.find('main').find('div').find_all('div')[1].\
            find('dl')
            kind = details_tb.find_all('dd')[1].text
            return kind
        return None

    def getMembers(self,html,dup_uids,uid_list,protein_dbids,dbid_dict,uid_dict):
        '''
        Summary of how to retrieve members of protein group that has a duplicate UniProt ID in data:
        1. get HTML of page DrugBank web page corresponding to protein group
        2. check UniProt ID of members. If UniProt ID is in @param dup_uids list, get the DrugBank ID that corresponds to a protein (not protein
        group)
        3. Look up DrugBank ID in @param dbid_dict to get the name of the member
        
        @param html is the HTML object that is to be parsed
        @param dup_uids is the list of Uniprot IDs that correspond to multiple entities (i.e. duplicate UniProt IDs) from the data
        @param uid_list is the list of unique UniProt IDs in the data
        @param protein_dbids is a list containing only the DrugBank IDs corresponding to proteins, not protein groups.
        @param dbid_dict is a dictionary containing DrugBank IDs from the data as keys and the corresponding UniProt IDs and names as values
        @param uid_dict is a dictionary containing non-duplicate DrugBank IDs from the data as keys and the corresponding DrugBank IDs and 
        names as values
        @return info_dict is a dictionary containing the DrugBank ID corresponding to a protein group that has a duplicate UniProt ID, its 
        name, and its info. about its members (name, DrugBank ID, UniProt ID) 
        '''
        if html:
            main_div = html.find('body').find('main').find('div')
            details_tb = main_div.find_all('div')[1].find('dl')
            name = main_div.find_all('div')[0].find('h1').text
            mlist = []
            info_dict = {}
            members = details_tb.find_all('dd')[3].find('table').\
            find('tbody').find_all('tr')
            for member in members:
                muid = member.find_all('td')[1].find('a').text
                if muid in uid_list: # if the member's UniProt ID is in the list of all UniProt IDs currently in the data
                    if muid in list(dup_uids.keys()): # if member has duplicate UniProt ID
                        dbid2check = dup_uids[muid]
                        for dbid in dbid2check:
                            if dbid in protein_dbids:
                                the_dbid = dbid
                                mdict = {}
                                mdict.update({'name':dbid_dict[the_dbid]['name'],\
                                                'drugbank_id':the_dbid,\
                                                'uniprot_id':muid
                                                })
                                mlist.append(mdict)
                    else: # member has non-duplicate UniProt ID
                        mdict = {}
                        mdict.update({'name':uid_dict[muid]['name'],\
                                        'drugbank_id':uid_dict[muid]['drugbank_id'],\
                                        'uniprot_id':muid
                                        })
                        mlist.append(mdict)
                    info_dict.update({'name':name,'drugbank_id':self.db_id,\
                                    'members':mlist})
            return info_dict
        return None