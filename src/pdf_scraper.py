import requests
from bs4 import BeautifulSoup
import openpyxl
import time
from scholarly import scholarly

def pubmed_search_ids(query):
    '''
    Parameters:
    ----------------
    query: the query of the pubmed search

    Returns:
    ----------------
    ids_list: list of ids that query search returns
    '''

    id_list = []
    # the url for PubMed query search
    search_ids_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?term="

    # the response from the url
    ids_response = requests.get(search_ids_url + query)

    # sleep for 0.4 seconds to not get blocked
    time.sleep(0.4)

    # the beautiful soup of the url
    ids_soup = BeautifulSoup(ids_response.text, features="xml")

    # the list of ids in xml format
    ids_xml_list = ids_soup.find_all('Id')

    for id_tag in ids_xml_list:
        if id_tag is not None:  
            # Extract the text content of the Id tag
            id_num = int(id_tag.text)
            id_list.append(id_num)

    return id_list  # Return id_list instead of ids_list

path = "../data/papers.xlsx"
def extract_articles(path):
    '''
    Parameters:
    ----------------
    path: the path of the excel search

    Returns:
    ----------------
    articles: dictionary of articles where the key is the sheet name and the value is a list of articles
    '''
    articles = []

    # Define variable to load the dataframe
    dataframe = openpyxl.load_workbook(path)

    # iterates through each sheet and adds each title into list
    for i in range(len(dataframe.sheetnames)):
        sheetname = dataframe.sheetnames[i]
        curr_sheet = dataframe._sheets[i]
        for i in range(1, curr_sheet.max_row):
            title = curr_sheet.cell(row = i, column = 8).value
            authors = curr_sheet.cell(row = i, column = 5).value
            if title is not None and authors is not None and title != 'TITLE' and authors != 'AUTHORS':
                pair = [title, authors]
                articles.append(pair)

    return articles

def extract_authors(path):
    '''
    Parameters:
    ----------------
    path: the path of the excel search

    Returns:
    ----------------
    authors: dictionary of titles where the key is the sheet name and the value is a list of titles of articles
    '''
    authors = {}

    # Define variable to load the dataframe
    dataframe = openpyxl.load_workbook(path)

    # iterates through each sheet and adds each title into list
    for i in range(len(dataframe.sheetnames)):
        sheetname = dataframe.sheetnames[i]
        curr_sheet = dataframe._sheets[i]
        authors[sheetname] = []
        for i in range(1, curr_sheet.max_row):
            title = curr_sheet.cell(row = i, column = 4).value
            if title is not None:
                titles[sheetname].append(title)
        
        # Gets rid of sheets without any titles
        if authors[sheetname] == []:
            del authors[sheetname]

    return authors


def pubmed_get_abstracts(id):
    '''
    Parameters:
    ----------------
    id: the id of the article

    Returns:
    ----------------
    abstract: the abstract of the article
    '''
    abstract = ""

    # the url of the article
    url = "https://pubmed.ncbi.nlm.nih.gov/" + str(id)

    # the response from the url
    response = requests.get(url)

    # sleep for 0.4 seconds to not get blocked
    time.sleep(0.4)

    # the beautiful soup of the url
    soup = BeautifulSoup(response.text, features="lxml")
    
    try:
        # finds the abstract content
        abstract_soup = soup.find("div", {"class": "abstract-content selected"})
        paragraph_soup = abstract_soup.find_all("p")

        for p in paragraph_soup:
            abstract = abstract + " " + p.get_text(strip=True)

        return abstract
    except:
        return None

def pubmed_get_free_fulltext_link(id):
    '''
    Parameters:
    ----------------
    id: the id of the article

    Returns:
    ----------------
    link: the link to the free full text of the article
    '''

    # the url of the article
    url = "https://pubmed.ncbi.nlm.nih.gov/" + str(id)

    # the response from the url
    response = requests.get(url)

    # sleep for 0.4 seconds to not get blocked
    time.sleep(0.4)

    # the beautiful soup of the url
    soup = BeautifulSoup(response.text, features="lxml")
    
    # finds the link
    link_soup = soup.find("a", {"class": "link-item pmc"})
    
    try:
        # finds the link
        link = link_soup['href']

    except:
        # no link
        return None

def get_full_text_from_search_results(search_results):
    batch_size = 100
    end = min(int(search_results["Count"]), start + batch_size)
    handle = Entrez.efetch(db="pmc", rettype="full", retmode="xml", 
                        retstart=start, retmax=batch_size,
                        webenv=search_results["WebEnv"], 
                        query_key=search_results["QueryKey"])
    data = handle.read()
    soup = BeautifulSoup(data, features="xml")
    body_soup = soup.find("body")
    p_soup = body_soup.find_all("p")
    for p in p_soup:
        full_text = full_text + p.text
    # fulltexts = open("../data/fulltext_test.txt", "w")
    # fulltexts.write(str(data))# + ",\n")
    # fulltexts.close()
    handle.close()
    return full_text



# ''' MAIN (abstracts)'''
# # getting the titles from the excel sheet
# path = "../data/papers.xlsx"
# titles = extract_titles(path)

# abstracts_downloaded = 0
# abstracts_not_downloaded = 0

# for title in titles.keys():
#     values = titles[title]

#     for values in titles.values():
#         for title in values:
#             id_list = pubmed_search_ids(title)
#             if len(id_list) == 1: # checking to see if only one query came back
#                 abstract = pubmed_get_abstracts(id_list[0])
#                 if abstract is not None:
#                     ## add the abstracts to a database

#                     try:
#                         abstracts = open("../data/abstracts.txt", "a")
#                         abstracts.write(abstract + ",\n")
#                         abstracts.close()

#                         abstracts_downloaded = abstracts_downloaded + 1

#                         if abstracts_downloaded < 10 or (abstracts_downloaded % 10 and abstracts_downloaded < 100) or (abstracts_downloaded % 100 and abstracts_downloaded < 1000) or (abstracts_downloaded % 500 and abstracts_downloaded < 10000) or abstracts_downloaded % 1000:
#                             print(str(abstracts_downloaded) + " abstracts downloaded.")
#                     except:
#                         abstracts_not_downloaded = abstracts_not_downloaded + 1
#                         print(str(abstracts_not_downloaded) + " abstracts not added.")

# print("Abstracts downloaded: " + str(abstracts_downloaded))
# print("Abstracts not downloaded: " + str(abstracts_not_downloaded))



''' MAIN (full texts) '''

from Bio import Entrez

Entrez.email = "eric7lee87@gmail.com" 


success = False


# getting the titles from the excel sheet
path = "../data/papers.xlsx"
articles = extract_articles(path)

articles_accessed = 0

for article in articles:
    articles_accessed += 1
    print("Articles accessed: " + str(articles_accessed))
    title = article[0]
    authors = article[1]
    authors_list = [line.strip() for line in authors.strip().split('\n\n')]
    author_query = authors_list[0]

    try:
        search_query = str(title) + '[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
        search_results = Entrez.read(Entrez.esearch(db="pubmed", term=search_query, retmax=10, usehistory="y"))
        time.sleep(1)
        success = False

        if int(search_results["Count"]) == 1: ## if searching only for title worked
            success = True

        elif int(search_results["Count"]) > 1: ## if too many queries popped up
            search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with author
            search_results = Entrez.read(Entrez.esearch(db="pubmed", term=search_query, retmax=10, usehistory="y"))
            time.sleep(1)
            if int(search_results["Count"]) > 1: ## searching with author did not work
                for i in range(1, len(authors_list)):
                    author_query = " OR "+ author_query
                search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with authors
                search_results = Entrez.read(Entrez.esearch(db="pubmed", term=search_query, retmax=10, usehistory="y"))
                time.sleep(1)
                if int(search_results["Count"]) > 1: ## searching with authors did not work
                    search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author] AND medline[sb] AND "open access"[filter]' ## try searching with filters
                    search_results = Entrez.read(Entrez.esearch(db="pubmed", term=search_query, retmax=10, usehistory="y"))
                    time.sleep(1)
                    if int(search_results["Count"]) > 1: ## searching with filter did not work
                        pass
                    elif int(search_results["Count"]) == 1: ## searching with filters worked
                        success = True
                else: ## searching with authors did work
                    success = True
            else: ## searching with author worked
                success = True

        elif int(search_results["Count"]) == 0: ## no search results
            pass
        
        if success:
            fulltext = get_full_text_from_search_results(search_results)
            print("Success!")
            fulltexts = open("../data/fulltexts.txt", "a")
            fulltexts.write(fulltext + ",\n")
            fulltexts.close()

    except Exception as e:
        continue




