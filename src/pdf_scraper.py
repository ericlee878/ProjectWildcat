import os
import requests
from bs4 import BeautifulSoup
import openpyxl
import time
import json
from fpdf import FPDF
from Bio import Entrez

Entrez.email = "eric7lee87@gmail.com" 

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
def extract_metadata(path):
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
            title = curr_sheet.cell(row = i, column = 8).value #H
            authors = curr_sheet.cell(row = i, column = 5).value #E
            pmid = curr_sheet.cell(row = i, column = 3).value #C
            date = curr_sheet.cell(row = i, column = 13).value #M
            abstract = curr_sheet.cell(row = i, column = 11).value #K
            mesh = curr_sheet.cell(row = i, column = 14).value #N
            keywords = curr_sheet.cell(row = i, column = 17).value #Q
            if mesh is not None:
                mesh = mesh.split()
            if keywords is not None:
                keywords = keywords.split()
            if title is not None:
                article = [title, authors, pmid, date, abstract, mesh, keywords]
                articles.append(article)

    return articles

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


def get_full_text_from_search_results(search_results):
    batch_size = 100
    for start in range(0, int(search_results["Count"]), batch_size):
        end = min(int(search_results["Count"]), start + batch_size)
        handle = Entrez.efetch(db="pmc", rettype="full", retmode="xml", 
                            retstart=start, retmax=batch_size,
                            webenv=search_results["WebEnv"], 
                            query_key=search_results["QueryKey"])
        data = handle.read()
        soup = BeautifulSoup(data, features="xml")
        # body_soup = soup.find("body")
        p_soup = soup.find_all("p")
        full_text = ""
        for p in p_soup:
            full_text = full_text + p.text
        return full_text


def create_pdf(article_dict, output_folder):
    filename = f"{article_dict['pmid']}.pdf"
    filepath = os.path.join(output_folder, filename)
    
    pdf = FPDF()
    pdf.add_page()
    
    # Set font for title
    pdf.set_font("Arial", 'B', size=16)
    
    # Add title
    pdf.multi_cell(0, 10, txt=article_dict['title'], align='C')
    pdf.ln(10)
    
    # Set font for normal text
    pdf.set_font("Arial", size=12)
    
    # Add authors
    pdf.multi_cell(0, 10, txt=f"Authors: {article_dict['authors']}", align='L')
    pdf.ln(5)
    
    # Add publication date
    pdf.multi_cell(0, 10, txt=f"Publication Date: {article_dict['publication_date']}", align='L')
    pdf.ln(10)
    
    # Add abstract
    if article_dict['abstract'] is not None:
        pdf.set_font("Arial", 'B', size=14)
        pdf.cell(0, 10, txt="Abstract:", ln=1)
        pdf.set_font("Arial", size=12)
        pdf.multi_cell(0, 10, txt=article_dict['abstract'])
        pdf.ln(10)
    else:
        pdf.set_font("Arial", 'B', size=14)
        pdf.cell(0, 10, txt="Abstract:", ln=1)
        pdf.set_font("Arial", size=12)
        pdf.multi_cell(0, 10, txt='Not Available')
        pdf.ln(10)
    
    # # Add full text if available
    # if article_dict['has_full_text'] == "yes" and article_dict['has_full_text'] is not None:
    #     pdf.add_page()
    #     pdf.set_font("Arial", 'B', size=14)
    #     pdf.cell(0, 10, txt="Full Text:", ln=1)
    #     pdf.set_font("Arial", size=12)
    #     pdf.multi_cell(0, 10, txt=article_dict['full_text'])
    # else:
    #     pdf.add_page()
    #     pdf.set_font("Arial", 'B', size=14)
    #     pdf.cell(0, 10, txt="Full Text:", ln=1)
    #     pdf.set_font("Arial", size=12)
    #     pdf.multi_cell(0, 10, txt='Not Available')
    

    pdf.output(filepath)
    print(f"Created PDF for article {article_dict['pmid']}")


def create_large_pdf(dict_list, output_folder, filename):
    filepath = os.path.join(output_folder, filename)
    
    pdf = FPDF()
    pdf.add_page()
    
    # Set font for Filename
    pdf.set_font("Arial", 'B', size=30)

    # Add Filename
    pdf.multi_cell(0, 10, txt=filename.capitalize(), align='C')
    pdf.ln(10)

    pdf_num = 0

    for d in dict_list:
        pdf.add_page()
        # Set font for title
        pdf.set_font("Arial", 'B', size=16)
        
        # Add title
        pdf.multi_cell(0, 10, txt=d['title'], align='C')
        pdf.ln(10)
        
        # Set font for normal text
        pdf.set_font("Arial", size=12)
        
        # Add authors
        pdf.multi_cell(0, 10, txt=f"Authors: {d['authors']}", align='L')
        pdf.ln(5)
        
        # Add publication date
        pdf.multi_cell(0, 10, txt=f"Publication Date: {d['publication_date']}", align='L')
        pdf.ln(10)
        
        # Add abstract
        if d['abstract'] is not None:
            pdf.set_font("Arial", 'B', size=14)
            pdf.cell(0, 10, txt="Abstract:", ln=1)
            pdf.set_font("Arial", size=12)
            pdf.multi_cell(0, 10, txt=d['abstract'])
            pdf.ln(10)
        else:
            pdf.set_font("Arial", 'B', size=14)
            pdf.cell(0, 10, txt="Abstract:", ln=1)
            pdf.set_font("Arial", size=12)
            pdf.multi_cell(0, 10, txt='Not Available')
            pdf.ln(10)

        if pdf.page_no() >= 198:
            print("New PDF page")
            pdf.output(filepath + "_" + str(pdf_num) + ".pdf")
            pdf = FPDF()
            pdf.add_page()
            pdf_num = pdf_num + 1
        
        # # Add full text if available
        # if d['has_full_text'] == "yes" and d['has_full_text'] is not None:
        #     pdf.add_page()
        #     pdf.set_font("Arial", 'B', size=14)
        #     pdf.cell(0, 10, txt="Full Text:", ln=1)
        #     pdf.set_font("Arial", size=12)
        #     pdf.multi_cell(0, 10, txt=d['full_text'])
        # else:
        #     pdf.add_page()
        #     pdf.set_font("Arial", 'B', size=14)
        #     pdf.cell(0, 10, txt="Full Text:", ln=1)
        #     pdf.set_font("Arial", size=12)
        #     pdf.multi_cell(0, 10, txt='Not Available')
        
    pdf.output(filepath + "_" + str(pdf_num) + ".pdf")


''' MAIN (full texts) '''

# from Bio import Entrez

# Entrez.email = "eric7lee87@gmail.com" 


# success = False


# # getting the titles from the excel sheet
# path = "../data/papers.xlsx"
# articles = extract_metadata(path)

# import os

# json_file_path = "../data/abstracts.json"

# if not os.path.exists(json_file_path):
#     with open(json_file_path, 'w') as f:
#         json.dump([], f)

# # Create output folder for PDFs
# output_folder = "..\\data\\full_texts_and_abstracts_pdf"

# if not os.path.exists(output_folder):
#     os.makedirs(output_folder)

# articles_accessed = 0

# for article in articles:
#     try:
#         articles_accessed += 1
#         title = article[0]
#         authors = article[1]
#         pmid = article[2]
#         date = article[3]
#         abstract = article[4]
#         authors_list = [line.strip() for line in authors.strip().split('\n\n')]
#         author_query = authors_list[0]

#         if title == 'TITLE':
#             continue

#         search_query = str(title) + '[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
#         search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#         time.sleep(1)
#         success = False

#         if int(search_results["Count"]) == 1: ## if searching only for title worked
#             success = True

#         elif int(search_results["Count"]) > 1: ## if too many queries popped up
#             search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with author
#             search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#             time.sleep(1)
#             if int(search_results["Count"]) > 1: ## searching with author did not work
#                 for i in range(1, len(authors_list)):
#                     author_query = " OR "+ author_query
#                 search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with authors
#                 search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#                 time.sleep(1)
#                 if int(search_results["Count"]) > 1: ## searching with authors did not work
#                     search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author] AND medline[sb] AND "open access"[filter]' ## try searching with filters
#                     search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#                     time.sleep(1)
#                     if int(search_results["Count"]) > 1: ## searching with filter did not work
#                         success = False
#                     elif int(search_results["Count"]) == 1: ## searching with filters worked
#                         success = True
#                 else: ## searching with authors did work
#                     success = True
#             else: ## searching with author worked
#                 success = True

#         elif int(search_results["Count"]) == 0: ## no search results
#             success = False

#         ## Finding abstracts
#         if abstract == "" or abstract is None:
#             id_list = pubmed_search_ids(title)
#             if len(id_list) == 1: # checking to see if only one query came back
#                 abstract = pubmed_get_abstracts(id_list[0])

#         if success and get_full_text_from_search_results(search_results) is not None:
#             fulltext = get_full_text_from_search_results(search_results)
#             # Create article dictionary
#             article_dict = {
#                 "pmid": pmid,
#                 "title": title,
#                 "authors": ", ".join(authors_list),
#                 "has_full_text": "yes",
#                 "abstract": abstract,
#                 "full_text": fulltext,
#                 "publication_date": date
#             }
#             # Read existing data
#             with open(json_file_path, 'r') as f:
#                 data = json.load(f)

#             # Append new article
#             data.append(article_dict)

#             # Write updated data back to file
#             with open(json_file_path, 'w') as f:
#                 json.dump(data, f, indent=2)

#             create_pdf(article_dict, output_folder)
#             print(f"Added article {articles_accessed} to JSON file. Success!")
#         else:
#             # Create article dictionary
#             article_dict = {
#                 "pmid": pmid,
#                 "title": title,
#                 "authors": ", ".join(authors_list),
#                 "has_full_text": "no",
#                 "abstract": abstract,
#                 "full_text": "",
#                 "publication_date": date
#             }
#             # Read existing data
#             with open(json_file_path, 'r') as f:
#                 data = json.load(f)

#             # Append new article
#             data.append(article_dict)

#             # Write updated data back to file
#             with open(json_file_path, 'w') as f:
#                 json.dump(data, f, indent=2)

#             create_pdf(article_dict, output_folder)
#             print(f"Added article {articles_accessed} to JSON file. No full text :(")
#     except Exception as e:
#         print("Error: " + str(e))





############################################## FULL TEXT ACCESSOR TESTING #################################################
# from Bio import Entrez

# Entrez.email = "eric7lee87@gmail.com" 

# search_query = 'Fibrin Sealants in Dura Sealing: A Systematic Literature Review' + '[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
# search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
# batch_size = 100
# for start in range(0, int(search_results["Count"]), batch_size):
#     end = min(int(search_results["Count"]), start + batch_size)
#     handle = Entrez.efetch(db="pmc", rettype="full", retmode="xml", 
#                         retstart=start, retmax=batch_size,
#                         webenv=search_results["WebEnv"], 
#                         query_key=search_results["QueryKey"])
#     data = handle.read()
#     soup = BeautifulSoup(data, features="xml")
#     # body_soup = soup.find("body")
#     p_soup = soup.find_all("p")
#     full_text = ""
#     for p in p_soup:
#         full_text = full_text + p.text
#     print(full_text)
#     # handle.close()

############ CODE THAT CREATES KEYWORD DICTIONARY ###########################

# articles_accessed = 0
# articles = extract_metadata(path)
# keyword_dict = {}
# for article in articles:
#     articles_accessed += 1
#     pmid = article[2]
#     mesh = article[5]
#     keywords = article[6]

#     if mesh is not None and keywords is not None:
#         combined = mesh + keywords
#         combined = list(set(combined))
#     elif mesh is None and keywords is None:
#         continue
#     elif mesh is None and keywords is not None:
#         combined = keywords
#     elif keywords is None and mesh is not None:
#         combined = mesh

#     if combined is not None:
#         for word in combined:
#             word = word.lower()
#             for letter in word:
#                 if letter in '*)]";,':
#                     word = word.replace(letter,'')
#                 elif letter in "[(/":
#                     if word.index(letter) == 0:
#                         word = word.replace(letter,'')
#                     else:
#                         word = word[0:word.index(letter)]
#             if word in keyword_dict.keys():
#                 if pmid not in keyword_dict[word]:
#                     keyword_dict[word].append(pmid)
#             else:
#                 keyword_dict[word] = [pmid]

# json_file_path = "../data/dictionary.json"

# if not os.path.exists(json_file_path):
#     with open(json_file_path, 'w') as f:
#         json.dump([], f)

# with open(json_file_path, 'w') as f:
#     json.dump(keyword_dict, f, indent=2)

############### CODE THAT DOWNLOADS PDFS ################

# # getting the titles from the excel sheet
# path = "../data/papers.xlsx"
# articles = extract_metadata(path)

# word = input("Input your word: ")

# output_folder = "..\\pdfs\\" + word
# if not os.path.exists(output_folder):
#     os.makedirs(output_folder)      
      
# dict_list = []

# articles_found = 0

# ## getting the relevant PMIDs from json

# file = open('../data/dictionary.json')
# data = json.load(file)
# pmids = data[word.lower()]
# for article in articles:
#     try:
#         title = article[0]
#         authors = article[1]
#         pmid = article[2]
#         date = article[3]
#         abstract = article[4]
#         authors_list = [line.strip() for line in authors.strip().split('\n\n')]
#         author_query = authors_list[0]
#         mesh = article[5]
#         keywords = article[6]

#         if pmid in pmids:

#             if title == 'TITLE':
#                 continue

#             search_query = str(title) + '[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
#             search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#             time.sleep(1)
#             success = False

#             if int(search_results["Count"]) == 1: ## if searching only for title worked
#                 success = True

#             elif int(search_results["Count"]) > 1: ## if too many queries popped up
#                 search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with author
#                 search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#                 time.sleep(1)
#                 if int(search_results["Count"]) > 1: ## searching with author did not work
#                     for i in range(1, len(authors_list)):
#                         author_query = " OR "+ author_query
#                     search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with authors
#                     search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#                     time.sleep(1)
#                     if int(search_results["Count"]) > 1: ## searching with authors did not work
#                         search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author] AND medline[sb] AND "open access"[filter]' ## try searching with filters
#                         search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#                         time.sleep(1)
#                         if int(search_results["Count"]) > 1: ## searching with filter did not work
#                             success = False
#                         elif int(search_results["Count"]) == 1: ## searching with filters worked
#                             success = True


#             ## Finding abstracts
#             if abstract == "" or abstract is None:
#                 id_list = pubmed_search_ids(title)
#                 if len(id_list) == 1: # checking to see if only one query came back
#                     abstract = pubmed_get_abstracts(id_list[0])

#             ## Create article dictionary
#             article_dict = {
#                 "pmid": pmid,
#                 "title": title,
#                 "authors": ", ".join(authors_list),
#                 "has_full_text": "no",
#                 "abstract": abstract,
#                 "full_text": "",
#                 "publication_date": date
#             } 

#             dict_list.append(article_dict)
#             articles_found = articles_found + 1
#             print(str(articles_found) + "/" + str(len(pmids)) + " articles found")

#     except Exception as e:
#         print("Error: " + str(e))

# print("Creating large pdf")
# create_large_pdf(dict_list, output_folder, word)







