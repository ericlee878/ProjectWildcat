import os
import requests
from bs4 import BeautifulSoup
import openpyxl
import time
import json
from fpdf import FPDF
from Bio import Entrez
import re
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT

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
    articles = {}

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
                articles[title] = [title, authors, pmid, date, abstract, mesh, keywords]

    return articles

def pubmed_search_articles(search):
    '''
    Parameters:
    ----------------
    search: the search query on pubmed

    Returns:
    ----------------
    id_list: the list of the urls
    '''
    id_list = []

    # the url of the article
    url = "https://pubmed.ncbi.nlm.nih.gov/?term=" + search + "&filter=simsearch3.fft"

    # the response from the url
    response = requests.get(url)

    # sleep for 0.4 seconds to not get blocked
    time.sleep(0.4)

    # the beautiful soup of the url
    soup = BeautifulSoup(response.text, features="lxml")

    try:
        # Find all 'a' tags with class 'docsum-title'
        a_tags = soup.find_all('a', class_='docsum-title')

        # Extract and print the data-article-id for each tag
        for i, tag in enumerate(a_tags, 1):
            data_article_id = tag.get('data-article-id')
            id_list.append(data_article_id)

        return id_list
    except:
        return None



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

def pubmed_extract_metadata(id):
    '''
    Parameters:
    ----------------
    id: the id of the article

    Returns:
    ----------------
    metadata: the metadata of the article
    '''
    metadata = []

    url = "https://pubmed.ncbi.nlm.nih.gov/" + str(id)

    # the response from the url
    response = requests.get(url)

    # sleep for 0.4 seconds to not get blocked
    time.sleep(0.4)

    # the beautiful soup of the url
    soup = BeautifulSoup(response.text, features="lxml")

    title = ""
    authors = []
    date = ""
    # Find the h1 tag with class "heading-title"
    h1_tag = soup.find('h1', class_='heading-title')

    if h1_tag:
        # Extract the text, strip whitespace
        title = h1_tag.get_text(strip=True)

    # Find all 'a' tags with class "full-name"
    author_tags = soup.find_all('a', class_='full-name')

    # Extract and print the author names
    for tag in author_tags:
        author_name = tag.get_text(strip=True)
        authors.append(author_name)

    # Find the span tag with class "cit"
    cit_span = soup.find('span', class_='cit')

    if cit_span:
        # Extract the text content
        date = cit_span.get_text(strip=True)
        match = re.match(r'(\d{4} \w{3} \d{2})', date)
        if match:
            date = match.group(1)

    metadata = [title, authors, date, id]
    return metadata



def get_full_text_from_search_results(search_results_xml):
    # batch_size = 100
    # for start in range(0, int(search_results["Count"]), batch_size):
    #     end = min(int(search_results["Count"]), start + batch_size)
    #     handle = Entrez.efetch(db="pmc", rettype="full", retmode="xml", 
    #                         retstart=start, retmax=batch_size,
    #                         webenv=search_results["WebEnv"], 
    #                         query_key=search_results["QueryKey"])
    #     data = handle.read()
        soup = BeautifulSoup(search_results_xml, features="xml")
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

    pdf.output(filepath)
    print(f"Created PDF for article {article_dict['pmid']}")


def create_large_pdf(dict_list, output_folder, filename, include_fulltexts):
    if include_fulltexts:
        filepath = os.path.join(output_folder, filename)
        
        doc = SimpleDocTemplate(filepath + ".pdf", pagesize=letter)
        styles = getSampleStyleSheet()
        story = []

        # Add filename
        story.append(Paragraph(filename.capitalize(), styles['Title']))
        story.append(Spacer(1, 12))

        for d in dict_list:
            # Add title
            story.append(Paragraph(d['title'], styles['Heading1']))
            story.append(Spacer(1, 12))

            # Add authors
            story.append(Paragraph(f"Authors: {d['authors']}", styles['Normal']))
            story.append(Spacer(1, 6))

            # Add publication date
            story.append(Paragraph(f"Publication Date: {d['publication_date']}", styles['Normal']))
            story.append(Spacer(1, 12))

            # Add PMID
            story.append(Paragraph(f"PMID: {d['pmid']}", styles['Normal']))
            story.append(Spacer(1, 12))

            # Add abstract
            story.append(Paragraph("Abstract:", styles['Heading2']))
            abstract_text = d['abstract'] if d['abstract'] is not None else 'Not Available'
            story.append(Paragraph(abstract_text, styles['Normal']))
            story.append(Spacer(1, 12))

            # Add fulltext
            story.append(Paragraph("Full Text:", styles['Heading2']))
            abstract_text = d['full_text'] if d['full_text'] is not None else 'Not Available'
            print(d['full_text'])
            story.append(Paragraph(abstract_text, styles['Normal']))
            story.append(Spacer(1, 12))

        doc.build(story)
    else:
        filepath = os.path.join(output_folder, filename)
        
        doc = SimpleDocTemplate(filepath + ".pdf", pagesize=letter)
        styles = getSampleStyleSheet()
        story = []

        # Add filename
        story.append(Paragraph(filename.capitalize(), styles['Title']))
        story.append(Spacer(1, 12))

        for d in dict_list:
            # Add title
            story.append(Paragraph(d['title'], styles['Heading1']))
            story.append(Spacer(1, 12))

            # Add authors
            story.append(Paragraph(f"Authors: {d['authors']}", styles['Normal']))
            story.append(Spacer(1, 6))

            # Add publication date
            story.append(Paragraph(f"Publication Date: {d['publication_date']}", styles['Normal']))
            story.append(Spacer(1, 12))

            # Add abstract
            story.append(Paragraph("Abstract:", styles['Heading2']))
            abstract_text = d['abstract'] if d['abstract'] is not None else 'Not Available'
            story.append(Paragraph(abstract_text, styles['Normal']))
            story.append(Spacer(1, 12))

        doc.build(story)

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
#     title = article[0]
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
#                 if title not in keyword_dict[word]:
#                     keyword_dict[word].append(title)
#             else:
#                 keyword_dict[word] = [title]

# json_file_path = "../data/title_dictionary.json"

# if not os.path.exists(json_file_path):
#     with open(json_file_path, 'w') as f:
#         json.dump([], f)

# with open(json_file_path, 'w') as f:
#     json.dump(keyword_dict, f, indent=2)

############### CODE THAT DOWNLOADS PDFS FROM EXCEL SHEET ################
def download_pdf_excel(word, include_fulltexts):
    # getting the titles from the excel sheet
    path = "../data/papers.xlsx"
    articles = extract_metadata(path)

    output_folder = "..\\pdfs\\" + word
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)      
        
    dict_list = []

    articles_found = 0

    ## getting the relevant PMIDs from json

    file = open('../data/title_dictionary.json')
    data = json.load(file)
    titles = data[word.lower()]
    
    for title in titles:
        if title in articles.keys():
            metadata = articles[title]
            authors = metadata[1]
            pmid = metadata[2]
            date = metadata[3]
            abstract = metadata[4]
            authors_list = [line.strip() for line in authors.strip().split('\n\n')]
            author_query = authors_list[0]
            mesh = metadata[5]
            keywords = metadata[6]
            fulltext = ""

        if title == 'TITLE':
            continue

        search_query = str(title) + '[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
        search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
        time.sleep(1)
        success = False

        if int(search_results["Count"]) == 1: ## if searching only for title worked
            success = True

        elif int(search_results["Count"]) > 1: ## if too many queries popped up
            search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with author
            search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
            time.sleep(1)
            if int(search_results["Count"]) > 1: ## searching with author did not work
                for i in range(1, len(authors_list)):
                    author_query = " OR "+ author_query
                search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with authors
                search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
                time.sleep(1)
                if int(search_results["Count"]) > 1: ## searching with authors did not work
                    search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author] AND medline[sb] AND "open access"[filter]' ## try searching with filters
                    search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
                    time.sleep(1)
                    if int(search_results["Count"]) > 1: ## searching with filter did not work
                        success = False
                    elif int(search_results["Count"]) == 1: ## searching with filters worked
                        success = True

        if success:
            fulltext = get_full_text_from_search_results(search_results)
            if fulltext is None:
                fulltext = ""
            print("came here")

        ## Finding abstracts
        if abstract == "" or abstract is None:
            id_list = pubmed_search_ids(title)
            if len(id_list) == 1: # checking to see if only one query came back
                abstract = pubmed_get_abstracts(id_list[0])

        ## Create article dictionary
        article_dict = {
            "pmid": pmid,
            "title": title,
            "authors": ", ".join(authors_list),
            "abstract": abstract,
            "full_text": fulltext,
            "publication_date": date
        }
        dict_list.append(article_dict)
        articles_found = articles_found + 1
        # print(str(articles_found) + "/" + str(len(id_list)) + " articles found")
        print("Article " + str(pmid) + " found.")

    print("Creating large pdf...")
    create_large_pdf(dict_list, output_folder, word, include_fulltexts)
    print("PDF created.")



    for article in articles:
        try:
            title = article[0]
            authors = article[1]
            pmid = article[2]
            date = article[3]
            abstract = article[4]
            authors_list = [line.strip() for line in authors.strip().split('\n\n')]
            author_query = authors_list[0]
            mesh = article[5]
            keywords = article[6]
            fulltext = ""

            if title in titles:

                if title == 'TITLE':
                    continue

                search_query = str(title) + '[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
                search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
                time.sleep(1)
                success = False

                if int(search_results["Count"]) == 1: ## if searching only for title worked
                    success = True

                elif int(search_results["Count"]) > 1: ## if too many queries popped up
                    search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with author
                    search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
                    time.sleep(1)
                    if int(search_results["Count"]) > 1: ## searching with author did not work
                        for i in range(1, len(authors_list)):
                            author_query = " OR "+ author_query
                        search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with authors
                        search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
                        time.sleep(1)
                        if int(search_results["Count"]) > 1: ## searching with authors did not work
                            search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author] AND medline[sb] AND "open access"[filter]' ## try searching with filters
                            search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
                            time.sleep(1)
                            if int(search_results["Count"]) > 1: ## searching with filter did not work
                                success = False
                            elif int(search_results["Count"]) == 1: ## searching with filters worked
                                success = True

                if success:
                    fulltext = get_full_text_from_search_results(search_results)
                    if fulltext is None:
                        fulltext = ""
                    print("came here")

                ## Finding abstracts
                if abstract == "" or abstract is None:
                    id_list = pubmed_search_ids(title)
                    if len(id_list) == 1: # checking to see if only one query came back
                        abstract = pubmed_get_abstracts(id_list[0])

                ## Create article dictionary
                article_dict = {
                    "pmid": pmid,
                    "title": title,
                    "authors": ", ".join(authors_list),
                    "abstract": abstract,
                    "full_text": fulltext,
                    "publication_date": date
                } 

                dict_list.append(article_dict)
                articles_found = articles_found + 1
                # print(str(articles_found) + "/" + str(len(id_list)) + " articles found")
                print("Article " + str(pmid) + " found.")

        except Exception as e:
            print("Error: " + str(e))

    print("Creating large pdf...")
    create_large_pdf(dict_list, output_folder, word, include_fulltexts)
    print("PDF created.")

############### CODE THAT DOWNLOADS PDFS NOT IN EXCEL SHEET ################
def download_pdf_not_excel(word, theme, include_fulltexts):
    # getting the titles from the excel sheet
    path = "../data/papers.xlsx"
    articles = extract_metadata(path)

    output_folder = "..\\pdfs\\" + word
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)      
        
    dict_list = []

    articles_found = 0

    ## getting the relevant PMIDs from json
    file = open('../data/dictionary.json')
    data = json.load(file)

    if word in data.keys():
        pmids = data[word.lower()]
        for article in articles:
            try:
                title = article[0]
                authors = article[1]
                pmid = article[2]
                date = article[3]
                abstract = article[4]
                authors_list = [line.strip() for line in authors.strip().split('\n\n')]
                author_query = authors_list[0]
                mesh = article[5]
                keywords = article[6]
                fulltext = ""

                if pmid in pmids:

                    if title == 'TITLE':
                        continue

                    search_query = str(title) + '[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
                    search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
                    time.sleep(1)
                    success = False

                    if int(search_results["Count"]) == 1: ## if searching only for title worked
                        success = True

                    elif int(search_results["Count"]) > 1: ## if too many queries popped up
                        search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with author
                        search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
                        time.sleep(1)
                        if int(search_results["Count"]) > 1: ## searching with author did not work
                            for i in range(1, len(authors_list)):
                                author_query = " OR "+ author_query
                            search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with authors
                            search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
                            time.sleep(1)
                            if int(search_results["Count"]) > 1: ## searching with authors did not work
                                search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author] AND medline[sb] AND "open access"[filter]' ## try searching with filters
                                search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
                                time.sleep(1)
                                if int(search_results["Count"]) > 1: ## searching with filter did not work
                                    success = False
                                elif int(search_results["Count"]) == 1: ## searching with filters worked
                                    success = True


                    ## Finding abstracts
                    if abstract == "" or abstract is None:
                        id_list = pubmed_search_ids(title)
                        if len(id_list) == 1: # checking to see if only one query came back
                            abstract = pubmed_get_abstracts(id_list[0])

                    if success:
                        fulltext = get_full_text_from_search_results(search_results)
                        if fulltext is None:
                            fulltext = ""
                        print("came here")

                    print(fulltext)

                    ## Create article dictionary
                    article_dict = {
                        "pmid": pmid,
                        "title": title,
                        "authors": ", ".join(authors_list),
                        "abstract": abstract,
                        "full_text": fulltext,
                        "publication_date": date
                    } 

                    dict_list.append(article_dict)
                    articles_found = articles_found + 1
                    print(str(articles_found) + "/" + str(len(pmids)) + " articles found")

            except Exception as e:
                print("Error: " + str(e))

    id_list = pubmed_search_articles(word + " " + theme)
    for id in id_list:
        if len(dict_list) >= 10:
            break
        abstract = pubmed_get_abstracts(id)
        metadata = pubmed_extract_metadata(id)
        title = metadata[0]
        authors = metadata[1]
        date = metadata[2]

        ##get fulltext here
        search_query = str(title) + '[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
        search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
        time.sleep(1)
        success = False

        if int(search_results["Count"]) == 1: ## if searching only for title worked
            success = True

        elif int(search_results["Count"]) > 1: ## if too many queries popped up
            search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with author
            search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
            time.sleep(1)
            if int(search_results["Count"]) > 1: ## searching with author did not work
                for i in range(1, len(authors_list)):
                    author_query = " OR "+ author_query
                search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with authors
                search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
                time.sleep(1)
                if int(search_results["Count"]) > 1: ## searching with authors did not work
                    search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author] AND medline[sb] AND "open access"[filter]' ## try searching with filters
                    search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
                    time.sleep(1)
                    if int(search_results["Count"]) > 1: ## searching with filter did not work
                        success = False
                    elif int(search_results["Count"]) == 1: ## searching with filters worked
                        success = True

            if success:
                fulltext = get_full_text_from_search_results(search_results)
                if fulltext is None:
                    fulltext = ""
                print("came here")

        ## Create article dictionary
        article_dict = {
            "pmid": id,
            "title": title,
            "authors": ", ".join(authors),
            "abstract": abstract,
            "full_text": fulltext,
            "publication_date": date
        } 

        dict_list.append(article_dict)
        articles_found = articles_found + 1
        # print(str(articles_found) + "/" + str(len(id_list)) + " articles found")
        print("Article " + str(id) + " found.")

    print("Creating large pdf...")
    create_large_pdf(dict_list, output_folder, word, include_fulltexts)
    print("PDF created.")

##########################  MAIN CODE THAT RUNS THE DOWNLOADING PDFS  ################################

# file = open('../data/dictionary.json')
# data = json.load(file)
# keyword = input("Enter your keyword:\n")
# theme = "Fibren Glue"

# if keyword not in data.keys() or len(data[keyword]) < 3:
#     ## NOT ENOUGH ARTICLES
#     response = input("Not enough articles available in database. Fill the rest with other articles from PMC?\n")
#     if response.lower() == 'yes' or response.lower() == 'y':
#         download_pdf_not_excel(keyword, theme, True)
#     else:
#         print("Sorry, not enough PDFs.")
# else:
#     download_pdf_excel(keyword, True)
    
#######################################################################################################


# full_text = ""
# search_query = 'Fibrin[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
# search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
# time.sleep(1)
# batch_size = 1
# for start in range(0, 1):#int(search_results["Count"]), batch_size):
#     end = min(int(search_results["Count"]), start + batch_size)
#     handle = Entrez.efetch(db="pmc", rettype="full", retmode="xml", 
#                         retstart=start, retmax=batch_size,
#                         webenv=search_results["WebEnv"], 
#                         query_key=search_results["QueryKey"])
#     data = handle.read()
#     print("came here")
#     soup = BeautifulSoup(data, features="xml")
#     fulltexts = open("../data/fulltext_test.txt", "w")
#     fulltexts.write(get_full_text_from_search_results(search_results))# + ",\n")
#     fulltexts.close()
#     body_soup = soup.find("body")
#     p_soup = body_soup.find_all("p")
#     for p in p_soup:
#         full_text = full_text + p.text
#     print(full_text)
#     handle.close()



# get_full_text_from_search_results(search_results)

output_folder = "..\\pdfs"
word = "test"

title = "Fibrin glue does not assist migration and proliferation of chondrocytes in collagenic membranes: an in vitro study"
authors = ["Filippo Migliorini", "Julia Prinz"]
pmid = "34572476"
date = "2021"
abstract = "In this study, the influence of two subfractions (with previously proven anti-cancer properties) isolated from wood rot fungus Cerrena unicolor on the formation of a fibrin clot was investigated in the context of potential use as fibrin glue and sealant enhancers and potential wound healing agents. With the use of ROTEM thromboelastometry, we demonstrated that, in the presence of fibrinogen and thrombin, the S6 fraction accelerated the formation of a fibrin clot, had a positive effect on its elasticity modulus, and enhanced the degree of fibrin cross-linking. The S5 fraction alone showed no influence on the fibrin coagulation process; however, in the presence of fibrin, it exhibited a decrease in anti-proliferative properties against the HT-29 line, while it increased the proliferation of cells in general at a concentration of 100 µg/mL. Both fractions retained their proapoptotic properties to a lesser degree. In combination with the S6 fraction in the ratio of 1:1 and 1:3, the fractions contributed to increased inhibition of the activity of matrix metalloproteinases (MMPs). This may suggest anti-metastatic activity of the combined fractions. In conclusion, the potential of the fractions isolated from the C. unicolor secretome to be used as a means of improving the wound healing process was presented. The potential for delivering agents with cytostatic properties introduced far from the site of action or exerting a pro-proliferative effect at the wound site with the aid of a fibrin sealant was demonstrated."
authors_list = ["Filippo Migliorini", "Julia Prinz"]
author_query = authors_list[0]
fulltext = ""



search_query = str(title) + '[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
time.sleep(1)
success = False

if int(search_results["Count"]) == 1: ## if searching only for title worked
    success = True

elif int(search_results["Count"]) > 1: ## if too many queries popped up
    print("Too many searches: " + str(search_results["Count"]))
    search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with author
    search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
    time.sleep(1)
    if int(search_results["Count"]) > 1: ## searching with author did not work
        print("Too many searches: " + str(search_results["Count"]))
        for i in range(1, len(authors_list)):
            author_query = " OR "+ author_query
        search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]' ## try searching with authors
        search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
        time.sleep(1)
        if int(search_results["Count"]) > 1: ## searching with authors did not work
            print("Still too many searches: " + str(str(search_results["Count"])))
            search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author] AND medline[sb] AND "open access"[filter]' ## try searching with filters
            search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
            time.sleep(1)
            if int(search_results["Count"]) > 1: ## searching with filter did not work
                success = False
            elif int(search_results["Count"]) == 1: ## searching with filters worked
                success = True
            batch_size = 100
            for start in range(0, int(search_results["Count"]), batch_size):
                end = min(int(search_results["Count"]), start + batch_size)
                handle = Entrez.efetch(db="pmc", rettype="full", retmode="xml", 
                                    retstart=start, retmax=batch_size,
                                    webenv=search_results["WebEnv"], 
                                    query_key=search_results["QueryKey"])
                data = handle.read()
                def extract_xref_p_values(text):
                    # Create a BeautifulSoup object
                    soup = BeautifulSoup(text, 'xml')
                    
                    # Find all xref tags with ref-type="bibr"
                    ps = soup.find_all('p')
                    
                    # Extract the text content (p value) from each xref
                    p_values = [p.text for p in ps]
                    
                    return p_values

                p_values = extract_xref_p_values(data)
                print("Extracted 'p' values:", p_values)

else:
    print("Article not found")

if success:
    fulltext = get_full_text_from_search_results(search_results)
    print("Is Success")

 ## Finding abstracts
if abstract == "" or abstract is None:
    id_list = pubmed_search_ids(title)
    if len(id_list) == 1: # checking to see if only one query came back
        abstract = pubmed_get_abstracts(id_list[0])

## Create article dictionary
article_dict = {
    "pmid": pmid,
    "title": title,
    "authors": ", ".join(authors_list),
    "abstract": abstract,
    "full_text": fulltext,
    "publication_date": date
}


print("Creating large pdf...")
create_large_pdf([article_dict], output_folder, word, True)
print("PDF created.")


