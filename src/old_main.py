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

# from Bio import Entrez

# Entrez.email = "eric7lee87@gmail.com" 


# search_success = 0
# search_fail = 0
# search_total = 0
# # getting the titles from the excel sheet
# path = "../data/papers.xlsx"
# articles = extract_articles(path)

# for sheet in articles.keys():
#     pairs = articles[sheet]
#     for pair in pairs:
#         title = pair[0]
#         authors = pair[1]
#         authors_list = [line.strip() for line in authors.strip().split('\n\n')]
#         author_query = authors_list[0]

#         if title == 'TITLE': ## skip the column that says title
#             continue

#         search_total = search_total + 1
#         print(title)

#         try:
#             search_query = str(title) + '[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
#             search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#             time.sleep(1)
#             if int(search_results["Count"]) == 1:
#                 search_success = search_success + 1
#                 print("Success (" + str(search_success) + "/" + str(search_total) + ")")
#                 print(title)
#                 break
#                 # batch_size = 100
#                 # for start in range(0, int(search_results["Count"]), batch_size):
#                 #     end = min(int(search_results["Count"]), start + batch_size)
#                 #     handle = Entrez.efetch(db="pubmed", rettype="full", retmode="xml", 
#                 #                         retstart=start, retmax=batch_size,
#                 #                         webenv=search_results["WebEnv"], 
#                 #                         query_key=search_results["QueryKey"])
#                 #     data = handle.read()
#                 #     # Process or save the data here
#                 #     fulltexts = open("../data/fulltexts.txt", "a")
#                 #     fulltexts.write(str(data) + ",\n")
#                 #     fulltexts.close()
#                 #     handle.close()
#                 #     print("here")
#             elif int(search_results["Count"]) > 1:
#                 search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]'
#                 search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#                 time.sleep(1)
#                 if int(search_results["Count"]) > 1:
#                     for i in range(1, len(authors_list)):
#                         author_query = " OR "+ author_query
#                     search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author]'
#                     search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#                     time.sleep(1)
#                     if int(search_results["Count"]) > 1:
#                         search_query = str(title) + '[Title] AND ' + str(author_query) + '[Author] AND medline[sb] AND "open access"[filter]'
#                         search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#                         time.sleep(1)
#                         if int(search_results["Count"]) > 1:
#                             print("Search resulted in too many results.")
#                         else:
#                             search_success = search_success + 1
#                             print("Success (" + str(search_success) + "/" + str(search_total) + ")")
#                     else:
#                         search_success = search_success + 1
#                         print("Success (" + str(search_success) + "/" + str(search_total) + ")")
#                 else:
#                     search_success = search_success + 1
#                     print("Success (" + str(search_success) + "/" + str(search_total) + ")")
#             elif int(search_results["Count"]) == 0:
#                 search_fail = search_fail + 1
#                 print("Search fail (" + str(search_fail) + "/" + str(search_total) + ")")
#                 # search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#                 # time.sleep(1)
#                 # if int(search_results["Count"]) == 0:
#                 #     search_results = Entrez.read(Entrez.esearch(db="nlmcatalog", term=search_query, retmax=10, usehistory="y"))
#                 #     time.sleep(1)
#                 #     if int(search_results["Count"]) == 0:
#                 #         search_fail = search_fail + 1
#                 #         print("Search fail (" + str(search_fail) + "/" + str(search_total) + ")")
#                 #     else:
#                 #         search_success = search_success + 1
#                 #         print("Success (" + str(search_success) + "/" + str(search_total) + ")")
#                 # else:
#                 #     search_success = search_success + 1
#                 #     print("Success (" + str(search_success) + "/" + str(search_total) + ")")
#         except Exception as e:
#             search_fail = search_fail + 1
#             print("Search fail (" + str(search_fail) + "/" + str(search_total) + ")")
#             print(f"Failed to process search for title: {title}. Error: {str(e)}")


# print("Search success #" + str(search_success))
# print("Search failure #" + search_results["Count"] + " search #" + str(search_fail))

################################## UNUSED FUNCTIONS ###############################
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


###########################################################################################################################################

# search_query = 'Fibrin glue patch for pacemaker lead perforation of the right ventricular free wall: A case report.[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
# search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
# batch_size = 100
# for start in range(0, int(search_results["Count"]), batch_size):
#     end = min(int(search_results["Count"]), start + batch_size)
#     handle = Entrez.efetch(db="pmc", rettype="full", retmode="xml", 
#                            retstart=start, retmax=batch_size,
#                            webenv=search_results["WebEnv"], 
#                            query_key=search_results["QueryKey"])
#     data = handle.read()
#     # Process or save the data here
#     fulltexts = open("../data/fulltexts.txt", "a")
#     fulltexts.write(str(data) + ",\n")
#     fulltexts.close()
#     handle.close()
#     print("here")



# url = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/17299597/unicode"
# response = requests.get(url, verify='4a18af8a15c4bf6286ff6f0658958009af7cee741ce4ba4ab41a9bc4dfaf6bdb')


# headers = {
#     'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.46 (KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36',
#     'From': 'eric7lee87@gmail.com',
#     'HTTP': 'Hydra/1.3.15'
# }
# headers = {
#     "Accept": "*/*",
#     "Accept-Encoding": "gzip, deflate, br, zstd",
#     "Accept-Language": "en-US,en;q=0.9",
#     "Connection": "keep-alive",
#     "Content-Length": "0",
#     "Cookie": "ar_debug=1",
#     "Host": "www.google-analytics.com",
#     "Origin": "https://pubmed.ncbi.nlm.nih.gov",
#     "Referer": "https://pubmed.ncbi.nlm.nih.gov/",
#     "Sec-Ch-Ua": '"Not/A)Brand";v="8", "Chromium";v="126", "Google Chrome";v="126"',
#     "Sec-Ch-Ua-Mobile": "?0",
#     "Sec-Ch-Ua-Platform": '"Windows"',
#     "Sec-Fetch-Dest": "empty",
#     "Sec-Fetch-Mode": "no-cors",
#     "Sec-Fetch-Site": "cross-site",
#     "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/538.36 (KHTML, like Gecko) Chrome/126.0.0.0 Safari/537.36"
# }
# r = requests.get('https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/17299597/unicode', headers=headers)

# from selenium import webdriver
# from selenium.webdriver.chrome.options import Options

# chrome_options = Options()
# #Set up user agent to avoid bot detection.
# user_agent = 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/60.0.3112.50 Safari/537.36'
# #Specify headless mode
# chrome_options.add_argument("--headless")
# #Add user agent
# chrome_options.add_argument(f'user-agent={user_agent}')
# DRIVER_PATH = "C://Program Files (x86)//Google//Chrome//Application//chrome.exe" 
# driver = webdriver.Chrome(DRIVER_PATH, options=chrome_options)
# driver.get("https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/17299597/unicode")
# page_source = driver.page_source
# print(page_source)










# full_text = ""
# search_query = 'Fibrin glue patch for pacemaker lead perforation of the right ventricular free wall: A case report[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
# search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
# time.sleep(1)
# batch_size = 100
# for start in range(0, int(search_results["Count"]), batch_size):
#     end = min(int(search_results["Count"]), start + batch_size)
#     handle = Entrez.efetch(db="pmc", rettype="full", retmode="xml", 
#                         retstart=start, retmax=batch_size,
#                         webenv=search_results["WebEnv"], 
#                         query_key=search_results["QueryKey"])
#     data = handle.read()
#     soup = BeautifulSoup(data, features="xml")
#     body_soup = soup.find("body")
#     p_soup = body_soup.find_all("p")
#     for p in p_soup:
#         full_text = full_text + p.text
#     # fulltexts = open("../data/fulltext_test.txt", "w")
#     # fulltexts.write(str(data))# + ",\n")
#     # fulltexts.close()
#     print(full_text)
#     handle.close()














##################### CODE THAT COUNTS NUMBER OF ENTREX SEARCHES ################################
''' MAIN (full texts) '''

# from Bio import Entrez

# Entrez.email = "eric7lee87@gmail.com" 


# search_success = 0
# search_fail = 0
# search_total = 0

# success_on_first = 0
# success_with_author = 0
# success_with_authors = 0
# success_with_filters = 0
# fail_with_no_results = 0
# fail_with_too_many_results = 0


# # getting the titles from the excel sheet
# path = "../data/papers.xlsx"
# articles = extract_articles(path)

# for article in articles:
#     title = article[0]
#     authors = article[1]
#     authors_list = [line.strip() for line in authors.strip().split('\n\n')]
#     author_query = authors_list[0]

#     search_total = search_total + 1
#     print(title)

#     try:
#         search_query = str(title) + '[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
#         search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
#         time.sleep(1)

#         if int(search_results["Count"]) == 1: ## if searching only for title worked
#             search_success = search_success + 1
#             success_on_first = success_on_first + 1
#             print("Success (" + str(search_success) + "/" + str(search_total) + ")")

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
#                         fail_with_too_many_results = fail_with_too_many_results + 1
#                         search_fail = search_fail + 1
#                         print("Search fail (" + str(search_fail) + "/" + str(search_total) + ")")
#                     elif int(search_results["Count"]) == 1: ## searching with filters worked
#                         search_success = search_success + 1
#                         success_with_filters = success_with_filters + 1
#                         print("Success (" + str(search_success) + "/" + str(search_total) + ")")
#                 else: ## searching with authors did work
#                     search_success = search_success + 1
#                     success_with_authors = success_with_authors + 1
#                     print("Success (" + str(search_success) + "/" + str(search_total) + ")")
#             else: ## searching with author worked
#                 search_success = search_success + 1
#                 success_with_author = success_with_author + 1
#                 print("Success (" + str(search_success) + "/" + str(search_total) + ")")

#         elif int(search_results["Count"]) == 0: ## no search results
#             search_fail = search_fail + 1
#             fail_with_no_results = fail_with_no_results + 1
#             print("Search fail (" + str(search_fail) + "/" + str(search_total) + ")")

#     except Exception as e:
#         search_fail = search_fail + 1
#         print("Search fail (" + str(search_fail) + "/" + str(search_total) + ")")
#         print('''f"Failed to process search for title: {title}.''' "Error: {str(e)}")
#         continue

#     print("Articles searched for: " + str(search_total))
#     print("Success: " + str(search_success) + "/" + str(search_total))
#     print("Failure: " + str(search_fail) + "/" + str(search_total))

#     print("Success with only title: " + str(success_on_first))
#     print("Success with title + author: " + str(success_with_author))
#     print("Success with title + authors: " + str(success_with_authors))
#     print("Success with title + authors + filter: " + str(success_with_filters))
#     print("Fail with no results: " + str(fail_with_no_results))
#     print("Fail with too many results: " + str(fail_with_too_many_results))
#     print("\n\n")





##################### CODE THAT COUNTS NUMBER OF ARTICLES EXTRACTED ################################


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
    articles_count = {}

    # Define variable to load the dataframe
    dataframe = openpyxl.load_workbook(path)

    # iterates through each sheet and adds each title into list
    for i in range(len(dataframe.sheetnames)):
        sheetname = dataframe.sheetnames[i]
        curr_sheet = dataframe._sheets[i]
        count = 0
        for i in range(1, curr_sheet.max_row):
            title = curr_sheet.cell(row = i, column = 8).value
            authors = curr_sheet.cell(row = i, column = 5).value
            if title is not None and authors is not None and title != 'TITLE' and authors != 'AUTHORS':
                pair = [title, authors]
                articles.append(pair)
                count += 1
        articles_count[sheetname] = count

    return articles_count


# def create_large_pdf(dict_list, output_folder, filename):
#     filepath = os.path.join(output_folder, filename)
    
#     pdf = FPDF()
#     pdf.add_page()
    
#     # Set font for Filename
#     pdf.set_font("Arial", 'B', size=30)

#     # Add Filename
#     pdf.multi_cell(0, 10, txt=filename.capitalize(), align='C')
#     pdf.ln(10)

#     pdf_num = 0

#     for d in dict_list:
#         pdf.add_page()
#         # Set font for title
#         pdf.set_font("Arial", 'B', size=16)
        
#         # Add title
#         pdf.multi_cell(0, 10, txt=d['title'], align='C')
#         pdf.ln(10)
        
#         # Set font for normal text
#         pdf.set_font("Arial", size=12)
        
#         # Add authors
#         pdf.multi_cell(0, 10, txt=f"Authors: {d['authors']}", align='L')
#         pdf.ln(5)
        
#         # Add publication date
#         pdf.multi_cell(0, 10, txt=f"Publication Date: {d['publication_date']}", align='L')
#         pdf.ln(10)
        
#         # Add abstract
#         if d['abstract'] is not None:
#             pdf.set_font("Arial", 'B', size=14)
#             pdf.cell(0, 10, txt="Abstract:", ln=1)
#             pdf.set_font("Arial", size=12)
#             pdf.multi_cell(0, 10, txt=d['abstract'])
#             pdf.ln(10)
#         else:
#             pdf.set_font("Arial", 'B', size=14)
#             pdf.cell(0, 10, txt="Abstract:", ln=1)
#             pdf.set_font("Arial", size=12)
#             pdf.multi_cell(0, 10, txt='Not Available')
#             pdf.ln(10)

#         if pdf.page_no() >= 198:
#             print("New PDF page")
#             pdf.output(filepath + "_" + str(pdf_num) + ".pdf")
#             pdf = FPDF()
#             pdf.add_page()
#             pdf_num = pdf_num + 1
        
#         # # Add full text if available
#         # if d['has_full_text'] == "yes" and d['has_full_text'] is not None:
#         #     pdf.add_page()
#         #     pdf.set_font("Arial", 'B', size=14)
#         #     pdf.cell(0, 10, txt="Full Text:", ln=1)
#         #     pdf.set_font("Arial", size=12)
#         #     pdf.multi_cell(0, 10, txt=d['full_text'])
#         # else:
#         #     pdf.add_page()
#         #     pdf.set_font("Arial", 'B', size=14)
#         #     pdf.cell(0, 10, txt="Full Text:", ln=1)
#         #     pdf.set_font("Arial", size=12)
#         #     pdf.multi_cell(0, 10, txt='Not Available')
#     pdf.output(filepath + "_" + str(pdf_num) + ".pdf")



