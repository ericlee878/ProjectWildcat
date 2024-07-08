##################### CODE THAT COUNTS NUMBER OF ENTREX SEARCHES ################################
''' MAIN (full texts) '''

from Bio import Entrez

Entrez.email = "eric7lee87@gmail.com" 


search_success = 0
search_fail = 0
search_total = 0

success_on_first = 0
success_with_author = 0
success_with_authors = 0
success_with_filters = 0
fail_with_no_results = 0
fail_with_too_many_results = 0


# getting the titles from the excel sheet
path = "../data/papers.xlsx"
articles = extract_articles(path)

for article in articles:
    title = article[0]
    authors = article[1]
    authors_list = [line.strip() for line in authors.strip().split('\n\n')]
    author_query = authors_list[0]

    search_total = search_total + 1
    print(title)

    try:
        search_query = str(title) + '[Title]'# AND ' + str(author_query) + '[Author]' # + ' AND medline[sb] AND "open access"[filter]'
        search_results = Entrez.read(Entrez.esearch(db="pmc", term=search_query, retmax=10, usehistory="y"))
        time.sleep(1)

        if int(search_results["Count"]) == 1: ## if searching only for title worked
            search_success = search_success + 1
            success_on_first = success_on_first + 1
            print("Success (" + str(search_success) + "/" + str(search_total) + ")")

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
                        fail_with_too_many_results = fail_with_too_many_results + 1
                        search_fail = search_fail + 1
                        print("Search fail (" + str(search_fail) + "/" + str(search_total) + ")")
                    elif int(search_results["Count"]) == 1: ## searching with filters worked
                        search_success = search_success + 1
                        success_with_filters = success_with_filters + 1
                        print("Success (" + str(search_success) + "/" + str(search_total) + ")")
                else: ## searching with authors did work
                    search_success = search_success + 1
                    success_with_authors = success_with_authors + 1
                    print("Success (" + str(search_success) + "/" + str(search_total) + ")")
            else: ## searching with author worked
                search_success = search_success + 1
                success_with_author = success_with_author + 1
                print("Success (" + str(search_success) + "/" + str(search_total) + ")")

        elif int(search_results["Count"]) == 0: ## no search results
            search_fail = search_fail + 1
            fail_with_no_results = fail_with_no_results + 1
            print("Search fail (" + str(search_fail) + "/" + str(search_total) + ")")

    except Exception as e:
        search_fail = search_fail + 1
        print("Search fail (" + str(search_fail) + "/" + str(search_total) + ")")
        print('''f"Failed to process search for title: {title}.''' "Error: {str(e)}")
        continue

    print("Articles searched for: " + str(search_total))
    print("Success: " + str(search_success) + "/" + str(search_total))
    print("Failure: " + str(search_fail) + "/" + str(search_total))

    print("Success with only title: " + str(success_on_first))
    print("Success with title + author: " + str(success_with_author))
    print("Success with title + authors: " + str(success_with_authors))
    print("Success with title + authors + filter: " + str(success_with_filters))
    print("Fail with no results: " + str(fail_with_no_results))
    print("Fail with too many results: " + str(fail_with_too_many_results))
    print("\n\n")
