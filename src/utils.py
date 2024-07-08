# importing required modules 
from pypdf import PdfReader 
import os


path = "C:\\Users\\leeh36\\OneDrive - Baxter\\Desktop\\MedExplainer\\data\\TISSEEL_Papers\\"

def pull_text_from_pdfs(path):
    '''
    Parameters:
    ----------------
    path: path of the directory where the pdfs are

    Returns:
    ----------------
    pdfs: list of pdfs where each element is a string that contains the entire text of a pdf 
    '''
    pdfs = []
    # iterates through all the files inside the path
    for i in os.listdir(path):
        # if the file is a pdf
        if i.endswith(".pdf"):
            # creating a pdf reader object 
            reader = PdfReader(path  + i) 

            curr_pdf = ""
            # for each page in the pdf
            for page in reader.pages:
                # extracting text from page 
                text = page.extract_text()
                curr_pdf = curr_pdf + text

            pdfs.append(curr_pdf)

    return pdfs