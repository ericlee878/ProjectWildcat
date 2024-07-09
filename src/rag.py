import os
import textwrap
import time
from pathlib import Path

from IPython.display import Markdown
from langchain.chains import RetrievalQA
from langchain.prompts import PromptTemplate
from langchain.retrievers import ContextualCompressionRetriever
from langchain.retrievers.document_compressors import FlashrankRerank
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain.vectorstores import Qdrant
from langchain_community.document_loaders import UnstructuredMarkdownLoader
from langchain_community.embeddings.fastembed import FastEmbedEmbeddings
from langchain_core.prompts import ChatPromptTemplate
from langchain_groq import ChatGroq
from llama_parse import LlamaParse


os.environ["GROQ_API_KEY"] = "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"


def print_response(response):
    response_txt = response["result"]
    for chunk in response_txt.split("\n"):
        if not chunk:
            print()
            continue
        print("\n".join(textwrap.wrap(chunk, 100, break_long_words=False)))

instruction = """The provided document is a paper about TISSEEL.
It contains many tables.
Try to be precise while answering the questions"""

parser = LlamaParse(
    api_key=os.environ["GROQ_API_KEY"],
    result_type="markdown",
    parsing_instruction=instruction,
    max_timeout=5000,
)

llama_parse_documents = await parser.aload_data("../data/TISSEEL_Papers/Angioli-2912.pdf")

time.sleep(5)

parsed_doc = llama_parse_documents[0]

print(Markdown(parsed_doc.text[:4096]))