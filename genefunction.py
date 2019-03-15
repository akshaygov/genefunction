import sys
import argparse
import requests
import nltk
import re
from bs4 import BeautifulSoup
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize, sent_tokenize
from string import punctuation
from nltk.probability import FreqDist


def removeparens(text):
	#removes text between parentheses

	new_txt = re.sub(r'\(.*?\)', '', text)
	return new_txt

def removeone(lst):
	#when tokenized text has one <
	return lst[0:lst.index("<")]

def removetwo(text):
	#when text has at least one pair of < and >
	new_txt = re.sub('<.*?>', '', text)
	return new_txt

def hasDigit(st):
	#if a word has a digit, then it should not be made lowercase (e.g. 'E3')
	for char in st:
		if (char.isdigit()):
			return True
	return False

def getID(id_s, gene):

	if (len(id_s) == 0):
		return "uncategorized"

	else:
		ID = ""
		for line in id_s:
			line = line.strip().split('\t')
			#Make sure gene names match
			if gene in line[len(line) - 1]:
				ID = line[0]
				break
		if ID == "":
			return "uncategorized"
		return ID

def getFunctionText(content):

	function_tag = content.find(id="function")
	function_text = function_tag.find("span", property="text")
	text = function_text.get_text()
	return text


def joinText(word_tokens):
	txt = ""

	for s in range(len(word_tokens) - 1):
		
		if (word_tokens[s+1] in list(punctuation)):
			txt += word_tokens[s]

		else:
			txt += word_tokens[s] + " "
	
	txt += word_tokens[len(word_tokens) - 1]
	
	return txt


def summaryIndex(text):
	rank = {}
	sent_tokens = sent_tokenize(text.lower())
	word_tokens = word_tokenize(text.lower())
	
	for j in range(len(sent_tokens)):
		rank[j] = 0
	
	word_freq = FreqDist(word_tokens)
	
	for i, sent in enumerate(sent_tokens):
		
		for word in word_tokenize(sent):

			if word in word_freq:
				rank[i] += word_freq[word]

	rank_list = list(rank.values())
	ind = rank_list.index(max(rank_list))

	return ind



def cleanList(out_list, gene):

	for i in range(1, len(out_list)):
		if (nltk.tag.pos_tag([out_list[i]])[0][1] == "PRP") and (nltk.tag.pos_tag([out_list[i-1]])[0][1] == "IN"):
			out_list[i] = "and"
			out_list[i-1] = " "

	if "." in out_list[len(out_list) - 1]:
		end = out_list[len(out_list) - 1]
		repl = end[0:len(end) - 1]
		out_list[len(out_list) - 1] = repl
		
	if gene.lower() in out_list or gene in out_list:
		out_list = out_list[1:]

	return out_list

def cleanTokens(word_tokens):

	if (word_tokens[0] == "and"):
		word_tokens.remove(word_tokens[0])

	if (word_tokens[len(word_tokens) - 1] == "."):
		word_tokens.remove(word_tokens[len(word_tokens) - 1])
		
	if (nltk.pos_tag(word_tokens)[0][1] == "VBZ" and nltk.pos_tag(word_tokens)[1][1] == "VBZ"):
		word_tokens.remove(word_tokens[0])

	word_tokens = word_tokens

def summarizeGene(gene):

	'''
	
	Obtains summarized gene function from UniProt database

	'''


	try:
		source_page = requests.get("https://www.uniprot.org/uniprot/?query={}+AND+organism:9606&sort=score&columns=id,entry name,genes&format=tab".format(gene))
		id_page = BeautifulSoup(source_page.content, 'html.parser')
		id_list = str(id_page).splitlines()


	except (requests.exceptions.ConnectionError, requests.exceptions.RequestException) as e:
		print("UniProt could not be accessed at this time. Please try again later.")
		sys.exit(1)

	
	#no results in UniProt for input

	if getID(id_list, gene) == "uncategorized":
		return "uncategorized"

	ID = getID(id_list, gene)


	content_page = requests.get("https://www.uniprot.org/uniprot/{}".format(ID))
	html_content = BeautifulSoup(content_page.content, 'html.parser')

	full_name = html_content.find("h1", property="name").get_text()
	
	try:

		text = getFunctionText(html_content)
		
		if "<" in text:
			text = text[0:text.index("<")]

		text = text[0:text.rfind('.') + 1]

		#cleaning up parentheses, word tokenizing, sentence tokenizing, removing stopwords
		no_parens_text = removeparens(text)

		word_tok = word_tokenize(no_parens_text)

		sent_tokens = sent_tokenize(no_parens_text)

		good_punct = [p for p in list(punctuation) if p != "."]
		stop_words = set(stopwords.words('english') + good_punct)

		
		word_tokens = [word for word in word_tok if word not in stop_words]
		new_text = ' '.join(word_tokens)
	
		if len(sent_tokens) > 1:
			
			out_ind = summaryIndex(new_text)

			out_text = sent_tokenize(no_parens_text)[out_ind]
			
			
			out_list = out_text.split()
			
			if (not hasDigit(out_list[0])):
				out_list[0] = out_list[0].lower()


			out_list = cleanList(out_list, gene)
			
			char_list = ' '.join(out_list)

			if "<" in char_list and ">" in char_list:
				char_list = removetwo(char_list)

				
			if "<" in char_list and ">" not in char_list:
				char_list = removeone(char_list)

				
				
			new_out_tokens = word_tokenize(char_list)
			
			cleanTokens(new_out_tokens)
			#if no valuable function information in UniProt, returns full name in database
			if "search" in new_out_tokens or "uniprotkb" in new_out_tokens or "chebi" in new_out_tokens:
				return full_name

			return joinText(new_out_tokens)


		else:

			if "search" in word_tok or "uniprotkb" in word_tok or "chebi" in word_tok:
				return full_name

			if "<" in word_tok and ">" in word_tok:
				str_word_tok = removetwo(' '.join(word_tok))
				word_tok = word_tokenize(str_word_tok)

			
			if "<" in word_tok and ">" not in word_tok:
				word_tok = word_tokenize(removeone(' '.join(word_tok)))
					
				
			if "." in word_tok[len(word_tok) - 1]:
				end = word_tok[len(word_tok) - 1]
				word_tok[len(word_tok) - 1] = end[0:len(end)-1]

			if (not hasDigit(word_tok[0])):
				word_tok[0] = word_tok[0].lower()
			
			return joinText(word_tok)
				

	except (ValueError, IndexError, AttributeError):
			

		if (gene in word_tokenize(full_name)):
			return "uncategorized"
		else:
			full_list = list(full_name)
			if (str(full_list[1]).lower() == str(full_list[1])):
				full_list[0] = full_list[0].lower()
			return ''.join(full_list)


def main():
	

	parser = argparse.ArgumentParser()
	parser.add_argument('-gene_name', type=str, help='A single gene whose function will be printed in the terminal (no input/output file required).')
	parser.add_argument('-input_file', type=str, help='A plain text file of newline-delimited gene names.')
	parser.add_argument('-output_file', type=str, help='The output plain text file of newline-delimited gene names with their functions.')
	args = parser.parse_args()
	
	#if no gene_name and both output_file and input_file arguments aren't given
	if (not (args.output_file is not None and args.input_file is not None) and args.gene_name is None):
		print("Both input and output file required, or single gene name.")
	#if no gene_name and both output_file and input_file arguments are given
	elif ((args.output_file is not None and args.input_file is not None) and (args.gene_name is None)):
		
		with open(args.input_file, 'r') as IN, open(args.output_file, 'w') as OUT:
			genes = list(set([line.rstrip('\n') for line in IN]))

			for g in genes:
				
				OUT.write(g.rstrip('\n') + ": " + summarizeGene(g.rstrip('\n')) + '\n')

	#if gene_name and both output_file and input_file arguments given
	elif ((args.gene_name is not None) and (not (args.output_file is None and args.input_file is None))):
		print("Single gene's name cannot be given as argument alongside input/output file arguments")

	else:
		print(str(args.gene_name) + ": " + summarizeGene(args.gene_name))


if __name__ == "__main__":
	main()
	

