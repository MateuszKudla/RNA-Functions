from Bio import SeqIO
from Bio.Blast import NCBIWWW
from os.path import exists
import sys


def fasta_exist(i, folder):
    if exists(folder + '\\fasta\\' + str(i) + 'seq.fasta'):
        return True
    return False


def xml_exist(i, folder):
    if exists(folder + '\\xml\\' + str(i) + 'seq.xml'):
        return True
    return False


def generate_xml(i, folder):
    if (xml_exist(i, folder)):
        return

    my_query = SeqIO.read(folder + '\\fasta\\' + str(i) + 'seq.fasta', format="fasta")
    result_handle = NCBIWWW.qblast("blastn", "nt", my_query.seq)
    blast_result = open(folder + '\\xml\\' + str(i) + "seq.xml", "w")
    blast_result.write(result_handle.read())
    blast_result.close()
    result_handle.close()

	
def execute(folder,):
    i = 1

    while (fasta_exist(i, folder_name)):
        print('Searching ' + str(i) + '...\t')
        generate_xml(i, folder_name)
        print('Done')
        i += 1


if __name__ == "__main__":
	print("\nFinding functions of RNA fragmnts step 2 - Create xml files")

	folder_name = input("Enter the name of the test folder: ")
	execute(folder_name)
