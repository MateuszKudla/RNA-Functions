from os.path import exists
from os import remove
from os import makedirs
import sys
import time


def remove_if_exist(name):
    if exists(name + '.fasta'):
        remove(name + '.fasta')


def create_folders(name):
    if not exists(name):
        makedirs(name)
        makedirs(name + '\\fasta')
        makedirs(name + '\\xml')


def save_file(folder, seq, i):
    save = open(folder + '\\fasta\\' + str(i) + 'seq.fasta', 'w')
    new_seq = ''

    save.write('>sequence|' + str(i) + '\n')

    for letter in seq:
        if len(new_seq) == 70:
            save.write(new_seq + '\n')
            new_seq = letter
        else:
            new_seq += letter

    save.write(new_seq + '\n\n')


def open_file(file, folder):
    default = open(folder + '\\SeqNumber.txt', 'w')
    i = 1
    f = open(file + '.txt', 'r')

    for line in f:
        seq = line.rstrip()
        save_file(folder, seq, i)
        default.write(str(i) + '. ' + line)
        i += 1

    print("Done")

	
def execute(file_name, folder_name):
    create_folders(folder_name)
    remove_if_exist(file_name)
    open_file(file_name, folder_name)

	
if __name__ == '__main__':
	print("\nFinding functions of RNA fragmnts step 1 - Create fasta files")
	
	folder_name = file_name = input("Enter the name of file contains sequences in txt format: ")

	execute(file_name, folder_name)