from Bio import SearchIO
from Bio import Entrez
from os.path import exists
import sys

email = ""
all_seq = []


class Result:
    def __init__(self, name, hit_start, hit_end, region):
        self.name = name
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.region = region

    def show_result(self):
        return self.name + " " + str(self.hit_start) + " " + str(self.hit_end) + " " + str(self.region)

    def get_log(self):
        return self.name + " " + str(self.hit_start) + " " + str(self.hit_end) + " " + self.region + "\n"

    def get_result(self):
        return self.name + " " + self.region + "\n"


class Protein:
    def __init__(self, id, hit_start, hit_end):
        self.result = []
        self.id = id
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.data_from_protein()

    def data_from_protein(self):
        Entrez.email = email  # Always tell NCBI who you are
        handle = Entrez.efetch(db="protein", id=self.id, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        for rec in record[0]["GBSeq_feature-table"]:
            key = rec["GBFeature_key"]
            l_from = []
            l_to = []

            for reg in rec["GBFeature_quals"]:
                region = reg["GBQualifier_value"]
                break

            for location in rec["GBFeature_intervals"]:
                if "other" not in region and region != "active":
                    l_from.append(location["GBInterval_from"])
                    l_to.append(location["GBInterval_to"])

            for x in range(0, len(l_from)):
                if check(self.hit_start, self.hit_end, l_from[x], l_to[x]):
                    if key != "source" and key != "CDS":
                        res = Result(key, l_from[x], l_to[x], region)
                        self.result.append(res)

    def get_result_number(self):
        return len(self.result)

    def show_all(self):
        print("CDS", self.id, "new ->", self.hit_start, self.hit_end)
        for res in self.result:
            if isinstance(res, Result):
                print("\t", res.show_result())

    def get_log(self):
        result = "CDS " + self.id + str(self.hit_start) + " " + str(self.hit_end) + "\n"
        for res in self.result:
            if isinstance(res, Result):
                result += "\t" + res.get_log()

        return result

    def get_result(self):
        result = ""
        count = 0
        for res in self.result:
            if isinstance(res, Result):
                result += "CDS " + res.get_result()
                if res.name != "Protein":
                    count += 1

        if count != 0:
            return result
        else:
            for res in self.result:
                if isinstance(res, Result):
                    if str(res.region).find("polyprotein") != -1:
                        result += "region not known\n"

        return result


class Nucleotide:
    def __init__(self, id, hit_start, hit_end):
        self.x = []
        self.id = id
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.data_from_ncbi()

    def data_from_ncbi(self):
        Entrez.email = email  # Always tell NCBI who you are
        handle = Entrez.efetch(db="nucleotide", id=self.id, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        for rec in record[0]["GBSeq_feature-table"]:
            key = rec["GBFeature_key"]
            l_from = []
            l_to = []
            protein_id = ""

            try:
                for location in rec["GBFeature_intervals"]:
                    l_from.append(location["GBInterval_from"])
                    l_to.append(location["GBInterval_to"])
            except:
                print(self.id)
                print(rec)
                print(sys.exc_info()[0])
                input("Press Enter to continue...")

            for r in range(0, len(l_from)):
                if int(l_from[r]) > int(l_to[r]):
                    x = l_to
                    l_to = l_from
                    l_from = x
                    break

            if key == "CDS":
                for test in rec["GBFeature_quals"]:
                    if test["GBQualifier_name"] == "protein_id":
                        protein_id = test["GBQualifier_value"]

            for x in range(0, len(l_from)):
                if key == "CDS":
                    seq_len = int((self.hit_end - self.hit_start) / 3)
                    new_start = int((self.hit_start - int(l_from[x])) / 3)
                    new_end = new_start + seq_len

                    if new_start < 0 and new_end > 0:
                        new_start = 1

                    if new_start > 0 and new_end > 0:
                        protein = Protein(protein_id, new_start, new_end)
                        self.x.append(protein)

                    break

                else:
                    if check(self.hit_start, self.hit_end, l_from[x], l_to[x]):
                        if key != "source":
                            r = Result(key, l_from[x], l_to[x], "")
                            self.x.append(r)
                            break

    def show_all(self):
        print("Data:", self.id, self.hit_start, self.hit_end)
        for res in self.x:
            if isinstance(res, Result):
                print(res.show_result())

            if isinstance(res, Protein):
                if res.get_result_number() > 0:
                    res.show_all()
        print()

    def get_log(self):
        result = "Data:\n" + self.id + " " + str(self.hit_start) + " " + str(self.hit_end) + "\nResults:\n"
        for res in self.x:
            if isinstance(res, Result):
                result += res.get_log()

            if isinstance(res, Protein):
                if res.get_result_number() > 0:
                    result += res.get_log()

        return result

    def get_result(self):
        result = ""
        for res in self.x:
            if isinstance(res, Result):
                if res.name != "gene":
                    result += res.get_result()

            if isinstance(res, Protein):
                if res.get_result_number() > 0:
                    result += res.get_result()

        if result == "":
            return "not known\n"
        else:
            return result


def red_seq(folder):
    file = open(folder + '\\' + 'SeqNumber.txt', 'r')

    for line in file:
        l = line.rstrip()
        all_seq.append(l);


def exist(i, folder):
    if exists(folder + '\\xml\\' + str(i) + 'seq.xml'):
        return True
    return False


def check(start, end, l_from, l_to):
    len = end - start
    if int(l_from) <= start and int(l_to) >= end:
        return True

    if int(l_from) <= start and int(l_to) + len >= end:
        return True

    return False


def calculate(i, folder):
    blast_qresult = SearchIO.read(folder + '\\xml\\' + str(i) + 'seq.xml', 'blast-xml')
    save_log = open(folder + '\\' + 'logResult.txt', 'a')
    save = open(folder + '\\' + 'Result.txt', 'a')

    blast_hsp = blast_qresult[0][0]
    hit_start = blast_hsp.hit_start + 1
    hit_end = blast_hsp.hit_end
    hit_id = blast_hsp.hit_id
    ids = hit_id.split("|")
    id_letter = ids[3]

    r = Nucleotide(id_letter, hit_start, hit_end)
    r.show_all()

    save_log.write(all_seq[i - 1] + '\n')
    save_log.write("========================\n")
    save_log.write(r.get_log())
    save_log.write("========================\n\n")

    number = all_seq[i - 1].split(".")[0]
    save.write('Sequence: ' + str(number) + '.\n')
    save.write(r.get_result() + "\n")

	
if __name__ == "__main__":
	print("\nFinding functions of RNA fragmnts step 3 - Final results")
	folder_name = input("Enter the name of the test folder: ")
	email = input("Enter Your email address: ")

	red_seq(folder_name)
	i = 1

	while (exist(i, folder_name)):
		number = all_seq[i - 1].split(".")[0]
		print('Calculating ' + str(number) + '...\t')
		calculate(i, folder_name)
		print('Done')
		i += 1