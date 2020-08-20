# insertMissingRows.py
# The goal of this script is to insert rows that were omitted from merging VCF files from different strains
# and to rename the strain columns with TaxonIDs

# To run this script, you will need Biopython on your local machine (or wherever you're running this script from)
#		For downloading the package, see http://biopython.org/wiki/Download
# 		For installation instructions, see http://biopython.org/DIST/docs/install/Installation.html
#
# INPUT FILEs: 1) Galaxy merged VCF file, 2) original VCF files to search, fullpath, 3) file with list of missing rows
# OUTPUT: 1) new merged VCF file with taxonIDs
#

from Bio import Entrez
import re
import sys, os
import getopt
import urllib2
import time
from dateutil.parser import parse

## database connection stuff -- for taxon IDs
from sqlalchemy import create_engine, and_
"""from ..models.models import DBSession, Straindbentity, DBentity
from ..data_helpers.data_helpers import get_eco_ids, get_output, SUBMISSION_VERSION

engine = create_engine(os.getenv('SQLALCHEMY_PROD_DB_URI'), pool_recycle=3600)
DBSession.configure(bind=engine)
"""

###############
## VARIABLES ##
###############
STRAINS = {
    'BY4742': 'S288C',
    'CEN.PK': 'CEN.PK',
    'CEN.PK2': 'CEN.PK',
    'D273-10B': 'D273-10B',
    'FL100': 'FL100',
    'JK9-3d': 'JK9-3d',
    'JK9': 'JK9-3d',
    'RM11-1A': 'RM11-1a',
    'RM11_1A': 'RM11-1a',
    'SEY6210': 'SEY6210',
    'SEY': 'SEY6210',
    '10560_6B': 'Sigma1278b',
    'SK1': 'SK1',
    'W303': 'W303',
    'W303-1a': 'W303',
    'W303-1A': 'W303',
    'W303-K6001': 'W303',
    'W303-1b': 'W303',
    'W303-1B': 'W303',
    'X2180': 'X2180-1A',
    'Y55': 'Y55'
}
CHR_ORD = ('chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII',
           'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV',
           'chrXV', 'chrXVI', 'chrM')

STRAINSTOTAXONS = {
    'BY4742': 'TAX:559292',
    'CEN.PK': 'TAX:889517',
    'D273': 'NTR:101',
    'FL100': 'TAX:947036',
    'JK9': 'NTR:104',
    'RM11_1A': 'TAX:285006',
    'SEY6210': 'NTR:107',
    'SEY': 'NTR:107',
    '10560_6B': 'TAX:658763',
    'SK1': 'TAX:580239',
    'SIGMA': 'TAX:658763',
    'W303': 'TAX:580240',
    'X2180': 'NTR:108',
    'Y55': 'NTR:112'
}

# BY4742	CEN.PK	D273	FL100	JK9	RM11_1A	SEY6210	SK1	SIGMA	W303	X2180	Y55
#
# #def getStrains(strainslist): ## convert strains to NCBI Taxon ID


def processInfile(filename):
    vcfDataObj = dict()
    strainsInVcf = list()
    colHeaders = list()

    with open(filename, 'r') as reader:
        for row in reader.readlines():
            ## CHECK FOR comments/headers at the beginning. Add to 'header'
            if re.match("^#", row):
                ## DO SOMETHING SPECIAL IF IT has the headers with strains ##
                if re.match("^#CHROM", row):
                    print "HEADER ROW: " + row
                    colHeaders = row.split("\t")
                    numCols = len(colHeaders)

                    if numCols > 9:
                        for each in colHeaders[9:numCols - 1]:
                            # print each + "==" + each.rstrip()
                            strainsInVcf.append(each.rstrip())

                if vcfDataObj.has_key('header'):
                    vcfDataObj['header'].append(row.rstrip())
                else:
                    vcfDataObj['header'] = [row.rstrip()]

            ## NON-COMMENT/HEADER rows ##
            ## MAKE AN OBJ -- CHROMOSOME, POS, REF - '
            else:
                dataRow = row.split("\t")
                tempObj = dict()

                for h in colHeaders:
                    #                    print "."
                    #                    print str(colHeaders.index(h)) + "." + h.strip(
                    #                    ) + ": " + dataRow[colHeaders.index(h)]
                    tempObj[h.rstrip()] = dataRow[colHeaders.index(h)].rstrip()
                    # print h.strip() + " indexed at " + str(colHeaders.index(h))

                    # chromosome, position, ref, alt,
                    chr = dataRow[0]
                    pos = int(dataRow[1])
                    ref = dataRow[3]
                    alt = dataRow[4]

                    if vcfDataObj.has_key(chr):
                        if vcfDataObj[chr].has_key(pos):
                            if vcfDataObj[chr][pos].has_key(ref):
                                if vcfDataObj[chr][pos][ref].has_key(alt):
                                    continue  #already there
                                else:
                                    vcfDataObj[chr][pos][ref][alt] = tempObj
                            else:
                                vcfDataObj[chr][pos][ref] = dict()
                                vcfDataObj[chr][pos][ref][alt] = tempObj
                        else:
                            vcfDataObj[chr][pos] = dict()
                            vcfDataObj[chr][pos][ref] = dict()
                            vcfDataObj[chr][pos][ref][alt] = tempObj

                    else:
                        print chr + ":" + str(pos)
                        # for pair in tempObj.items():
                        #     print(pair)
                        #  sys.exit()
                        vcfDataObj[chr] = dict()
                        vcfDataObj[chr][pos] = dict()
                        vcfDataObj[chr][pos][ref] = dict()
                        vcfDataObj[chr][pos][ref][alt] = tempObj

    return (vcfDataObj, colHeaders, strainsInVcf)


def getMissingRows(missingrowsfile, dirToSearch, headers):
    missingDataObj = dict()
    # open file
    # for each line, split and then grep #
    # make an object of each row returning # #associate with the strain/file #

    with open(missingrowsfile, 'r') as reader:
        for row in reader.readlines():
            searchTerms = row.split("\t")
            while ("" in searchTerms):
                searchTerms.remove("")

            print "|".join(searchTerms)

            command = "grep '^" + searchTerms[
                0] + "\t' " + dirToSearch + "*.vcf | grep '\t" + searchTerms[
                    1] + "\t' | grep '\t" + searchTerms[3] + "\t'"

            #print command
            missingLines = os.popen(command)  ## list

            for each in missingLines:
                #      print each

                oriVcf, dataRow = each.split(".vcf")

                data = dataRow.split("\t")  # get datarow

                headcommand = "grep '#CHROM' " + oriVcf + ".vcf"  # grep the column headers from the file
                headrow = os.popen(headcommand)

                for one in headrow:
                    fileheadrow = one.split("\t")

                tempObj = dict()

                for col in fileheadrow:
                    #                    print str(fileheadrow.index(col)) + ": " + col
                    # making a temp obj of the data row

                    tempObj[col.rstrip()] = data[fileheadrow.index(
                        col)].lstrip(":").rstrip()

                chrom = data[0].lstrip(":")
                pos = int(data[1])
                ref = data[3]
                alt = data[4]
                strainnum = len(data) - 1
                strbkd = fileheadrow[strainnum].rstrip()

                # check if chromosome key
                if missingDataObj.has_key(chrom):
                    if missingDataObj[chrom].has_key(pos):  #check pos
                        if missingDataObj[chrom][pos].has_key(
                                ref):  # check ref
                            if missingDataObj[chrom][pos][ref].has_key(
                                    alt):  #check alt
                                if missingDataObj[chrom][pos][ref][
                                        alt].has_key(strbkd):
                                    continue  # skip if strain already there
                                else:  # just add strain
                                    missingDataObj[chrom][pos][ref][alt][
                                        strbkd] = data[len(data) - 1].rstrip()
                            else:  # no strain, no alt
                                missingDataObj[chrom][pos][ref][alt] = tempObj
                        else:  # no strain, no alt, no ref
                            missingDataObj[chrom][pos][ref] = dict()
                            missingDataObj[chrom][pos][ref][alt] = tempObj
                    else:
                        # no strain, no alt, no ref, no pos
                        missingDataObj[chrom][pos] = dict()
                        missingDataObj[chrom][pos][ref] = dict()
                        missingDataObj[chrom][pos][ref][alt] = tempObj
                else:  # no chromosome, pos, ref, alt, strain
                    missingDataObj[chrom] = dict()
                    missingDataObj[chrom][pos] = dict()
                    missingDataObj[chrom][pos][ref] = dict()
                    missingDataObj[chrom][pos][ref][alt] = tempObj

    return missingDataObj


def merge_data(vcfData, missing, strains):
    mergeDataObj = vcfData

    for chrom in missing.keys():
        for pos in missing[chrom].keys():
            for ref in missing[chrom][pos].keys():
                for alt in missing[chrom][pos][ref].keys():
                    newRowObj = missing[chrom][pos][ref][alt]

                    for one in strains:  # adding '.' for all strains without info
                        if newRowObj.has_key(one.rstrip()):
                            #  print 'strain ' + one + " already for " + ",".join(
                            #       [chrom, pos, ref, alt])
                            continue
                        else:
                            #     if one == 'BY4742':
                            #        print 'missing: ' + chrom + ',' + pos + ',' + ref + ',' + alt + ' adding . for ' + one
                            newRowObj[one.rstrip()] = "."

                    if mergeDataObj[chrom].has_key(pos):
                        if mergeDataObj[chrom][pos].has_key(ref):
                            if mergeDataObj[chrom][pos][ref].has_key(alt):
                                continue
                            else:
                                mergeDataObj[chrom][pos][ref][alt] = newRowObj
                        else:
                            mergeDataObj[chrom][pos][ref] = dict()
                            mergeDataObj[chrom][pos][ref][alt] = newRowObj
                    else:
                        mergeDataObj[chrom][pos] = dict()
                        mergeDataObj[chrom][pos][ref] = dict()
                        mergeDataObj[chrom][pos][ref][alt] = newRowObj

    return (mergeDataObj)


def write_lines(fileout, wlines):
    print 'printing to:' + fileout
    writer = open(fileout, 'w')
    # writer.write(u'')
    for line in wlines:
        #     print line
        if re.match(
                "^#CHROM", line
        ):  ##IF it is the header, replace strain names with taxon IDs
            newheader = list()
            headerrow = line.split("\t")
            for col in headerrow:
                # print col + " for header"
                if STRAINSTOTAXONS.has_key(col.rstrip()):
                    #      print "taxon:" + STRAINSTOTAXONS[col.rstrip()]
                    newheader.append(STRAINSTOTAXONS[col.rstrip()])
                else:
                    newheader.append(col.rstrip())
            writer.write("\t".join(newheader) + "\n")
        else:
            writer.write(line + '\n')

    writer.close()


def main(argv):
    if len(argv) > 2:
        opts, args = getopt.getopt(argv, 'i:m:l:o:',
                                   ['inputfile=', 'missing=', 'list=', 'out='])
        print 'getting opts' + str(argv)

        #   outpath = './'

        for opt, arg in opts:
            arg = re.sub('=', '', arg)
            print opt + ":" + arg

            if opt in ("-i", "--inputfile"):
                inputfile = arg
            if opt in ("-l", "--list"):
                searchdir = arg
            if opt in ("-m", "--missing"):
                missingRows = arg
            if opt in ("-o", "--out"):
                outfile = arg

    ### here's where we run the functions
        headerList = [
            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
            'FORMAT', 'BY4742', 'CEN.PK', 'D273', 'FL100', 'JK9', 'RM11_1A',
            'SEY6210', 'SK1', 'SIGMA', 'W303', 'X2180', 'Y55'
        ]

        (data, headerList, strains) = processInfile(inputfile)
        #  print "strains:" + "*".join(strains)
        missingdata = getMissingRows(missingRows, searchdir, headerList)

        alldata = merge_data(data, missingdata, strains)

        print outfile

        write_lines(outfile, data['header'])

        writer = open(outfile, 'a')

        for each in CHR_ORD:
            #
            #    positionsList = alldata[each].keys()
            #    print(positionsList)

            sortedpos = sorted(alldata[each].keys())
            #  print(sortedpos)
            #sys.exit()

            for pos in sortedpos:
                for ref in alldata[each][pos].keys():
                    for alt in alldata[each][pos][ref].keys():
                        rowdata = alldata[each][pos][ref][alt]
                        tempRow = list()
                        for col in headerList:
                            if rowdata.has_key(col.rstrip()):
                                #print col.rstrip() + ":" + rowdata[col.rstrip()]
                                if (rowdata[col.rstrip()] == 'chrM'):
                                    tempRow.append('chrMt')
                                else:
                                    tempRow.append(rowdata[col.rstrip()])
                            else:
                                tempRow.append(".")

                        writer.write("\t".join(tempRow) + "\n")
        #  if each != 'chrMt':
        #      writer.write("\n")  #newline at end of each chromosome
        writer.close()

    #  if (len(result) > 0):
    #output_obj = get_output(result)

    #         file_name = 'src/data_dump/SGD' + SUBMISSION_VERSION + 'basicGeneInformation.json'
    #         json_file_str = os.path.join(root_path, file_name)
    #          with open(json_file_str, 'w+') as res_file:
    #              res_file.write(json.dumps(output_obj))
    else:
        print 'Usage: insertMissingRows.py -i <input file> -l <where to find original VCF files> -m <file with missing rows> -o <outfile path>'
        sys.exit(2)


'''except getopt.GetoptError, err:
        print str(err)
        print 'Usage: insertMissingRows.py -i <input file> -l <list of original VCFs> -m <file with missing rows> -o <outfile path>'
        sys.exit(2)
'''

if __name__ == "__main__":
    main(sys.argv[1:])
# main()
