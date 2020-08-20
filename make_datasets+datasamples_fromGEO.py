# The goal of this script is to generate the three tables of metadata for our expression data directly from GEO...(without use of Janos's script)
# 	-File
#	-Dataset
#	-DataSample
#
#
# To run this script, you will need Biopython on your local machine (or wherever you're running this script from)
#		For downloading the package, see http://biopython.org/wiki/Download
# 		For installation instructions, see http://biopython.org/DIST/docs/install/Installation.html
#
# INPUT: two dates--the range of PDATs you want GEO to look in for yeast datasets
# OUTPUT: two tsv files--the DATASET table and the DATASAMPLE table

from Bio import Entrez
import re
import sys
import getopt
import urllib2
import time
from dateutil.parser import parse

# pass in e-mail as argument
#Entrez.email = raw_input('Please enter you e-mail address:') # sys.argv[1] #
MAX_SEARCH = 5000
MAX_SUMMARIES = 10
WAIT_TIME = 4


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# 												VARIABLES
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# =======
# STRAINS
# =======

strains = {
    'S288C': 'S288C',		'S288c': 'S288C',		's288c': 'S288C',
    'BY4741': 'S288C',		'BY4742': 'S288C',		'BY4743': 'S288C',
    'BY4744': 'S288C',		'BY4745': 'S288C',	's288C': 'S288C',
    'DBY12020': 'S288C',		'DBY12021': 'S288C',
    'CEN.PK': 'CEN.PK',		'CEN.PK2': 'CEN.PK',
    'D273-10B': 'D273-10B',	'FL100': 'FL100',		'JK9-3d': 'JK9-3d',
    'RM11-1a': 'RM11-1a',	'SEY6210': 'SEY6210',	'Sigma1278b': 'Sigma1278b',
    'SK1': 'SK1',
    'W303': 'W303',			'W303-1a': 'W303',			'W303-1A': 'W303',
    'W303-K6001': 'W303',	'W303-1b': 'W303',			'W303-1B': 'W303',
    'X2180-1A': 'X2180-1A',	'Y55': 'Y55'}


# ================
# OBI TERMS--ASSAY
# ================


geo2obi_assay = {
    'Other': 'OBI:0000070',
    'Genome binding/occupancy profiling by array': 'OBI:0001954',
    'Genome binding/occupancy profiling by genome tiling array': 'OBI:0001419',
    'Genome binding/occupancy profiling by high throughput sequencing': 'OBI:0000716',
    'Genome variation profiling by array': 'OBI:0001393',
    'Genome variation profiling by genome tiling array': 'OBI:0002030',
    # 'Genome variation profiling by genome tiling array':'OBI:0001274',
    'Genome variation profiling by SNP array': 'OBI:0002031',
    'Expression profiling by array': 'OBI:0001463',
    'Expression profiling by genome tiling array': 'OBI:0001235',
    'Expression profiling by high throughput sequencing': 'OBI:0001271',
    'Methylation profiling by high throughput sequencing': 'OBI:0001266',
    'Non-coding RNA profiling by high throughput sequencing': 'OBI:0001271',
    'Third-party reanalysis': 'OBI:0200000'}


obi2name_assay = {
    'NTR(GSE14217)': 'synthetic lethality analysis with microarrays',
    'NTR(GSE1753)': 'synthetic lethality analysis with microarrays',
    'NTR(GSE1838)': 'synthetic lethality analysis with microarrays',
    'NTR(GSE20961)': 'competitive growth assay analysis with microarrays',
    'NTR(GSE46145)': 'competitive growth assay analysis with microarrays',
    'OBI:0000070': 'assay',
    'OBI:0000438': 'DNA sequence variation detection',
    'OBI:0000626': 'DNA sequencing',
    'OBI:0000716': 'ChIP-seq assay',
    'OBI:0001235': 'transcription profiling by tiling array assay',
    'OBI:0001247': 'genotyping by high throughput sequencing assay',
    'OBI:0001248': 'ChIP-chip assay',
    'OBI:0001266': 'DNA methylation profiling by high throughput sequencing assay',
    'OBI:0001271': 'RNA-seq assay',
    'OBI:0001274': 'genotyping by array assay',
    'OBI:0001318': 'proteomic profiling by array assay',
    'OBI:0001332': 'DNA methylation profiling by array assay',
    'OBI:0001393': 'comparative genomic hybridization by array assay',
    'OBI:0001419': 'ChIP-chip by tiling array assay',
    'OBI:0001463': 'transcription profiling by array assay',
    'OBI:0001858': 'cross-linking immunoprecipitation high-throughput sequencing assay',
    'OBI:0001924': 'micrococcal nuclease digestion followed by high throughput sequencing assay',
    'OBI:0001954': 'ChIP assay',
    'OBI:0002030': 'genotyping by tiling array',
    'OBI:0002031': 'genotyping by SNP array',
    'OBI:0002036': 'array based nucleic acid structure mapping assay',
    'OBI:0200000': 'data transformation'}

# ====================
# OBI TERMS--BIOSAMPLE
# ====================

geo2obi_biosample = {
    'nuclear RNA': 'OBI:0000862',
    'nuclear RNA extract': 'OBI:0000862',
    'polyA RNA': 'OBI:0000869',
    'polyA RNA extract': 'OBI:0000869',
    'cytoplasmic RNA': 'OBI:0000876',
    'cytoplasmic RNA extract': 'OBI:0000876',
    'protein extract': 'OBI:0000894',
    'total RNA': 'OBI:0000895',
    'total RNA extract': 'OBI:0000895',
    'genomic DNA': 'OBI:0001051',
    'DNA extract': 'OBI:0001051',
    'other': 'OBI:0000423'}		# this is handling all the "else" cases

obi2name_biosample = {
    'OBI:0000862': 'nuclear RNA extract',
    'OBI:0000869': 'polyA RNA extract',
    'OBI:0000876': 'cytoplasmic RNA extract',
    'OBI:0000894': 'protein extract',
    'OBI:0000895': 'total RNA extract',
    'OBI:0001051': 'DNA extract',
    'OBI:0000423': 'extract'}


# ======
# TAXONS
# ======

yeasts = {
    'Saccharomyces cerevisiae': '4932',			'Saccharomyces bayanus': '4931',
    'Saccharomyces kluyveryii': '4934',			'Saccharomyces kudriavzevii': '114524',
    'Saccharomyces mikatae': '114525',			'Saccharomyces paradoxus': '27291',
    'Saccharomyces uvarum': '230603',			'Saccharomyces castellii': '27288',
    'Candida albicans': '5476',					'[Candida] glabrata': '5478',
    'Debaryomyces hansenii': '4959',				'Kluyveromyces lactis': '28985',
    'Kluyveromyces waltii': '4914',				'Lachancea kluyveri': '4934',
    'Lachancea waltii': '4914',					'Naumovozyma castellii': '27288',
    'Schizosaccharomyces japonicus': '4897',		'Schizosaccharomyces octosporus': '4899',
    'Schizosaccharomyces pombe': '4896',			'Yarrowia lipolytica': '4952'}

# strains for saccharomyces cerevisiae all map to 4932
sc_strains = {
    'Saccharomyces cerevisiae S288c': '4932',	'Saccharomyces cerevisiae BY4741': '4932',
    'Saccharomyces cerevisiae S288C': '4932',	'Saccharomyces cerevisiae SK1': '4932',
    'Saccharomyces cerevisiae Sigma1278b': '4932', 'Saccharomyces cerevisiae W303': '4932'}

# yeast crosses--we can't just search for "Saccharomyces cerevisiae" since it'll include some of these guys too...
yeast_crosses = {
    'Saccharomyces cerevisiae x Saccharomyces kudriavzevii': '332112',
    'Saccharomyces cerevisiae x Saccharomyces paradoxus': '595493',
    'Saccharomyces cerevisiae x Saccharomyces uvarum': '489137',
    'Saccharomyces paradoxus x Saccharomyces uvarum': '1911234'}

# non-yeast orgs
various_orgs = {
    'Homo sapiens': '9606',						'Mus musculus': '10090',
    'Caenorhabditis elegans': '6239',			'Drosophila melanogaster': '7227',
    'Neurospora crassa': '5141',					'Escherichia coli': '562',
    'Pseudomonas putida': '303',					'Synechococcus sp. WH 8102': '84588'}

# viruses--we can't just search for "Saccharomyces cerevisiae" since it'll include some of these guys too...
viruses = {
    'Lachancea kluyveri NRRL Y-12651': '4934',
    'Saccharomyces cerevisiae killer virus M1': '12450',
    'Saccharomyces 20S RNA narnavirus': '186772',
    'Saccharomyces 23S RNA narnavirus': '198599',
    'Saccharomyces cerevisiae virus L-A': '11008',
    'Saccharomyces cerevisiae virus L-BC (La)': '42478'}

# strains of other orgs and other weird cases
iffy_ones = {
    'Neurospora crassa OR74A': '5141',
    'synthetic construct': '32630',
    'Synthetic plasmid': '1440148',
    'Cryptococcus neoformans var. grubii': '178876'}

taxons2id = {}
taxons2id.update(yeasts)
taxons2id.update(sc_strains)
taxons2id.update(yeast_crosses)
taxons2id.update(various_orgs)
taxons2id.update(viruses)
taxons2id.update(iffy_ones)


# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# 												METHODS
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# generate_dataset_datasample_tables()
# build_series()
# build_samples()


# generate Dataset and DataSample tables directly from GEO...(without use of Janos's script)
# @param date_start Search GEO records from after this month. Format YYYY/MM.
# @param date_end Search GEO records from before this month. Format YYYY/MM.

def generate_dataset_datasample_tables(date_start, date_end, ext_str, path):

    gse2dataset_lns = {}
    gsm2datasample_lns = {}

    # Build Search with the following criteria:
    #	-TYPE is gse
    #	-ORGANISM is Saccharomyces cerevisiae
    #	-PUBLIC DATE is between 'date_start' and  'date_end'
    gse_term = 'Saccharomyces cerevisiae[ORGN]&%s:%s[PDAT]&gse[ETYP]' % (
        date_start, date_end)
    gse_idlist = search_gds(gse_term)
    print gse_idlist

    # Retrieve Summaries for GSE search results ids
    gse_summaries = retrieve_summaries(gse_idlist)

    # For each GSE series object....
    for gse_summary in gse_summaries:

        print '<><><><><>%s<><><><><>' % ('GSE'+gse_summary['GSE'])
        print u'(%s)' % gse_summary[u'Id']
# 		for k,v in gse_summary.items():
# 			print (k,v)

        # Build Search to get GSM-typed objects associated with current GSE accn
        gsm_term = '(%s[ACCN]&gsm[ETYP])' % ('GSE'+gse_summary['GSE'])
        gsm_idlist = search_gds(gsm_term)

        # CHECK that number of samples match metadata
        gse_accn = 'GSE%s' % gse_summary['GSE']
        if len(gsm_idlist) != gse_summary['n_samples']:
            print 'CURATOR ACTION: check samples--number of samples found with Entrez search does not match GEO-recorded number of samples for series (%s)' % gse_accn
            print '\tEntrezSearch=%d, GEOrecord=%d' % (
                len(gsm_idlist), gse_summary['n_samples'])

        # Retrieve Summaries for GSM search results ids
        gsm_summaries = retrieve_summaries(gsm_idlist)

        # Add GSM entries to the DataSample Table
        gsm_lns, gse_channel_ct = build_samples(gsm_summaries, gse_summary)

        # Add the GSE string for the DataSet Table
        gse_dict = build_series(gse_summary, gse_channel_ct)

        gse2dataset_lns.update(gse_dict)

        # Don't add to sample table if this series is a superseries
        if is_superseries(gse_summary) == 1:
            continue
        gsm2datasample_lns.update({gse_accn: gsm_lns})

    # Add parent series to dataset_lns
    gse2dataset_lns = add_parent_series(gse2dataset_lns)

    # Dict of lines to List and sort
    dataset_lns = []
    datasample_lns = []
    sorted_gse = gse2dataset_lns.keys()
    sorted_gsm = gsm2datasample_lns.keys()

    for gse in sorted_gse:
        dataset_lns.append(gse2dataset_lns[gse])
    for gsm in sorted_gsm:
        datasample_lns.append(gsm2datasample_lns[gsm])

    # issues b/c unicode formatting
# 	dataset_lns = gse2dataset_lns.values()
# 	datasample_lns = gsm2datasample_lns.values()
# 	dataset_lns.sort()
# 	datasample_lns.sort()

    # Write lines to file
    write_lines(path+'/dataset'+ext_str.replace("/", "-")+'.tsv', dataset_lns)
    write_lines(path+'/datasample'+ext_str.replace("/", "-") +
                '.tsv', datasample_lns)

    #print_tables(dataset_lns, datasample_lns)


# ======
# SERIES
# ======

# Collect GSEs and GSMs separately...
#	-add GSEs to Dataset lists
#	-check GSM list for any samples related to GSEs not in first list
#	-print the missed GSEs and add to gse_list...


# This bundles the data collection for the DataSeries table--given a GSE summary, return the formatted string for adding to the DataSeries table
# @param	gse_summary	dictionary	Entrez summary for a given GSE object
# @return:	gse_ln		string		Formatted row for Dataset Table
# called in the generate_dataset_datasample_tables() method
def build_series(gse_summary, channel_ct):  # (gse_id)
    global taxons2id
    gse_ln = ''

    # used for pmids, lab_name, and lab_location
    pmid_lst = [str(i) for i in gse_summary['PubMedIds']]

    # DataSample Table Info
    gse_accn = 'GSE%s' % gse_summary['GSE']
    title = gse_summary['title']  # unicode line
    source = 'GEO'
    ref_type = 'GEO'
    date_public = parse(gse_summary['PDAT']).strftime('%Y-%m-%d')

    parent_dataset = ''
    if is_superseries(gse_summary) == 1:
        parent_dataset = '*'

    # Use helper method to find appropriate obi term
    assay_name, assay_id = get_assay(gse_summary)

    # Use taxon dictionary to translate taxon names to ids
    t_ids = []
    gse_taxons = gse_summary['taxon'].split(';')
    for t in gse_taxons:
        try:
            id = taxons2id[t.strip()]
            # avoids listing same taxon twice (e.g. "Saccharomyces cerevisiae;Saccharomyces cerevisiae SK1")
            if id not in t_ids:
                t_ids.append(id)
        except KeyError:
            print 'MISSING TAXON: %s' % t.strip()
            print 'CURATOR ACTION: update the taxon dictionary in this script to include the above missing taxon'
            t_ids.append('*')
    taxon_id = '|'.join(t_ids)

    sample_ct = str(gse_summary['n_samples'])
    is_in_spell = '0'
    is_in_browser = '0'
    series_summary = gse_summary['summary']  # unicode line

    # call helper method that retrieves from db=pubmed
    lab_name, lab_location = get_lab_info(pmid_lst, gse_accn)

    # Leave blank--manually filled in by curators
    keywords = ''

    pmids = '|'.join(pmid_lst)

    # Leave blank--manually filled in by curators
    display_name = ''

    # 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s' % gse_accn
    obj_url = gse_summary['FTPLink']
    url_type = 'GEO'

    gse_ln = u'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        gse_accn, title, source, gse_accn, ref_type, date_public, parent_dataset, assay_name, assay_id, taxon_id, channel_ct, sample_ct, is_in_spell, is_in_browser, series_summary, lab_name, lab_location, keywords, pmids, display_name, obj_url, url_type)

    print gse_ln
    return {gse_accn: gse_ln}


# =======
# SAMPLES
# =======


# This bundles the data collection for the DataSamples table into strings formatted according to specifications
# @param:	gsm_summaries	listofdictionaries	List of Entrez summaries for the samples of a given GSE object
# @param:	gse_summary		dictionary			Entrez summaries for the given GSE object
# @return:	gse_lns 		listofstrings		Formatted rows for DataSample Table associated with a given GSE object
# called in the generate_dataset_datasample_tables() method
def build_samples(gsm_summaries, gse_summary):

    gsm2lns = {}
    gsm_lines = ''

    gse_channels = []
    gse_channel_ct = '1'

    # Two ways of assigning sample_order:
    #		A) by the order gse_summary returns them
    #		B) by ascending GSM-values
    # Edith has chosen the B) option...
    # first, sort samples into ascending order
    sample_list = []
    for sample in gse_summary['Samples']:
        sample_list.append(sample['Accession'][3:])
    sample_list.sort()

    gse_accn = 'GSE%s' % gse_summary['GSE']
    for gsm in gsm_summaries:
        #print gsm
        # 		for km,vm in gsm.items():
        # 			print (km,vm)

        gsm_display_name = gsm['title']
        gsm_description = gsm['summary']
        gsm_accn = gsm['Accession']
        gsm_ref_type = 'GEO'

        # Use helper method to find appropriate obi term
        gsm_biosample_name, gsm_biosample_id = get_biosample(gsm_accn)

        # STRAIN_INFO
        gsm_strain = get_strain(
            gsm_accn, gsm_description, gse_summary['taxon'])
        #print 'Strain:%s'	% gsm_strain
        # SAMPLE_ORDER--use placeholder for now: '****' will be replaced after loop finishes
        #gsm_sample_order = str( sample_list.index(gsm_accn[3:]) + 1 )
        gsm_sample_order = '****'

        gsm_ln = u'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
            gse_accn, gsm_display_name, gsm_description, gsm_accn, gsm_ref_type, gsm_biosample_name, gsm_biosample_id, gsm_strain, gsm_sample_order)
        gsm2lns.update({gsm_accn: gsm_ln})

        #print '<><><><><><><><><><><><><><><>'
        if 'channel 2' in gsm_description:
            gse_channels.append(gsm_description)

    gsm_accessions = gsm2lns.keys()
    gsm_accessions.sort()
    order_count = 1
    for gsm in gsm_accessions:
        gsm_lines = gsm_lines + \
            gsm2lns[gsm].replace('\t****\n', '\t%d\n' % order_count)
        order_count += 1

    print gsm_lines

    # return 1 or 2 based on channel count
    if len(gse_channels) == len(gsm_summaries):
        print 'all have channel=2'
        return gsm_lines, '2'
    elif len(gse_channels) == 0:
        print 'all have channel=1'
        return gsm_lines, '1'

    # if conflicting 1 and 2 within a series, we will check if it's a superseries (just set as 2), otherwise tag for curator review
    if is_superseries(gse_summary) == 1:
        return gsm_lines, '2'

    print 'CURATOR REVIEW: this series may contain samples with differing numbers of channels (%s)' % gse_accn
    print '%d of %d' % (len(gse_channels), len(gsm_summaries))

    return gsm_lines, 'CURATOR REVIEW'


# =======
# HELPERS
# =======

# write_lines()
# print_tables()

# Helper method for writing lines to a file. Formats to unicode so we can use the special characters in our database.
# @param	fileout	string			filename to write strings to
# @param	wlines	listofStrings	strings to write to the file
# @return			void
# called by generate_dataset_datasample_tables
def write_lines(fileout, wlines):
    writer = open(fileout, 'w')
    writer.write(u'')
    for line in wlines:
        print line
        writer.write(line.encode('utf8'))
    writer.close()

# De-bugging/display purposes
# Print results so I don't have to open the file to see obvious errors...


def print_tables(dataset_lns, datasample_lns):
    print '<><><><><><><><><><> DATASET TABLE <><><><><><><><><><>'
    for l in dataset_lns:
        print l
    print '<><><><><><><><><> DATASAMPLE TABLE <><><><><><><><><><>'
    for l in datasample_lns:
        print l


# ==============
# ENTREZ HELPERS
# ==============


# search_gds()
# retrieve_summaries()

# Retrieve the uids for the search term input for use with the retrieve_summaries method (see below).
# @param: search_term string Entrez search string used
# @return: idlist list Uids from results under the 'IdList' field
def search_gds(search_term):
    search_handle = Entrez.esearch(
        db='gds', retmax=MAX_SEARCH, term=search_term)
    idlist = Entrez.read(search_handle)['IdList']
    print '# results found:%d' % len(idlist)

    # Throw Warning if number of results meets the default cap.
    if len(idlist) == MAX_SEARCH:
        print 'Warning: search results max reached. possibly more results not captured by this script.\n\t Try running with a smaller date range'
    return idlist

# Formats idlist and runs eutils to retrieve object summaries, only the designated maximum at a time (MAX_SUMMARIES). Slowest part about this script. This runs very slow because of the frequent programmed delays but this method prevents the api from throwing ""HTTP Error 503: Service Temporarily Unavailable" errors.
# @param: idlist list Uids from search results
# @return: summaries list Summary objects returned by Entrez


def retrieve_summaries(idlist):
    summaries = []
    idlist_str = [str(i) for i in idlist]
    idlist_len = len(idlist_str)
    iterations = len(idlist)/MAX_SUMMARIES + 1

    for i in range(iterations):
        end = idlist_len-(i*MAX_SUMMARIES)
        beg = idlist_len-(i*MAX_SUMMARIES)-MAX_SUMMARIES
        if beg < 0:
            beg = 0

        idstring = ','.join(idlist_str[beg:end])
        if idstring == '':
            continue

        try:
            summary_handle = Entrez.esummary(db='gds', id=idstring)
            summaries.extend(Entrez.read(summary_handle))
        except urllib2.HTTPError:
            print "Error fetching", idlist
            # we have angered the API! Try waiting longer?
            time.sleep(WAIT_TIME)
            try:
                summary_handle = Entrez.esummary(db='gds', id=idstring)
                summaries.extend(Entrez.read(summary_handle))
            except:
                print "Error fetching(second try)", idlist
                # we have angered the API! Try waiting longer?
                time.sleep(WAIT_TIME)
                try:
                    summary_handle = Entrez.esummary(db='gds', id=idstring)
                    summaries.extend(Entrez.read(summary_handle))
                except:
                    print "Error fetching(third try)", idlist
                    # wait and try once more...
                    time.sleep(WAIT_TIME)
                    summary_handle = Entrez.esummary(db='gds', id=idstring)
                    summaries.extend(Entrez.read(summary_handle))
        time.sleep(WAIT_TIME)

    return summaries

# @param	string	GSE or GSM accession
# @return	string 	HTML text from the GEO page


def get_html(geo_accn):
    opened_url = urllib2.urlopen(
        'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + geo_accn)
    html_txt = opened_url.read()
    return html_txt

# =============
# "GET" HELPERS
# =============

# get_lab_info()
# get_assay()
# get_biosample()
# get_strain()
# is_superseries()


# Retrieve the information fields that can be found in PubMed.
# @param pmid_lst list All pmids associated with a GSE summary entry cast as string-types
# @return lab_name str LastAuthor of PubMed citation's AuthorList. Field in Dataset Table.
# @return lab_location str Affiliation of the LastAuthor object in the PubMed citation. Field in Dataset Table.
# Important notes:
#		-IGNORING MULTIPLE PMIDS FOR NOW...for GSE with multiple PMIDs, we are just taking the lab location and lab name info from one of them.
#		-for lab name, I am assuming last name in author list of efetch results is the same as last author results for esummary
#		-Alternative for those without pmids: scrape html for 'Contributors' section...see https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76117
#		-Taking first affiliation listed for now...
#		-Formatting of lab_location is still a little raw and may require some further adjustment/trimming...
def get_lab_info(pmid_lst, gse_accn):

    lab_name = u''
    lab_location = u''
    # if no publication associated, can't fill in lab_location and lab_name
    if len(pmid_lst) == 0:
        # return u'',u''
        html_txt = get_html(gse_accn)

        # LAB_LOCATION
        lab_location = re.findall(
            '<tr bgcolor="#eeeeee" valign="top"><td nowrap>Organization name</td>\n<td style="text-align: justify">(.+)<br></td>', html_txt)[0].decode('utf-8')

    else:
        if len(pmid_lst) > 1:
            print 'CURATOR NOTE: multiple PMIDs for this series so this script is arbitrarily using the first one (PMIDs=%s)' % ','.join(
                pmid_lst)

        # grab first pmid
        pmid = pmid_lst[0]
        # get record from PubMed
        summary_handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
        record = Entrez.read(summary_handle)['PubmedArticle'][0]

        # LAB_NAME
        # build lab_name field using LastName and Initials of the last AuthorList entry. Formatted to be consistent with previously collected metadata
        last_name = record['MedlineCitation']['Article']['AuthorList'][-1]['LastName']
# 		print last_name
# 		print type(last_name)
        initials = record['MedlineCitation']['Article']['AuthorList'][-1]['Initials']
# 		print initials
# 		print type(initials)
        lab_name = u'%s %s' % (last_name, initials)
# 		print lab_name
# 		print type(lab_name)

        # CODE: may need to cap string size for lab affiliation.

        # LAB_LOCATION
        #	If no affiliations, set as empty list.
        #	If 1 or greater affiliations, set as first affiliation listed.
        #	If greater than 1 affiliations, also print warning to user about multiple affiliations listed so curator can fix if needed.
        lab_location = u''
        no_affiliations = len(
            record['MedlineCitation']['Article']['AuthorList'][-1]['AffiliationInfo'])
        if no_affiliations >= 1:
            lab_location = record['MedlineCitation']['Article']['AuthorList'][-1]['AffiliationInfo'][0]['Affiliation']
            if no_affiliations > 1:
                print 'CURATOR NOTE: multiple affiliations for this author so this script is arbitrarily using the first one\n\t(PMID=%s)' % pmid
# 	print lab_name
# 	print type(lab_name)
#
# 	print lab_location
# 	print type(lab_location)
#
    return lab_name, lab_location


# GEO to OBI translation for assay information in two columns of the Datasample files
# @param	string	gse_summary				is the GSE summary
# @return	tuple	(obi_names, obi_ids)	should return a tuple
#				obi_names	is a string of |-delimited assay term names and
#				obi_ids		is a string of |-delimited assay term ids
# called in the build_series()
def get_assay(gse_summary):
    geo_ln = gse_summary['gdsType']

    global geo2obi_assay
    global obi2name_assay

    obi_ids = []
    obi_names = []

    # for each geo string that describes an assay (;-delimited)...
    geo_assays = geo_ln.split(';')
    for string in geo_assays:
        str = string.strip()
        # if it is in the geo/obi assay dictionaries...
        if str in geo2obi_assay.keys():
            # translate it to an obi term and add the respective term name and ids to the obi_ids and obi_names lists
            id = geo2obi_assay[str]
            obi_ids.append(id)
            obi_names.append(obi2name_assay[id])
        # if the geo assay string was not recognized, print message for how it should be added to the geo2obi_assay dictionaries at the top of this script
        else:
            print 'CURATOR ACTION: Update GEO assay dictionary to include the following assay--%s' % str
            obi_ids.append('CURATOR REVIEW')
            obi_names.append('CURATOR REVIEW:'+str)

    return '|'.join(obi_names), '|'.join(obi_ids)


# GEO webpage html scraping for biosample information in two columns of the Datasample files
# @param string gsm_accn is the GSM accession string
# @return tuple (obi_name, obi_id) should return a tuple of the biosample term name and id
# called in build_samples()
def get_biosample(gsm_accn):

    global geo2obi_biosample
    global obi2name_biosample

    # retrieve html text from webpage
    html_txt = get_html(gsm_accn)
    # grab 50 characters following "Extracted molecule" which should be where the biosample has been formatted into the html
    #html_chunk = html_txt.split('Extracted molecule</td>\n')[1][:50]

    # Use RegEx to find extracted molecule chunks
    extracted_molecule = ','.join(re.findall(
        '<td nowrap>Extracted molecule</td>\n<td>(.+)</td>', html_txt))

    # store all possible biosamples (strings that match any of the known geo2obi_biosample keys)
    biosamples_column = []
    for mol in extracted_molecule.split(','):
        if mol in geo2obi_biosample.keys():
            biosamples_column.append(geo2obi_biosample[mol])
    # unique the items in biosamples_column
    biosamples_column = list(set(biosamples_column))

    # if one biosample is found, simply return tuple of obi_term name and obi_term id...
    if len(biosamples_column) == 1:
        return obi2name_biosample[biosamples_column[0]], biosamples_column[0]
    # if no biosamples were found, curator may need to add/update the geo2obi_biosample and obi2name_biosample dictionaries at the top of this script
    elif len(biosamples_column) == 0:
        print '\nCURATOR ACTION: Biosample Dictionary in need of update <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
        print extracted_molecule
        return '', ''
    # if  multiple biosample hits, curator will have to go in and figure out what it is and what to do with it...(typo in record or maybe fix script to handle case?)
    else:									# if len(potential_biosamples)>1:
        name = 'CURATOR REVIEW:%s <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<multiple "extracted molecules" found' % (
            ','.join(biosamples_column))
        return name, 'CURATOR REVIEW'


# This helper method uses several approaches to try and determine the appropriate strain to assign the datasample this includes (A)taking a look at the taxon to see if it specifies a specific strain and (B)searching the description for known strain names (based on strains dictionary at the top of the script)
# @param	string		gsm_accn
# @param	string		gsm_description
# @param	string		gse_taxon
# @return	either (1) a string of the strain name, (2) a "CURATOR REVIEW" message with multiple possible strains (for curator to check manually), or (3) a "CURATOR REVIEW" message declaring none found (curator can double-check)
# called in the build_samples() method
def get_strain(gsm_accn, gsm_description, gse_taxon):

    global strains
    gsm_strains = []

    # (A)-prep
    taxons = gse_taxon.split(';')
    # (B)-prep
    description_tokens = gsm_description.replace('_', ' ').split()

    for k in strains.keys():
        # (A) use strain taxon info
        if 'Saccharomyces cerevisiae ' + k in taxons and strains[k] not in gsm_strains:
            gsm_strains.append(strains[k])

        # (B) use gsm esummary
        # this method for identifying strains may require further development...
        if k in description_tokens and strains[k] not in gsm_strains:
            gsm_strains.append(strains[k])

    if len(gsm_strains) == 1:
        return gsm_strains[0]
    # if multiple found
    elif len(gsm_strains) > 1:
        review_str = 'CURATOR REVIEW: multiple found in description summary/gse taxons.--Check if strain is %s' % (
            ' or '.join(gsm_strains))
        print review_str
        return review_str
    # if none found
    return ''


# This helper method just returns a binary 0/1--if the GSE is a superseries or not...
# called in generate_dataset_datasample_tables(), build_series(), and build_samples()
def is_superseries(gse_summary):
    gse_accn = 'GSE%s' % gse_summary['GSE']
    if 'SuperSeries' in gse_summary['summary']:
        # 		print 'CURATOR ACTION: SuperSeries found so manual edits required %s' % (gse_accn)
        return 1
    return 0

    # html---scrape for...
    #		(1) superseries parents
    #		(2) for lab location and name for series w/o pmid associated
    #		(3) improve strain search by adding "characteristics" section (format "strain: %s" % strain_name)
    #


def add_parent_series_via_superseries(dataset_lns):

    parent2sub = {}

    for gse_accn in dataset_lns.keys():
        # retrieve html text from webpage
        html_txt = get_html(gse_accn)

        # Grab chunk of text related to Superseries relationships
        start = html_txt.find(
            '<tr valign="top"><td colspan="2">This SuperSeries is composed of the following SubSeries:</td></tr>')
        end = html_txt.find(
            '<tr valign="top"><td colspan="2"><strong>Relations</strong></td></tr>')
        if start == -1:
            continue
        html_chunk = html_txt[start:end]

        html_chunk.split('geoaxema_recenter)">GSE')[1:]
        print html_chunk

        subseries = list(set(re.findall('GSE[0-9]{3,6}', html_chunk)))
        parent2sub.update({gse_accn: subseries})

    for parent in parent2sub.keys():

        for subseries in parent2sub[parent]:

            dataset_tokens = dataset_lns[subseries].split('\t')
            dataset_tokens[6] = parent
            dataset_lns[subseries] = '\t'.join(dataset_tokens)

    return dataset_lns

# HTML scrape for superseries information and input into appropriate subseries column
# @param dict gse2dataset_lns
# @return dict the updated version with parent series added to the appropriate


def add_parent_series(dataset_lns):
    # Generate subseries to parents dictionary
    sub2parent = {}
    for gse_accn in dataset_lns.keys():
        # retrieve html text from webpage
        html_txt = get_html(gse_accn)

        # Grab chunk of text related to SuperSeries/SubSeries relationships (and ignore if there isn't any)
        html_chunk = ''.join(re.findall(
            '<tr valign="top"><td colspan="2">This SubSeries is part of SuperSeries:</td></tr>(.+)<tr valign="top"><td colspan="2"><strong>Relations</strong></td></tr>', html_txt, re.DOTALL))
        # grab any GSEs within this chunk and assign in sub2parent dictionary
        parents = list(set(re.findall('GSE[0-9]{2,6}', html_chunk)))

        # if series doesn't have a parent, skip
        if len(parents) == 0:
            continue
        # update dictionary of subseries to add parents for
        sub2parent.update({gse_accn: parents})

    # Update the dataset_lns dictionary
    for subseries in sub2parent.keys():
        dataset_tokens = dataset_lns[subseries].split('\t')
        dataset_tokens[6] = '|'.join(sub2parent[subseries])
        dataset_lns[subseries] = '\t'.join(dataset_tokens)

    return dataset_lns

# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


def retrieve_display_names(pmids):

    pmid2display_names = {}
    idstring = ','.join(pmids)

    summary_handle = Entrez.efetch(db='pubmed', id=idstring, retmode='xml')
    records = Entrez.read(summary_handle)['PubmedArticle']

    for record in records:
        pmid = record['MedlineCitation']['PMID']
        author = record['MedlineCitation']['Article']['AuthorList'][0]['LastName']
        year = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']

        display_name = '%s_%s_PMID_%s' % (author, year, pmid)

        pmid2display_names.update({pmid: display_name.decode('utf8')})

    return pmid2display_names


# Based on this test, I have determined that the maximum number of summaries I can retrieve at once using the idlist is 80.
# Thus, I am editing this script so the helper method only retrieves 80 at a time but returns a full list of summaries so that I don't have to adjust the whole script (just the helper)
# Update: the maximum number changes as you use it...I have switched to retrieving 5 or so at a time (see MAX_SUMMARIES constant)
def test_max():
    term = '(GSE71490[ACCN]&gsm[ETYP])'
    idlist = search_gds(term)
    idlist_str = [str(i) for i in idlist]

    cur_len = 60
    for i in range(500):
        print 'length: %d' % cur_len
        idstring = ','.join(idlist_str[:cur_len])
        print idstring
        summary_handle = Entrez.esummary(db='gds', id=idstring)
        summaries = Entrez.read(summary_handle)
        print 'found summaries!'
        cur_len += 1
    return summaries


#esearch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds...'

#efetch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds...'


# def test_range(some_int):
#
# 	lst = []
# 	for i in range(some_int):
# 		lst.append( i )
# 	print 'length: %d\n' % len(lst)
# 	print lst
#
# 	summaries = []
# 	lst = [str(i) for i in lst]
#
# 	iterations = len(lst)/80 + 1
#
# 	lst_len = len(lst)
# 	for i in range(iterations):
# 		print 'i=%d' % i
#
# 		end = lst_len-(i*80)
# 		beg = lst_len-(i*80)-80
#
# 		if beg < 0:
# 			beg = 0
# 		id_lst = lst[beg:end]
#
# 		print id_lst
# 		print len(id_lst)
# 		summaries.extend(id_lst)
# 	print summaries
# 	print len(summaries)

# def test_get_lab_info():
# 	get_lab_info([],'GSE76716')
# 	get_lab_info([],'GSE86279')
# 	get_lab_info([],'GSE76715')
# 	get_lab_info(['26895426'])
    # ...For dealing with case of empty AffiliationInfo list:
# 	get_lab_info(['25425491'])

# def test_add_parent_series():
# 	dataset_lns = {'GSE76778':'0\t1\t2\t3\t4\t5\t6\t7\t8', 'GSE76715':'0\t1\t2\t3\t4\t5\t6\t7\t8',
# 	'GSE76716':'0\t1\t2\t3\t4\t5\t6\t7\t8', 'GSE93959':'0\t1\t2\t3\t4\t5\t6\t7\t8',
# 	'GSE86279':'0\t1\t2\t3\t4\t5\t6\t7\t8', 'GSE84066':'0\t1\t2\t3\t4\t5\t6\t7\t8',
# 	'GSE84065':'0\t1\t2\t3\t4\t5\t6\t7\t8'}
#
# 	dataset_lns = add_parent_series(dataset_lns)
# 	print dataset_lns

def main(argv):

    try:
        if len(argv) > 2:
            opts, args = getopt.getopt(
                argv, 's:l:e:o', ['startdate=', 'lastdate=', 'email=', 'outpath='])
            print 'getting opts' + str(argv)

            outpath = './'

            for opt, arg in opts:
                arg = re.sub('=', '', arg)

                if opt in ("-s", "--startdate"):
                    start = arg
                if opt in ("-l", "--lastdate"):
                    enddate = arg
                if opt in ("-e", "--email"):
                    Entrez.email = arg
                if opt in ("-o", "--outpath"):
                    outpath = arg

            out_ext = start+"to"+enddate

            #Entrez.email = #raw_input('Please enter you e-mail address:') # sys.argv[1] #
        #    print out_ext + "," + Entrez.email + "," + start + \
         #       "," + enddate + "," + start.replace("-", "/"

            generate_dataset_datasample_tables(start.replace(
                "-", "/"), enddate.replace("-", "/"), out_ext, outpath)
            # generate_dataset_datasample_tables(
            #    start, enddate, out_ext, outpath)

        else:
            print 'Usage: make_datasets+datasamples_fromGEO.py -s <YYYY-MM> -l <YYYY-MM> -e <email address> -o <outfile path>'
            sys.exit(2)

    except getopt.GetoptError, err:
        print str(err)
        print 'Usage: make_datasets+datasamples_fromGEO.py -s <YYYY-MM> -l <YYYY-MM> -e <email address> -o <outfile path>'
        sys.exit(2)


# 	t0 = time.time()
# 	generate_dataset_datasample_tables('2017/01','2017/01','2017-01')
# 	t1 = time.time()
# 	print 'TIME(Jan):'
# 	print t1-t0
# 	generate_dataset_datasample_tables('2017/02','2017/02','2017-02')
# 	t2 = time.time()
# 	print 'TIME(Feb):'
# 	print t2-t1
# 	generate_dataset_datasample_tables('2017/03','2017/03','2017-03')
# 	t3 = time.time()
# 	print 'TIME(Mar):'
# 	print t3-t2
# 	generate_dataset_datasample_tables('2017/04','2017/04','2017-04')
# 	t4 = time.time()
# 	print 'TIME(Apr):'
# 	print t4-t3
# 	generate_dataset_datasample_tables('2017/05','2017/05','2017-05')
# 	t5 = time.time()
# 	print 'TIME(May):'
# 	print t5-t4

# TIME(Jan):	# (25)
# 568.598527908
# TIME(Feb):	# (15)
# 378.935008049
# TIME(Mar):	# (15)
# 306.390367031
# TIME(Apr):	# (20)
# 732.77582407
# TIME(May):	# (26)
# 574.97438097

# 	t0 = time.time()
# 	generate_dataset_datasample_tables('2016/01','2016/01','2016-01')
# 	t1 = time.time()
# 	print 'TIME(Jan-2016):'
# 	print t1-t0
# 	generate_dataset_datasample_tables('2016/02','2016/02','2016-02')
# 	t2 = time.time()
# 	print 'TIME(Feb-2016):'
# 	print t2-t1
# 	generate_dataset_datasample_tables('2016/03','2016/03','2016-03')
# 	t3 = time.time()
# 	print 'TIME(Mar-2016):'
# 	print t3-t2
# 	generate_dataset_datasample_tables('2016/04','2016/04','2016-04')
# 	t4 = time.time()
# 	print 'TIME(Apr-2016):'
# 	print t4-t3
# 	generate_dataset_datasample_tables('2016/05','2016/05','2016-05')
# 	t5 = time.time()
# 	print 'TIME(May-2016):'
# 	print t5-t4
# 	generate_dataset_datasample_tables('2016/06','2016/06','2016-06')
# 	t6 = time.time()
# 	print 'TIME(Jun-2016):'
# 	print t6-t5


# TIME(Jan-2016):	# (16)
# 386.633100986
# TIME(Feb-2016):	# (15)
# 406.331851006
# TIME(Mar-2016):	# (23)
# 674.77325511
# TIME(Apr-2016):	# (16)
# 383.377298832
# TIME(May-2016):	# (18)
# 572.716678143
# TIME(Jun-2016):	# (20)
# 607.956492901


# 	t0 = time.time()
# 	generate_dataset_datasample_tables('2016/07','2016/07','2016-07')
# 	t1 = time.time()
# 	print 'TIME(Jul-2016):'
# 	print t1-t0
# 	generate_dataset_datasample_tables('2016/08','2016/08','2016-08')
# 	t2 = time.time()
# 	print 'TIME(Aug-2016):'
# 	print t2-t1
# 	generate_dataset_datasample_tables('2016/09','2016/09','2016-09')
# 	t3 = time.time()
# 	print 'TIME(Sep-2016):'
# 	print t3-t2
# 	generate_dataset_datasample_tables('2016/10','2016/10','2016-10')
# 	t4 = time.time()
# 	print 'TIME(Oct-2016):'
# 	print t4-t3
# 	generate_dataset_datasample_tables('2016/11','2016/11','2016-11')
# 	t5 = time.time()
# 	print 'TIME(Nov-2016):'
# 	print t5-t4
# 	generate_dataset_datasample_tables('2016/12','2016/12','2016-12')
# 	t6 = time.time()
# 	print 'TIME(Dec-2016):'
# 	print t6-t5

# TIME(Jul-2016):	# (4)
# 64.4038999081
# TIME(Aug-2016):	# (19)
# 976.838187933
# TIME(Sep-2016):	# (23)
# 945.410880089
# TIME(Oct-2016):	# (17)
# 584.166773081
# TIME(Nov-2016):	# (24)
# 708.761294842
# TIME(Dec-2016):	# (20)
# 780.703518152

# 	t0 = time.time()
# 	generate_dataset_datasample_tables('2016/01','2016/02','2016-01to02')
# 	t1 = time.time()
# 	print 'TIME(JanFeb-2016):'
# 	print t1-t0
# 	generate_dataset_datasample_tables('2016/03','2016/04','2016-03to04')
# 	t2 = time.time()
# 	print 'TIME(MarApr-2016):'
# 	print t2-t1
# 	generate_dataset_datasample_tables('2016/05','2016/06','2016-05to06')
# 	t3 = time.time()
# 	print 'TIME(MayJun-2016):'
# 	print t3-t2
# 	generate_dataset_datasample_tables('2016/07','2016/08','2016-07to08')
# 	t4 = time.time()
# 	print 'TIME(JulAug-2016):'
# 	print t4-t3
# 	generate_dataset_datasample_tables('2016/09','2016/10','2016-09to10')
# 	t5 = time.time()
# 	print 'TIME(SepOct-2016):'
# 	print t5-t4
# 	generate_dataset_datasample_tables('2016/11','2016/12','2016-11to12')
# 	t6 = time.time()
# 	print 'TIME(NovDec-2016):'
# 	print t6-t5


# TIME(JanFeb-2016):	# (31)
# 797.660248995
# TIME(MarApr-2016):	# (39)
# 1062.00900316
# TIME(MayJun-2016):	# (38)
# 1188.43708491
# TIME(JulAug-2016):	# (23)
# 935.717015982
# TIME(SepOct-2016):	# (40)
# 1490.66893101
# TIME(NovDec-2016):	# (44)
# 1420.69402695


# 	t0 = time.time()
# 	generate_dataset_datasample_tables('2016/01','2016/03','2016-01to03')
# 	t1 = time.time()
# 	print 'TIME(JanMar-2016):'
# 	print t1-t0
# 	generate_dataset_datasample_tables('2016/04','2016/06','2016-04to06')
# 	t2 = time.time()
# 	print 'TIME(AprJun-2016):'
# 	print t2-t1
# 	generate_dataset_datasample_tables('2016/07','2016/09','2016-07to09')
# 	t3 = time.time()
# 	print 'TIME(JulSep-2016):'
# 	print t3-t2
# 	generate_dataset_datasample_tables('2016/10','2016/12','2016-10to12')
# 	t4 = time.time()
# 	print 'TIME(OctDec-2016):'
# 	print t4-t3

# TIME(JanMar-2016):	# (54)
# 1397.26973104
# (54)


# 	generate_dataset_datasample_tables('2016/07','2016/12','2016-07to12')
# 	generate_dataset_datasample_tables('2016/01','2016/06','2016-01to06')
# 	generate_dataset_datasample_tables('2015/07','2015/12','2015-07to12')
# 	generate_dataset_datasample_tables('2015/01','2015/06','2015-01to06')

# 	200076058--series w/ a superseries...
# 	build_entries('200076058')
# 	test_range(79)
# 	test_max()


if __name__ == "__main__":
    main(sys.argv[1:])
# main()

    #print 'done'
