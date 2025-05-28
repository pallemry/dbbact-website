import tempfile
import urllib.parse
from collections import defaultdict
from io import TextIOWrapper
import os
import json
import requests
import re
from io import BytesIO
import base64
import hashlib

import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import scipy.stats

from flask import Blueprint, request, render_template, make_response, redirect, url_for, Markup, render_template_string, send_from_directory, current_app, session, send_file, current_app

from .utils import debug, get_fasta_seqs, get_dbbact_server_address, get_dbbact_server_color
from . import enrichment
import calour as ca
import dbbact_calour.dbbact

Site_Main_Flask_Obj = Blueprint('Site_Main_Flask_Obj', __name__, template_folder='templates')


# the dbbact rest-api server address
dbbact_server_address = get_dbbact_server_address()


# robots.txt file for google/other web crawlers
@Site_Main_Flask_Obj.route("/robots.txt")
def robots_txt():
    return render_template('robots.txt')


# Sitemap for google/other web crawlers
@Site_Main_Flask_Obj.route("/sitemap")
@Site_Main_Flask_Obj.route("/sitemap/")
@Site_Main_Flask_Obj.route("/sitemap.xml")
def sitemap():
    """
        Route to dynamically generate a sitemap of your website/application.
        lastmod and priority tags omitted on static pages.
        lastmod included on dynamic content such as blog posts.

        based on:
        https://gist.github.com/Julian-Nash/aa3041b47183176ca9ff81c8382b655a
    """
    from urllib.parse import urlparse

    host_components = urlparse(request.host_url)
    host_base = host_components.scheme + "://" + host_components.netloc

    # Static routes with static content
    static_urls = list()
    # skip_urls = ['/sequence_annotations', '/ontology_info','/user_info','/annotation_info','/annotation_seqs','/species_info']
    for rule in current_app.url_map.iter_rules():
    #     skip_addr = False
    #     for url in skip_urls:
    #         if str(rule).startswith(url):
    #             skip_addr = True
    #             break
    #     if skip_addr:
    #         continue
        if "GET" in rule.methods and len(rule.arguments) == 0:
            url = {
                "loc": f"{host_base}{str(rule)}"
            }
            static_urls.append(url)

    # Dynamic routes with dynamic content
    dynamic_urls = list()

    xml_sitemap = render_template("sitemap.xml", static_urls=static_urls, dynamic_urls=dynamic_urls, host_base=host_base)
    response = make_response(xml_sitemap)
    response.headers["Content-Type"] = "application/xml"

    return response


@Site_Main_Flask_Obj.route('/', methods=['POST', 'GET'])
def landing_page():
    '''
    Redirect to the main search page
    '''
    # TODO: fix to non hard-coded
    return redirect('main')


@Site_Main_Flask_Obj.route('/main', methods=['POST', 'GET'])
def main_html():
    debug(1, 'main page')
    """Title: the main dbBact page and search tool
    URL: dbbact.org/main
    Method: GET
    """
    # get the dbbact statistics from the dbbact rest-api server
    debug(2, 'getting stats from %s' % dbbact_server_address)
    NumAnnotation = 0
    NumSequences = 0
    NumSequenceAnnotation = 0
    NumExperiments = 0
    dbbact_api_server_type = 'unknown'
    try:
        httpRes = requests.get(dbbact_server_address + '/stats/stats')
        if httpRes.status_code == 200:
            jsonRes = httpRes.json()
            # NumOntologyTerms = jsonRes.get("stats").get('NumOntologyTerms')
            NumAnnotation = jsonRes.get("stats").get('NumAnnotations')
            NumSequences = jsonRes.get("stats").get('NumSequences')
            NumSequenceAnnotation = jsonRes.get("stats").get('NumSeqAnnotations')
            NumExperiments = jsonRes.get("stats").get('NumExperiments')
            dbbact_api_server_type = jsonRes.get("stats").get("Database")
    except:
        pass
    # NumOntologyTerms = 0

    webPage = render_template('searchpage.html',
                              alert_text=get_alert_text(),
                              header_color=get_dbbact_server_color(api_server_type=dbbact_api_server_type),
                              numAnnot=(str(NumAnnotation).replace('.0', '')),
                              numSeq=(str(NumSequences).replace('.0', '')),
                              numExp=(str(NumExperiments).replace('.0', '')),
                              numSeqAnnot=(str(NumSequenceAnnotation).replace('.0', '')))
    return webPage


@Site_Main_Flask_Obj.route('/enrichment', methods=['POST', 'GET'])
def test_enrichment():
    '''
    Redirect to enrichment page
    '''
    # TODO: fix to non hard-coded
    httpRes = requests.get(dbbact_server_address + '/stats/stats')
    # NumOntologyTerms = 0
    NumAnnotation = 0
    NumSequences = 0
    NumSequenceAnnotation = 0
    NumExperiments = 0
    if httpRes.status_code == 200:
        jsonRes = httpRes.json()
        # NumOntologyTerms = jsonRes.get("stats").get('NumOntologyTerms')
        NumAnnotation = jsonRes.get("stats").get('NumAnnotations')
        NumSequences = jsonRes.get("stats").get('NumSequences')
        NumSequenceAnnotation = jsonRes.get("stats").get('NumSeqAnnotations')
        NumExperiments = jsonRes.get("stats").get('NumExperiments')

    webPage = render_template('enrichment.html',
                              numAnnot=(str(NumAnnotation).replace('.0', '')),
                              numSeq=(str(NumSequences).replace('.0', '')),
                              numExp=(str(NumExperiments).replace('.0', '')),
                              numSeqAnnot=(str(NumSequenceAnnotation).replace('.0', '')))
    return webPage


def build_res_html(success, expId, isNewExp, annotId, additional=None):
    successStr = ''
    expIdStr = ''
    existingStr = ''
    annotIdStr = ''
    debugStr = ''

    if success is True:
        successStr = 'Operation succeed'
    else:
        successStr = 'Operation failed'

    if expId == -1:
        expIdStr = 'NA'
        existingStr = ''
    else:
        # expIdStr = "<a href='http://127.0.0.1:5000/exp_info/" + str(expId) + "'>" + str(expId) + "</a>"
        expIdStr = "<a href='%s'>" % url_for('.exp_info', expid=expId) + + str(expId) + "</a>"
        if isNewExp is True:
            existingStr = '(new)'
        else:
            existingStr = '(existing)'

    if annotId == -1:
        annotIdStr = 'NA'
    else:
        annotIdStr = str(annotId)

    if additional is not None:
        debugStr = "Error information : " + additional

    webStr = render_header(title='Error') + render_template('add_data_results.html', title_str=successStr, new_or_existing_str=existingStr, annotation_id=annotIdStr, exp_id=expIdStr, debug_info=debugStr)
    return webStr


@Site_Main_Flask_Obj.route('/add_data_results', methods=['POST', 'GET'])
def add_data_results():
    """
    Title: Add data processing
    URL: site/add_data_results
    Method: POST
    """
    webPageTemp = ''
    webpage = ''

    if 'fastaFileTb' in request.files:
        debug(1, 'Fasta file uploaded, processing it')
        try:
            file1 = request.files['fastaFileTb']
            textfile1 = TextIOWrapper(file1)
            seqs1 = get_fasta_seqs(textfile1)
        except:
            webpage = build_res_html(False, -1, isNewExp, annotId,'Could not open fasta file')
            return(webpage, 400)
        if seqs1 is None:
            webpage = build_res_html(False, -1, False, -1, 'Invalid fasta file')
            return(webpage, 400)
    else:
        webpage = build_res_html(False, -1, False, -1, 'Missing fasta file')
        return(webpage, 400)

    # Prepare all exp data in array
    methodName = request.form.get('methodNameTb')
    if methodName is None or len(methodName.strip()) == 0:
        methodName = 'website'

    descName = request.form.get('descNameTb')
    if descName is None or len(descName.strip()) == 0:
        descName = 'na'
    # print(">>>>>>>>>>>>>><<<<<<<<<<<<<<<<<method name" + methodName)

    # Exp list
    hiddenExpName = request.form.get('hiddenExpName')
    hiddenExpValue = request.form.get('hiddenExpValue')
    # Ont list
    hiddenOntName = request.form.get('hiddenOntName')
    hiddenOntType = request.form.get('hiddenOntType')
    hiddenOntDetType = request.form.get('hiddenOntDetType')

    hiddenRegionStr = request.form.get('hiddenRegion')
    if hiddenRegionStr is None or len(hiddenRegionStr.strip()) == 0:
        webpage = build_res_html(False, -1, False, -1, 'Invalid region value')
        return(webpage, 400)

    if hiddenOntType is None or len(hiddenOntType.strip()) == 0:
        webpage = build_res_html(False, -1, False, -1, 'Invalid Input (error code: -1)')
        return(webpage, 400)

    # in case one of the parameters is missing
    if hiddenExpName is None or len(hiddenExpName.strip()) == 0 or hiddenExpValue is None or len(hiddenExpValue.strip()) == 0 or hiddenOntName is None or len(hiddenOntName.strip()) == 0 or hiddenOntDetType is None or len(hiddenOntDetType.strip()) == 0:
        webpage = build_res_html(False, -1, False, -1, 'Invalid Input (error code: -2)')
        return(webpage, 400)

    expDataNameArr = hiddenExpName.split(';')  # split string into a list
    expDataValueArr = hiddenExpValue.split(';')  # split string into a list

    ontDataNameArr = hiddenOntName.split(';')  # split string into a list
    annotationTypeArr = hiddenOntType.split(';')  # split string into a list
    annotationDetTypeArr = hiddenOntDetType.split(';')  # split string into a list

    # ### Strings to return
    # {{title_str}} - Operation failed / Operation completed successfully
    # {{exp_str}}
    # {{anot_str}}
    resTitleStr = ''
    resExpStr = ''
    resAnotStr = ''
    newExpFlag = False

    '''
        #####################################################
        # Get expirement id or -1 if doesn't exist
        #####################################################
    '''
    rdata = {}
    rdata['nameStrArr'] = expDataNameArr
    rdata['valueStrArr'] = expDataValueArr
    httpRes = requests.get(dbbact_server_address + '/experiments/get_id_by_list', json=rdata)
    if httpRes.status_code == 200:
        jsonRes = httpRes.json()
        expId = jsonRes.get("expId")
        errorCode = jsonRes.get("errorCode")
        errorText = jsonRes.get("errorText")
        if expId < 0:
            # identification appears in more than one expirement
            if errorCode == -2:
                webpage = build_res_html(False, -1, False, -1, 'More than one experiments was found')
                return(webpage, 400)
            # expirement was not found, try to find
            elif errorCode == -1:
                rdataExp = {}
                test = []

                for i in range(len(expDataNameArr)):
                    test.append((expDataNameArr[i], expDataValueArr[i]))

                rdataExp = {'expId': -1, 'private': False, 'details': test}

                httpRes = requests.post(dbbact_server_address + '/experiments/add_details', json=rdataExp)
                if httpRes.status_code == 200:
                    jsonRes = httpRes.json()
                    expId = jsonRes.get("expId")
                    newExpFlag = True
                else:
                    webpage = build_res_html(False, -1, False, -1, 'Failed to create new expirement')
                    return(webpage, 400)
    else:
        webpage = build_res_html(False, -1, False, -1, 'Failed to get expirement id')
        return(webpage, 400)

    # ####################################################
    # Add sequences if they are missing
    # ####################################################
    rdata = {}
    rdata['sequences'] = seqs1
    rdata['primer'] = hiddenRegionStr

    httpRes = requests.post(dbbact_server_address + '/sequences/add', json=rdata)
    if httpRes.status_code == 200:
        jsonRes = httpRes.json()
        seqList = jsonRes.get("seqIds")
        if len(seqList) != len(seqs1):
            webpage = build_res_html(False, expId, newExpFlag, -1, 'Failed to retrieve all sequneces IDs')
            return(webpage, 400)
    else:
        webpage = build_res_html(False, expId, newExpFlag, -1, 'Failed to retrieve sequneces IDs')
        return(webpage, 400)
    #####################################################      
    
    #####################################################      
    # Get ontologies list
    #####################################################      
    rdata = {}
    rdata['ontologies'] = ontDataNameArr

    httpRes = requests.post(dbbact_server_address + '/ontology/get',json=rdata)
    if httpRes.status_code == 200:
        jsonRes = httpRes.json()
        ontList = jsonRes.get("ontIds")
        if len(ontDataNameArr) != len(ontList) :
            webpage = build_res_html(False, expId, newExpFlag, -1 , 'Failed to retrieve ontologies IDs')
            return(webpage, 400)
    else:
        webpage = build_res_html(False, expId, newExpFlag, -1 , 'Failed to retrieve ontologies IDs')
        return(webpage, 400)
    # ####################################################
    # for i in range(len(ontList)):
    #    webpage += str(ontList[i]) + "<br>"

    annotationListArr = []

    for i in range(len(ontDataNameArr)):
        annotationListArr.append((annotationDetTypeArr[i],ontDataNameArr[i]))

    rannotation = {}
    rannotation['expId'] = expId
    rannotation['sequences'] = seqs1
    rannotation['region'] = hiddenRegionStr
    rannotation['annotationType'] = annotationTypeArr
    rannotation['method'] = methodName
    rannotation['agentType'] = 'DBBact website submission'
    rannotation['description'] = descName
    rannotation['annotationList'] = annotationListArr
    # Everything is ready to add the data
    httpRes = requests.post(dbbact_server_address + '/annotations/add', json=rannotation)
    if httpRes.status_code == 200:
        jsonRes = httpRes.json()
        annotId = jsonRes.get("annotationId")
        webpage = build_res_html(True, expId, newExpFlag, annotId)
        return(webpage, 400)
    else:
        webpage = build_res_html(False, expId, newExpFlag, -1, 'Failed to add annotations')
        return(webpage, 400)

    return ''


@Site_Main_Flask_Obj.route('/enrichment_results', methods=['POST', 'GET'])
def enrichment_results():
    """
    Title: Search results page
    URL: site/search_results
    Method: POST
    """
    webPageTemp = ''

    if request.method == 'GET':
        exampleStr = request.args['example']
    else:
        exampleStr = request.form['example']

    # first file
    if exampleStr != 'true':
        if 'seqs1' in request.files:
            debug(1, 'Fasta file uploaded, processing it')
            file1 = request.files['seqs1']
            print(file1)
            # textfile1 = file1.read()
            # print(textfile1)
            if hasattr(file1, '_file'):
                textfile1 = TextIOWrapper(file1._file)
            else:
                textfile1 = TextIOWrapper(file1)
            seqs1 = get_fasta_seqs(textfile1)
            if seqs1 is None:
                webPageTemp = render_header(title='Error') + render_template('error_page.html', error_str='Error: Uploaded file1 not recognized as fasta')
                return(webPageTemp, 400)
        else:
            webPageTemp = render_header(title='Error') + render_template('error_page.html', error_str='Error: Missing fasta file name 1')
            return(webPageTemp, 400)
    else:
        # only used for example query
        with open("dbbact_website/enrichment_example/seqs-fec.fa", "r") as myfile:
            textfile1 = myfile.readlines()
            seqs1 = get_fasta_seqs(textfile1)
            if seqs1 is None:
                webPageTemp = render_header(title='Error') + render_template('error_page.html', error_str='Error: Uploaded file1 not recognized as fasta')
                return(webPageTemp, 400)
    debug(2, 'Loaded %d sequences for seqs1 file' % len(seqs1))

    # second file
    if exampleStr != 'true':
        if 'seqs2' in request.files:
            debug(1, 'Fasta file uploaded, processing it')
            file2 = request.files['seqs2']
            if hasattr(file2, '_file'):
                textfile2 = TextIOWrapper(file2._file)
            else:
                textfile2 = TextIOWrapper(file2)
            seqs2 = get_fasta_seqs(textfile2)
            if seqs2 is None:
                webPageTemp = render_header(title='Error') + render_template('error_page.html', error_str='Error: Uploaded file2 not recognized as fasta')
                return(webPageTemp, 400)
        else:
            webPageTemp = render_header(title='Error') + render_template('error_page.html', error_str='Error: Missing fasta file name')
    else:
        # only used for example query
        with open("dbbact_website/enrichment_example/seqs-sal.fa", "r") as myfile:
            textfile2 = myfile.readlines()
            seqs2 = get_fasta_seqs(textfile2)
            if seqs2 is None:
                webPageTemp = render_header(title='Error') + render_template('error_page.html', error_str='Error: Uploaded file1 not recognized as fasta')
                return(webPageTemp, 400)
    debug(2, 'Loaded %d sequences for seqs2 file' % len(seqs2))

    webpage = render_header()
    # webpage = render_template('info_header.html')
    for term_type in ['term', 'annotation']:
        webpage += "<h2>%s enrichment</h2>" % term_type
        webpage += '(negative (red) LOWER in fasta file 1, positive (blue) HIGHER in fasta file 1)<br>'
        webpage += render_template('enrichment_results.html')

        debug(2, 'looking for enriched %s' % term_type)
        err, terms, pval, odif = enrichment.enrichment(seqs1, seqs2, term_type=term_type)
        # err, terms, pval, odif = enrichment.calour_enrichment(seqs1, seqs2, term_type=term_type)
        debug(2, 'found %d enriched' % len(terms))

        if err:
            webPageTemp = render_header(title='Error') + render_template('error_page.html', error_str=err)
            return(webPageTemp, 400)
        for idx, cterm in enumerate(terms):
            if odif[idx] < 0:
                ccolor = 'red'
            else:
                ccolor = 'blue'
            webpage += '<tr><td><span style="color:%s">%s</span></td>' % (ccolor, cterm)
            webpage += '<td>%f</td>' % odif[idx]
            webpage += '<td>%f</td>' % pval[idx]
            webpage += "</tr>"
            webpage += '</span>'
        webpage += "</table>"
    return webpage


@Site_Main_Flask_Obj.route('/search_results', methods=['POST', 'GET'])
def search_results():
    """
    Title: Search results page
    URL: dbbact.org/search_results
    Method: GET / POST
    """
    webPageTemp = ''
    if request.method == 'GET':
        sequence = request.args['sequence']
    else:
        sequence = request.form['sequence']

    # if we have a fasta file attached, process it
    if sequence == '':
        if 'fasta file' in request.files:
            debug(2, 'Fasta file uploaded, processing it')
            # save the file locally, get the sequences and remove it
            file = request.files['fasta file']
            filepos = os.path.join(tempfile.gettempdir(), tempfile.gettempprefix())
            debug(2, 'fasta tmp file is %s' % filepos)
            file.save(filepos)
            with open(filepos) as textfile:
                seqs = get_fasta_seqs(textfile)
                debug(2,'sequences: %s' % seqs)
            os.remove(filepos)
            # if no sequences in fasta file - return error message
            if seqs is None:
                webPageTemp = render_header(title='Error') + render_template('error_page.html', error_str='Error: Uploaded file not recognized as fasta')
                return(webPageTemp, 400)
            # return the webpage for the group annotations
            err, webpage = draw_sequences_annotations_compact(seqs)
            return webpage

    # if it is short, try if it is a greengenesID/ontology term/taxonomy
    if len(sequence) < 50:
        # # if number assume it is a greengenes id
        # if sequence.isdigit():
        #     # try gg first
        #     err, webPage = get_gg_info(sequence)
        #     if not err:
        #         debug(2, 'get info for greengenes id %s' % sequence)
        #         return webPage

        #     debug(2, 'get info for greengenesid %s' % sequence)
        #     webPage = sequence_annotations(sequence)
        #     return webPage

        # try is it an ontology term
        err, webPage = get_ontology_info(sequence)
        if not err:
            debug(2, 'get info for ontology term %s' % sequence)
            return webPage
        # or maybe a taxonomy term
        err, webPage = get_taxonomy_info(sequence)
        if not err:
            debug(2, 'get info for taxonomy %s' % sequence)
            return webPage
        # maybe a SILVA species name?
        if len(sequence.split(' ')) > 1:
            err, webPage = get_species_info(sequence)
            if not err:
                debug(2, 'get species info')
                return webPage
        # or maybe based on qiime2 hash string
        err, webPage = get_hash_info(sequence)
        if not err:
            debug(2, 'get info for qiime2 hash %s' % sequence)
            return webPage
        # or maybe based on silva id
        err, webPage = get_silva_info(sequence)
        if not err:
            debug(2, 'get info for silva id %s' % sequence)
            return webPage
        # so we can't find it
        debug(2, 'search sequence/term/etc not found for %s' % sequence)
        return error_message('Not found', 'Keyword "%s" was not found in dbBact ontology '
                             'or taxonomy.' % sequence)

    # 50 < length < 100 - it's a sequence but not long enough
    if len(sequence) < 100:
        debug(2, 'sequence too short len=%d' % len(sequence))
        webPageTemp = render_header(title='Error') + render_template('error_page.html', error_str='Sequences must be at least 100bp long.')
        return(webPageTemp, 400)

    # so it's a legit sequence - let's get the annotations for it
    debug(2, 'guess it is a sequence')
    webPage = sequence_annotations(sequence)
    return webPage


@Site_Main_Flask_Obj.route('/sequences_wordcloud', methods=['POST'])
def sequences_wordcloud():
    '''show the info page (wordcloud, term table, etc.) for a set of sequences

    Parameters
    ----------
    sequences: list of str
        the sequences to get the annotations for
    ignore_exp: list of int, optional
        dbbact experiment ids to exclude from the analysis
    '''
    alldat = request.get_json()
    sequences = alldat.get('sequences')
    ignore_exp = alldat.get('ignore_exp')

    err, webpage = draw_sequences_annotations_compact(sequences, ignore_exp=ignore_exp)
    return webpage


@Site_Main_Flask_Obj.route('/sequence_annotations/<string:sequence>')
def sequence_annotations(sequence):
    '''Get annotations for a given sequence and return the info page about it
    (including species info, wordcloud, term table, etc.)

    Parameters
    ----------
    sequence: str
        the sequence to look for (ACGT string)
    '''

    # first we check whether the sequence is in dbbact - if not, we'll try to trim it
    trim_msg = ''
    found_seq = test_if_sequence_exists(sequence)

    # indicates if trimming was needed (in order to show the trim message)
    show_trim_msg = False

    # indicates if only a <100% identity sequence was found (in order to show the mismatch message)
    found_only_mismatch = False

    # if we didn't find any match, try to trim the sequence from known primers
    if not found_seq:
        trimmed, trim_msg = trim_primers_from_sequence(sequence)
        if trimmed is not None:
            show_trim_msg = True
            sequence = trimmed
            found_seq = True
    if not found_seq:
        close_seqs = get_close_sequences(sequence)
        if len(close_seqs) > 0:
            found_only_mismatch = True
            found_seq = True
            if len(close_seqs) > 1:
                num_mismatches = [cseq['num_mismatches'] for cseq in close_seqs]
                seqs = [cseq['sequence'] for cseq in close_seqs]
                err, webpage = draw_sequences_annotations_compact(seqs, inexact_match=True, num_mismatches=num_mismatches)
                return webpage
            else:
                sequence = close_seqs[0]['sequence']
                num_mismatches = close_seqs[0]['num_mismatches']
        else:
            # we didn't find any match, so lets try if the sequence is left trimmed
            pass

    # Get the taxonomy for the sequence
    rdata = {}
    rdata['sequence'] = sequence
    taxStr = "na"
    httpResTax = requests.get(dbbact_server_address + '/sequences/get_taxonomy_str', json=rdata)
    if httpResTax.status_code == requests.codes.ok:
        taxStr = httpResTax.json().get('taxonomy')

    # Get the species and taxonomies based on 100% matching to the whole sequence database (i.e. silva)
    species = []
    num_species_match = 0
    httpResTax = requests.get(dbbact_server_address + '/sequences/get_whole_seq_taxonomy', json=rdata)
    if httpResTax.status_code == requests.codes.ok:
        species = httpResTax.json().get('species')
        ids = httpResTax.json().get('ids')
    else:
        return httpResTax
    species_details = ''
    species_dict = defaultdict(int)
    species_ids = defaultdict(list)
    for cspecies, cwsid in zip(species, ids):
        if cspecies == '':
            continue
        species_dict[cspecies] += 1
        species_ids[cspecies].append(cwsid)
        num_species_match += 1
    species_dict = dict(sorted(species_dict.items(), key=lambda item: item[1], reverse=True))
    for ck, cv in species_dict.items():
        species_details += '<tr><td><a href=%s>%s</a></td><td>%s</td><td>' % (url_for('.species_info', species=ck), ck, cv)
        for idx, ccid in enumerate(species_ids[ck]):
            if idx > 5:
                break
            species_details += '<a href=%s target="_blank">%s</a>, ' % ('https://www.arb-silva.de/browser/ssu-138.1/' + ccid, ccid)
        if idx > 5:
            species_details += ',...'
        species_details += '</td></tr>'

    # Create the results page
    # the sequence info part, with species details
    webPage = render_header(title='dbBact sequence annotation')

    if found_only_mismatch:
        webPage += '<h2>Exact sequence was not found in dbBact</h2>'
        webPage += '<p>However, a single sequence with a high similarity to your sequence was found in dbBact. '
        webPage += 'Showing results for this sequence</p>'

    webPage += render_template('seqinfo.html', sequence=sequence.upper(), taxonomy=taxStr, species_details=species_details, num_species_match=num_species_match)

    if show_trim_msg:
        if found_seq:
            webPage += '<h2>Sequence found after trimming primers: %s</h2>' % trim_msg
            webPage += '<a href=%s target="_blank">Click here</a> for details how to prepare your 16S sequence for querying dbBact<br>' % (url_for('.faq_prepare_sequences'))
        else:
            pass

    # Get the annotations for the sequence
    httpRes = requests.get(dbbact_server_address + '/sequences/get_annotations', json=rdata)
    if httpRes.status_code != requests.codes.ok:
        # problem with the annotations per sequence
        debug(6, "sequence annotations Error code:" + str(httpRes.status_code))
        webPage = render_header(title='dbBact sequence annotation')
        webPage += "Failed to get annotations for sequence:\n%s" % Markup.escape(sequence)
        webPage += render_template('footer.html')
    else:
        # draw the annotations info
        annotations = httpRes.json().get('annotations')
        webPage += render_sequence_annotations(annotations, sequence=sequence)
        webPage += render_template('footer.html')
    return webPage


def render_sequence_annotations(annotations, sequence=None):
    '''Draw the webpage summarizing the annotations and taxonomy of a given sequence

    Parameters
    ----------
    annotations: 
    sequence: str or None
        the sequence the annotations are for
    '''
    webPage = ''
    if len(annotations) == 0:
        webPage += '<br><br><h1>No annotations for sequence found in dbBact</h1>'
        webPage += '<h2>Are you using >100bp sequences?</h2>'
        webPage += '<h2><br><a href=%s target="_blank">Click here</a> for details how to prepare your 16S sequence for querying dbBact</h2>' % (url_for('.faq_prepare_sequences'))
    else:
        for cannotation in annotations:
            cannotation['website_sequences'] = [0]
        annotations = sorted(annotations, key=lambda x: x.get('num_sequences', 0), reverse=False)
        term_info = get_term_info_for_annotations(annotations)
        webPage += draw_annotation_details(annotations, term_info=term_info, show_relative_freqs=True, sequences=[sequence])

    return webPage


def get_annotations_terms(annotations, get_low=True):
    '''
    Get a list of terms present in the annotations

    Parameters
    ----------
    annotations : list of annotations

    Returns
    -------
    terms : list of str
    list of terms from all annotations (each term appears once)
    '''
    terms = set()
    for cannotation in annotations:
        details = cannotation['details']
        for cdet in details:
            ctype = cdet[0]
            cterm = cdet[1]
            if ctype == 'low':
                cterm = '-' + cterm
            terms.add(cterm)
    terms = list(terms)
    return terms


def draw_sequences_annotations(seqs):
    '''Draw the webpage for annotations for a list of sequences

    Parameters
    ----------
    seqs : list of str
        list of DNA sequences sequences to get annotations for

    Returns
    -------
    err : str
        the error encountered or '' if ok
    webpage : str
        the webpage for the annotations of these sequences
    '''
    res = requests.get(get_dbbact_server_address() + '/sequences/get_list_annotations',
                       json={'sequences': seqs})
    if res.status_code != 200:
        msg = 'error getting annotations for sequences : %s' % Markup.escape(res.content)
        debug(6, msg)
        return msg, msg
    seqannotations = res.json()['seqannotations']
    if len(seqannotations) == 0:
        msg = 'no sequences found'
        return msg, msg
    annotations = []
    for cseqannotation in seqannotations:
        if len(cseqannotation) == 0:
            continue
        for cannotation in cseqannotation:
            annotations.append(cannotation)

    webPage = render_header()
    webPage += '<h2>Annotations for sequence list:</h2>'
    webPage += draw_annotation_details(annotations, sequences=seqs)
    webPage += render_template('footer.html')
    return '', webPage


def draw_sequences_annotations_compact(seqs, ignore_exp=[], draw_only_details=False, inexact_match=None, num_mismatches=None):
    '''Draw the webpage for annotations for a set of sequences

    Parameters
    ----------
    seqs : list of str sequences (ACGT)
    ignore_exp : list of int (optional)
        list of experiment ids to ignore when calculating the score. None to include all experiments
    draw_only_details: bool, optional
        True to plot only the annotations part (no header/footer)
        False to draw complete page
    inexact_match: bool, optional
        If True, results are for inexact matches to a given original sequence
    num_mismatches: list of int or None, optional
        if not None and inexcat_match is True, the number of mismatches for each sequence

    Returns
    -------
    err : str
        the error encountered or '' if ok
    webpage : str
        the webpage for the annotations of these sequences
    '''
    # get the compact annotations for the sequences
    debug(2, 'draw_sequences_annotations_compact for %d sequences. ignore_exp=%s' % (len(seqs), ignore_exp))
    res = requests.get(get_dbbact_server_address() + '/sequences/get_fast_annotations',
                       json={'sequences': seqs})
    if res.status_code != 200:
        msg = 'error getting annotations for sequences : %s' % Markup.escape(res.content)
        debug(6, msg)
        return msg, msg

    res = res.json()
    annotations = res['annotations']
    seqannotations = res['seqannotations']
    if len(seqannotations) == 0:
        msg = 'None of the %d sequences were found in dbBact. Are these >100bp long 16S sequences?\nNote dbBact is populated mostly by EMP V4 (515F) amplicon sequences.' % len(seqs)
        debug(3, msg)
        return msg, msg
    term_info = res['term_info']

    if draw_only_details:
        webPage = ''
    else:
        webPage = render_header()
    if not inexact_match:
        webPage += '<h2>Annotations for %d sequences</h2>' % len(seqs)
    else:
        webPage += '<h2>No exact match found</h2>'
        seqs_details = ''
        for i, cseq in enumerate(seqs):
            if num_mismatches is not None:
                cnum_mismatch = num_mismatches[i]
            else:
                cnum_mismatch = 'NA'
            seqs_details += ('<tr><td><a href=%s target="_blank">%s</a></td><td>%s</td></tr>' % (url_for('.sequence_annotations', sequence=cseq),  Markup.escape(cseq.upper()), cnum_mismatch))
        webPage += render_template('seq-list-collapse.html', title='Inexact match sequences', seq_count=len(seqs), details=seqs_details)
    webPage += draw_group_annotation_details(annotations, seqannotations, term_info=term_info, ignore_exp=ignore_exp, sequences=seqs)
    if not draw_only_details:
        webPage += render_template('footer.html')
    return '', webPage


def getannotationstrings(cann, use_descriptions=True):
    """
    get a nice string summary of a curation

    input:
    cann : dict from /sequences/get_annotations (one from the list)
    use_descriptions : bool
        True to include the free text description in the output

    output:
    cdesc : str
        a short summary of each annotation
    """
    cdesc = ''
    if use_descriptions:
        if cann['description']:
            cdesc += cann['description'] + ' ('
    if cann['annotationtype'] == 'diffexp':
        chigh = []
        clow = []
        call = []
        for cdet in cann['details']:
            if cdet[0] == 'all':
                call.append(cdet[1])
                continue
            if cdet[0] == 'low':
                clow.append(cdet[1])
                continue
            if cdet[0] == 'high':
                chigh.append(cdet[1])
                continue
        cdesc += ' high in '
        for cval in chigh:
            cdesc += cval + ' '
        cdesc += ' compared to '
        for cval in clow:
            cdesc += cval + ' '
        cdesc += ' in '
        for cval in call:
            cdesc += cval + ' '
    elif cann['annotationtype'] == 'isa':
        cdesc += ' is a '
        for cdet in cann['details']:
            cdesc += 'cdet,'
    elif cann['annotationtype'] == 'contamination':
        cdesc += 'contamination'
    else:
        cdesc += cann['annotationtype'] + ' '
        for cdet in cann['details']:
            cdesc = cdesc + ' ' + cdet[1] + ','

    if len(cdesc) >= 1 and cdesc[-1] == ',':
        cdesc = cdesc[:-1]

    if use_descriptions:
        if cann['description']:
            cdesc += ')'
    return cdesc


@Site_Main_Flask_Obj.route('/annotation_info/<int:annotationid>')
def annotation_info(annotationid):
    """
    get the information about an annotation
    input:
    annotationid : int
        the annotationid to get the info for
    """
    # get the experiment info for the annotation
    rdata = {}
    rdata['annotationid'] = annotationid
    # get the experiment annotations
    res = requests.get(get_dbbact_server_address() + '/annotations/get_annotation', params=rdata)
    if res.status_code != 200:
        message = Markup('Annotation ID <b>%d</b> was not found.' % annotationid)
        webPageTemp = render_header(title='Error') + render_template('error_page.html', error_str='Not found')
        return(webPageTemp, 400)

    annotation = res.json()
    # get the experiment details
    rdata = {}
    expid = annotation['expid']
    rdata['expId'] = expid
    webPage = render_header(title='Annotation %s' % annotationid)
    webPage += render_template('annotinfo.html', annotationid=annotationid)
    res = requests.get(dbbact_server_address + '/experiments/get_details', json=rdata)
    if res.status_code == 200:
        webPage += draw_experiment_info(expid, res.json()['details'])
    else:
        message = Markup('Error getting experiment details.')
        return(render_header(title='Not found') +
               render_template('error.html', title='Not found',
                               message=message) +
               render_template('footer.html'), 400)
    webPage += '<h2>Annotations Details</h2>'
    webPage += draw_annotation_table([annotation])

    webPage += render_template('annotdetail.html')
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('description', Markup.escape(annotation['description']))
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('type', Markup.escape(annotation['annotationtype']))
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('num_sequences', Markup.escape(annotation['num_sequences']))
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('region', Markup.escape(annotation['primer']))
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('flags', len(annotation['flags']))
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('agent', Markup.escape(annotation['agent']))
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('method', Markup.escape(annotation['method']))
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('date', Markup.escape(annotation['date']))
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('username', Markup.escape(annotation['username']))
    webPage += '<tr><td>%s</td><td>%s</td></tr>' % ('private', Markup.escape(annotation['private']))
    annotationdetails = annotation['details']

    for cad in annotationdetails:
        webPage += "<tr>"
        webPage += '<td>' + str(cad[0]) + '</td>'
        webPage += '<td><a href=' + urllib.parse.quote('../ontology_info/' + str(cad[1])) + '>' + str(cad[1]) + '</a></td></tr>'

    webPage += '</table>'

    if len(annotation['flags']) > 0:
        webPage += list_flags(annotation['flags'])
    else:
        webPage += '<br>Annotation not flagged as potentially problematic<br>'

    webPage += draw_flag_annotation_button(annotationid)

    webPage += '<h2>Review status</h2>'
    review_status = annotation.get('review_status', -1)
    if review_status == 0:
        webPage += '<br>Pending review (annotation has not been reviewed yet by the dbBact team)<br>'
    elif review_status == 1:
        webPage += '<br>Accepted (annotation has been reviewed by the dbBact team and accepted)<br>'
    else:
        webPage += '<br>NA<br>'

    webPage += '<h2>Sequences</h2>'
    webPage += draw_download_fasta_button(annotationid)

    # add the ontology parent terms for the annotation
    webPage += '<h2>Ontology terms (including parents)</h2>'
    res = requests.get(get_dbbact_server_address() + '/annotations/get_annotation_ontology_parents', json={'annotationid': annotationid})
    if res.status_code != 200:
        debug(6, 'no ontology parents found for annotationid %d' % annotationid)
        parents = []
    else:
        parents = res.json().get('parents')
        debug(1, 'found %d parent groups for annotationid %d' % (len(parents), annotationid))
    webPage += '<div style="margin: 20px"><blockquote style="font-size: 1em;">'
    for ctype, cparents in parents.items():
        cparents = list(set(cparents))
        webPage += '<p>%s: ' % ctype
        for cparentname in cparents:
            webPage += '<a href=' + urllib.parse.quote('../ontology_info/' + str(cparentname)) + '>' + str(Markup.escape(cparentname)) + '</a> '
        webPage += '</p>'
    webPage += '</blockquote></div>'
    webPage += render_template('footer.html')
    return webPage


def list_flags(flags):
    '''return an html table with the flags

    Parameters
    ----------
    flags: list of dict returned to the 'flags' field of an annotation

    Returns
    -------
    html table with list of flags
    '''
    webPage = render_template('flag_info.html')
    for cflag in flags:
        webPage += '<tr>'
        webPage += '<td>%s</td>' % cflag['reason']
        webPage += '<td>%s</td>' % cflag['status']
        webPage += '<td>%s</td>' % cflag['userid']
        webPage += '</tr>'
    webPage += '</table>'
    return webPage


def draw_flag_annotation_button(annotationid):
    '''Draw a button to flag the annotation as suspicious

    Parameters
    ----------
    annotationid : int
        the annotationid to flag

    Returns
    -------
    webPage : str
        html for the download button with the link to the fasta file download page
    '''
    webPage = '<div style="margin: 20px"><button class="btn btn-default" onclick="location.href=\'%s\';"><i class="glyphicon glyphicon glyphicon-thumbs-down"></i> Flag as suspicious</button></div>' % url_for('.annotation_flag', annotationid=annotationid)
    return webPage


@Site_Main_Flask_Obj.route('/annotation_flag/<int:annotationid>')
def annotation_flag(annotationid):
    '''webpage used to request input from the user to flag an annotation'''
    webpage = render_header()
    webpage += render_template('annotation_flag.html', annotationid=annotationid, back_addr=url_for('.annotation_info', annotationid=annotationid))
    webpage += render_template('footer.html', header_color=get_dbbact_server_color())
    return webpage


@Site_Main_Flask_Obj.route('/annotation_flag_submit', methods=['POST', 'GET'])
def annotation_flag_submit():
    """
    result of the annotation_flag form. Will flag the annotation as anonimous user
    """
    reason = request.form['reason']
    annotationid = request.form['annotationid']

    data = {}
    data['annotationid'] = annotationid
    data['reason'] = reason
    httpRes = requests.post(dbbact_server_address + '/annotations/add_annotation_flag', json=data)
    if httpRes.status_code == 200:
        webpage = render_header(title='Password Recovery')
        webpage += 'Annotation %s has been flagged and will be reviewed by the dbbact team<br>' % annotationid
        webpage += 'Thank you for your input<br><br>'
        webpage += '<a href=%s>Back to annotation %d details</a>' % (annotationid, url_for('.annotation_info', annotationid=annotationid))
    else:
        webpage = render_template('done_fail.html', mes='Failed to flag annotation', error=httpRes.text)

    webpage += render_template('footer.html', header_color=get_dbbact_server_color())
    return webpage


def draw_download_fasta_button(annotationid):
    '''
    Draw a button with a link to download the fasta sequences of the annotation

    Parameters
    ----------
    annotationid : int
        the annotationid for which to download the sequences

    Returns
    -------
    webPage : str
        html for the download button with the link to the fasta file download page
    '''
    webPage = '<div style="margin: 20px"><button class="btn btn-default" onclick="location.href=\'%s\';"><i class="glyphicon glyphicon-download-alt"></i> Download FASTA</button></div>' % url_for('.annotation_seq_download', annotationid=annotationid)
    return webPage


def draw_download_button(sequences=None):
    '''
    Draw a button with a link to download the fscores

    Parameters
    ----------
    sequences: list of str or None, optional
        if not None, draw the button, linking to the download f-scores function for these sequences
        if None, don't draw the button

    Returns
    -------
    webPage : str
        html for the download button with the link to the fasta file download page
    '''
    if sequences is None:
        return ''

    webPage = '<div style="margin: 20px">'
    webPage += '<form action="%s" method=post>' % url_for('.download_fscores_sequences_form')
    webPage += '<input type="hidden" id="sequences" name="sequences" value=%s>' % ','.join(sequences)
    webPage += '<input type="submit" value="Download scores">'
    webPage += '</form>'
    webPage += '</div>'
    # webPage = '<div style="margin: 20px"><button class="btn btn-default" onclick="location.href=\'%s\';"><i class="glyphicon glyphicon-download-alt"></i> Download scores</button></div>' % url_for('.download_fscores_sequence', sequences=','.join(sequences))
    return webPage


@Site_Main_Flask_Obj.route('/ontology_info/<string:term>')
def ontology_info(term):
    """
    get the information all studies containing an ontology term (exact or as parent)
    input:
    term : str
        the ontology term to look for
    """
    err, webpage = get_ontology_info(term)
    return webpage


def get_ontology_info(term, show_ontology_tree=False, show_associated_seqs=True):
    """
    get the information all studies containing an ontology term (exact or as parent)
    input:
    term : str
        the ontology term to look for
    show_ontology_tree: bool, optional
        if True, show the term tree graph using cytoscape.js
    show_associated_seq: bool, optional
        if True, show the top term positive/negative-associated sequences
    """
    # get the term annotations
    res = requests.get(get_dbbact_server_address() + '/ontology/get_annotations', params={'term': term, 'get_children': 'true'})
    if res.status_code != 200:
        msg = 'error getting annotations for ontology term %s: %s' % (Markup.escape(term), res.content)
        debug(6, msg)
        return msg, msg
    annotations = res.json()['annotations']
    if len(annotations) == 0:
        debug(1, 'ontology term %s not found' % Markup.escape(term))
        return 'term not found', 'term not found'

    for cannotation in annotations:
        cannotation['website_sequences'] = [0]

    webPage = render_header(title='dbBact taxonomy')
    webPage += '<h1>Summary for ontology term: %s</h1>\n' % Markup.escape(term)
    webPage += 'Number of annotations with term: %d<br>' % len(annotations)

    # plot the top positive/negative associated sequences with the term
    if show_associated_seqs:
        webPage += get_term_seq_scores(term, annotations)

    if show_ontology_tree:
        webPage += draw_term_info(term)
    webPage += '<h2>Annotations:</h2>'
    webPage += draw_annotation_details(annotations)
    webPage += render_template('footer.html')
    return '', webPage


@Site_Main_Flask_Obj.route('/annotations_list')
def annotations_list():
    debug(1, 'annotations_list')
    res = requests.get(get_dbbact_server_address() + '/annotations/get_all_annotations')
    if res.status_code != 200:
        msg = 'error getting annotations list: %s' % res.content
        debug(6, msg)
        return msg, msg
    webPage = render_header(title='dbBact annotation list')
    webPage += '<h2>dbBact Annotation List</h2>'
    annotations = res.json()['annotations']
    for cannotation in annotations:
        cannotation['website_sequences'] = [-1]
    annotations = sorted(annotations, key=lambda x: x.get('date', 0), reverse=True)
    webPage += draw_annotation_details(annotations,include_ratio=False)
    return webPage


@Site_Main_Flask_Obj.route('/experiments_list')
def experiments_list():
    err, webpage = get_experiments_list()
    if err:
        webPageTemp = render_header(title='Error') + render_template('error_page.html', error_str=err)
        return(webPageTemp, 400)
    return webpage


def get_experiments_list():
    '''Get the list of experiments in the database and the details about each one
    Parameters
    ----------

    Returns
    -------
    webpage : str
        the webpage for the experiment list
    '''
    # get the experiments list
    debug(1, 'get_experiments_list')
    res = requests.get(get_dbbact_server_address() + '/experiments/get_experiments_list')
    if res.status_code != 200:
        msg = 'error getting experiments list: %s' % res.content
        debug(6, msg)
        return msg, msg
    explist = res.json().get('explist', [])
    if len(explist) == 0:
        msg = 'no experiments found.'
        debug(3, msg)
        return msg, msg
    webPage = render_header(title='dbBact experiment List')
    webPage += render_template('explist.html')
    for cexp in explist:
        cid = cexp[0]
        cexpname = ''
        for cdetail in cexp[1]:
            cname = cdetail[0]
            cval = cdetail[1]
            if cname != 'name':
                continue
            if len(cexpname) > 0:
                cexpname += '<br>'
            cexpname += cval
        webPage += '<tr><td><a href=%s>%d</a></td>' % (url_for('.experiment_info', expid=cid), cid)
        webPage += '<td>%s</td>' % Markup.escape(cexpname)
        # webPage += '<td><a href=exp_info/' + str(cid) + ">" + str(cid) + "</a></td>"
        # webPage += '<td>' + cval + '</td>'
        webPage += "</tr>"
    webPage += "</table>"
    webPage += render_template('footer.html')
    return '', webPage


@Site_Main_Flask_Obj.route('/taxonomy_info/<string:taxonomy>')
def taxonomy_info(taxonomy):
    '''
    get the information all studies containing any bacteria with taxonomy as substring

    Parameters
    ----------
    taxonomy : str
        the partial taxonomy string to look for

    Returns
    -------
    err : str
        empty ('') if found, none empty if error encountered
    webPage : str
        the html of the resulting table
    '''
    err, webpage = get_taxonomy_info(taxonomy)
    return webpage


def get_taxonomy_info(taxonomy):
    '''
    get the information all studies containing any bacteria with taxonomy as substring

    Parameters
    ----------
    taxonomy : str
        the partial taxonomy string to look for

    Returns
    -------
    err : str
        empty ('') if found, none empty if error encountered
    webPage : str
        the html of the resulting table
    '''
    # get the taxonomy annotations
    debug(2, 'get_taxonomy_info for %s' % taxonomy)
    res = requests.get(get_dbbact_server_address() + '/sequences/get_taxonomy_annotations', json={'taxonomy': taxonomy})
    if res.status_code != 200:
        msg = 'error getting taxonomy annotations for %s: %s' % (Markup.escape(taxonomy), res.content)
        debug(6, msg)
        return msg, msg
    tax_seqs = res.json()['seqids']
    annotations_counts = res.json()['annotations']
    debug(2, 'found %d taxonomy annotations for %d sequences for taxonomy %s' % (len(annotations_counts), len(tax_seqs), taxonomy))
    if len(annotations_counts) == 0:
        msg = 'no annotations found for taxonomy %s' % Markup.escape(taxonomy)
        debug(1, msg)
        return msg, msg

    # convert to list of annotations with counts as a key/value
    annotations = []
    for cann in annotations_counts:
        cannotation = cann[0]
        cannotation['website_sequences'] = [-1] * cann[1]
        annotations.append(cannotation)

    annotations = sorted(annotations, key=lambda x: x.get('num_sequences', 0), reverse=False)
    annotations = sorted(annotations, key=lambda x: len(x.get('website_sequences', [])), reverse=True)

    # get all dbbact sequences containing the taxonomy
    res = requests.get(get_dbbact_server_address() + '/sequences/get_taxonomy_sequences', json={'taxonomy': taxonomy})
    if res.status_code != 200:
        msg = 'error getting taxonomy sequences for %s: %s' % (Markup.escape(taxonomy), res.content)
        debug(6, msg)
        return msg, msg
    seqs = res.json()['sequences']
    debug(2, 'found %d sequences for taxonomy %s' % (len(seqs), taxonomy))

    # add the list of bacterial sequences with the taxonomy
    tax_seq_list = ''
    # sort so highest total experiments sequence is first
    seqs = sorted(seqs, key=lambda k: k['total_experiments'], reverse=True)
    for cseqinfo in seqs:
        cseqinfo['seq'] = cseqinfo['seq'].upper()
        tax_seq_list += "<tr>"
        tax_seq_list += '<td>' + str(cseqinfo['total_experiments']) + '</td>'
        tax_seq_list += '<td>' + str(cseqinfo['total_annotations']) + '</td>'
        tax_seq_list += '<td>' + cseqinfo['taxonomy'] + '</td>'
        tax_seq_list += '<td><a href=%s>%s</a></td>' % (url_for('.sequence_annotations', sequence=cseqinfo['seq']), Markup.escape(cseqinfo['seq']))
        tax_seq_list += '</tr>'

    webPage = render_header(title='dbBact ontology')
    webPage += render_template('taxinfo.html', taxonomy=taxonomy, seq_count=len(tax_seqs), details=tax_seq_list)

    sequences = [x['seq'] for x in seqs]
    webPage += draw_annotation_details(annotations, sequences=sequences)
    webPage += render_template('footer.html')
    return '', webPage


@Site_Main_Flask_Obj.route('/species_info/<string:species>')
def species_info(species):
    '''
    get the information on all sequences matching SILVA sequences with the given species name

    Parameters
    ----------
    species : str
        the SILVA species name (i.e. 'akkermansia muciniphila') to match

    Returns
    -------
    err : str
        empty ('') if found, none empty if error encountered
    webPage : str
        the html of the resulting table
    '''
    err, webpage = get_species_info(species)
    return webpage


def get_species_info(species):
    '''
    get the information all sequences and studies containing any bacteria with SILVA species name

    Parameters
    ----------
    species : str
        the species name to search for

    Returns
    -------
    err : str
        empty ('') if found, none empty if error encountered
    webPage : str
        the html of the resulting table
    '''
    # get the sequences
    debug(2, 'Get species info')
    res = requests.get(get_dbbact_server_address() + '/sequences/get_species_seqs', json={'species': species})
    if res.status_code != 200:
        msg = 'error getting species sequences for %s: %s' % (Markup.escape(species), res.content)
        debug(6, msg)
        return msg, msg
    ids = res.json()['ids']
    debug(2, 'Got %d sequences matching the species %s' % (len(ids), species))
    if len(ids) == 0:
        msg = "No sequences found for species %s" % species
        return msg, msg
    res = requests.get(get_dbbact_server_address() + '/sequences/get_info', json={'seqids': ids})
    debug(2, 'got info')
    if res.status_code != 200:
        msg = 'error getting sequence info: %s' % (res.content)
        debug(6, msg)
        return msg, msg
    seqs = res.json()['sequences']
    ok_seqs = []
    for cseq in seqs:
        if 'total_experiments' in cseq:
            ok_seqs.append(cseq)
    seqs = ok_seqs

    # add the list of bacterial sequences with the taxonomy
    tax_seq_list = ''
    # sort so highest total experiments sequence is first
    seqs = sorted(seqs, key=lambda k: k['total_experiments'], reverse=True)
    # true seqs stores a list of the actual sequences - for plotting the wordcloud etc.
    true_seqs = []
    for cseqinfo in seqs:
        cseqinfo['seq'] = cseqinfo['seq'].upper()
        tax_seq_list += "<tr>"
        tax_seq_list += '<td>' + str(cseqinfo['total_experiments']) + '</td>'
        tax_seq_list += '<td>' + str(cseqinfo['total_annotations']) + '</td>'
        tax_seq_list += '<td>' + cseqinfo['taxonomy'] + '</td>'
        tax_seq_list += '<td><a href=%s>%s</a></td>' % (url_for('.sequence_annotations', sequence=cseqinfo['seq']), Markup.escape(cseqinfo['seq']))
        tax_seq_list += '</tr>'
        true_seqs.append(cseqinfo['seq'])

    webPage = render_header(title='dbBact ontology')
    webPage += render_template('taxinfo.html', taxonomy=species, seq_count=len(ids), details=tax_seq_list)
    err, seq_annotations_compact_part = draw_sequences_annotations_compact(true_seqs, ignore_exp=[], draw_only_details=True)
    if not err:
        webPage += seq_annotations_compact_part
    else:
        webPage += '<br>Sequence details error encountered<br>'
    # webPage += draw_annotation_details(annotations)
    webPage += render_template('footer.html')
    return '', webPage


def get_hash_info(hash_str):
    '''
    get the information about a sequence based on its qiime2 hash

    Parameters
    ----------
    hash_str : string
        sequence represented by its hash

    Returns
    -------
    err : str
        empty ('') if found, none empty if error encountered
    webPage : str
        the html of the resulting table
    '''
    # get the hash annotations
    res = requests.get(get_dbbact_server_address() + '/sequences/get_hash_annotations', json={'hash': hash_str})
    if res.status_code != 200:
        msg = 'error getting hash annotations for %s: %s' % (Markup.escape(hash_str), res.content)
        debug(6, msg)
        return msg, msg
    seq_strs = res.json()['seqstr']
    hash_seqs = res.json()['seqids']
    annotations_counts = res.json()['annotations']
    if len(annotations_counts) == 0:
        msg = 'no annotations found for hash %s' % Markup.escape(hash_str)
        debug(1, msg)
        return msg, msg

    # convert to list of annotations with counts as a key/value
    annotations = []
    for cann in annotations_counts:
        cannotation = cann[0]
        cannotation['website_sequences'] = [-1] * cann[1]
        annotations.append(cannotation)

    annotations = sorted(annotations, key=lambda x: x.get('num_sequences', 0), reverse=False)
    annotations = sorted(annotations, key=lambda x: len(x.get('website_sequences', [])), reverse=True)

    seq_web = ''
    for seq in seq_strs:
        if len(seq_web) > 0:
            seq_web += '<br>'
        seq_web += seq.upper()

    webPage = render_header(title='dbBact ontology')
    webPage += render_template('hashinfo.html', hash_place_holder=hash_str, seq_names_place_holder=seq_web)
    webPage += draw_annotation_details(annotations, sequences=seq_strs)
    webPage += render_template('footer.html')
    return '', webPage


def get_silva_info(silva_str):
    '''
    get the information about a sequence based on its silva id

    Parameters
    ----------
    silva_str : str
        sequence represented by its silva id (should be XXXX.N.NNN or XXXX)

    Returns
    -------
    err : str
        empty ('') if found, none empty if error encountered
    webPage : str
        the html of the resulting table
    '''
    # trim the .XXX.YYY part from the silva id if present
    parts = silva_str.split('.')
    if len(parts) > 1:
        silva_str = '.'.join(parts[:-2])

    # get the annotations
    res = requests.get(get_dbbact_server_address() + '/sequences/get_annotations', json={'sequence': silva_str, 'dbname': 'silva'})
    if res.status_code != 200:
        msg = 'error getting silva annotations for %s: %s' % (Markup.escape(silva_str), res.content)
        debug(6, msg)
        return msg, msg

    annotations = res.json().get('annotations')

    if len(annotations) == 0:
        webPage = render_header(title='dbBact ontology')
        webPage += 'No annotations found for silva sequence %s' % Markup.escape(silva_str)
        webPage += render_template('footer.html')
        return 'no annotations found', webPage

    # get the matching dbbact sequences for the silva id
    res = requests.get(get_dbbact_server_address() + '/sequences/getid', json={'sequence': silva_str, 'dbname': 'silva'})
    ids = res.json().get('seqId')
    res = requests.get(get_dbbact_server_address() + '/sequences/get_info', json={'seqids': ids})

    silva_seq_strs = [x['seq'] for x in res.json()['sequences']]
    silva_seq_tax = [x['taxonomy'] for x in res.json()['sequences']]
    silva_seq_ids = ids
    if len(annotations) == 0:
        msg = 'no annotations found for silva id %s' % Markup.escape(silva_str)
        debug(1, msg)
        return msg, msg

    seq_web = ''
    total_num_seqs = len(silva_seq_strs)
    for idx, seq in enumerate(silva_seq_strs):
        if len(seq_web) > 0:
            seq_web += '\n'
        seq_web += seq.upper()
        seq_web += '\n'
        seq_web += '(' + silva_seq_tax[idx].lower() + ')\n'

    webPage = render_header(title='dbBact ontology')
    webPage += render_template('silvainfo.html', silva_place_holder=silva_str, seq_names_place_holder=seq_web, number_seqs_place_holder=total_num_seqs)
    webPage += render_sequence_annotations(annotations)
    # webPage += draw_annotation_details(annotations)
    webPage += render_template('footer.html')
    return '', webPage


@Site_Main_Flask_Obj.route('/exp_info/<int:expid>')
def experiment_info(expid):
    """
    get the information about a given study dataid
    input:
    dataid : int
        The dataid on the study (DataID field)

    output:
    info : list of (str,str,str)
        list of tuples for each entry in the study:
        type,value,descstring about dataid
        empty if dataid not found
    """

    # get the experiment details
    webPage = render_header()
    res = requests.get(dbbact_server_address + '/experiments/get_details', json={'expId': expid})
    if res.status_code == 200:
        webPage += draw_experiment_info(expid, res.json()['details'])
    else:
        message = 'Experiment ID:%s was not found.' % expid
        return(render_header(title='Not found') +
               render_template('error.html', title='Not found',
                               message=message) +
               render_template('footer.html'), 400)

    # get the experiment annotations
    res = requests.get(dbbact_server_address + '/experiments/get_annotations', json={'expId': expid})
    annotations = res.json()['annotations']
    for cannotation in annotations:
        cannotation['website_sequences'] = [-1]
    annotations = sorted(annotations, key=lambda x: x.get('num_sequences', 0), reverse=True)
    webPage += '<h2>Annotations for experiment:</h2>'
    webPage += draw_annotation_details(annotations, include_word_cloud=False, include_ratio=False)
    webPage += render_template('footer.html')
    return webPage


def draw_experiment_info(expid, exp_details):
    '''
    Draw the table with all experiment details

    Parameters
    ----------
    expid : int
        the experiment id
    exp_details : list of (str,str)
        list of experiment detail name and type ('details' from REST API /experiments/get_details/ )

    Returns
    -------
    webPage : str
        the html of the experiment info table
    '''
    webPage = ""
    webPageTemp = render_template('expinfo.html', expid=expid)

    lastExp = ""
    firstExp = ""
    lastExp = ""
    for cres in exp_details:
        if cres[0].upper() == "NAME":
            firstExp += "<tr>"
            firstExp += '<td>' + cres[0] + '</td>'
            firstExp += '<td>' + cres[1] + '</td><tr>'
        elif cres[0].upper().find("MD5") > -1:
            lastExp += "<tr>"
            lastExp += '<td>' + cres[0] + '</td>'
            lastExp += '<td>' + cres[1] + '</td><tr>'
        else:
            webPage += "<tr>"
            webPage += '<td>' + cres[0] + '</td>'
            webPage += '<td>' + cres[1] + '</td><tr>'
    webPage = webPageTemp + firstExp + webPage + lastExp
    webPage += '</table>'
    return webPage


@Site_Main_Flask_Obj.route('/annotation_seqs/<int:annotationid>')
def annotation_seqs(annotationid):
    '''
    get the information about all sequences in a given annotation
    input:
    annotationid : int
        The annotation for which to show the sequence info

    returns
    -------
    webPage : str
        the html page for the annotation sequences
    '''

    # get the annotation details
    res = requests.get(dbbact_server_address + '/annotations/get_annotation', json={'annotationid': annotationid})
    if res.status_code != 200:
        msg = 'Error encountered when getting info for annotation ID %d: %s' % (annotationid, res.content)
        debug(6, msg)
        return(render_header(title='Not found') +
               render_template('error.html', title='Not found',
                               message=msg) +
               render_template('footer.html'), 600)
    annotation = res.json()
    shortdesc = getannotationstrings(annotation)
    webPage = render_header()
    webPage += render_template('annotseqs.html', annotationid=annotationid)
    webPage += '<div style="margin: 20px;"><blockquote style="font-size: 1em;"><p>%s</p></blockquote></div>\n' % Markup.escape(shortdesc)

    webPage += '<h2>Download</h2>'
    webPage += draw_download_fasta_button(annotationid)

    # get the sequence information for the annotation
    res = requests.get(dbbact_server_address + '/annotations/get_full_sequences', json={'annotationid': annotationid})
    if res.status_code != 200:
        msg = 'Error encountered when getting sequences for annotation ID %d: %s' % (annotationid, res.content)
        debug(6, msg)
        return(render_header(title='Not found') +
               render_template('error.html', title='Not found',
                               message=msg) +
               render_template('footer.html'), 600)
    sequences = res.json()['sequences']
    webPage += draw_sequences_info(sequences)

    return webPage


def draw_sequences_info(sequences):
    '''write the table entries for each sequence (sequence, total counts etc.)
    '''
    webPage = render_template('seqlist.html')
    # sort the sequences based on taxonomy
    sequences = sorted(sequences, key=lambda x: x.get('taxonomy', ''))
    for cseqinfo in sequences:
        cseqinfo['seq'] = cseqinfo['seq'].upper()
        webPage += "<tr>"
        webPage += '<td>' + cseqinfo['taxonomy'] + '</td>'
        webPage += '<td><a href=%s>%s</a></td>' % (url_for('.sequence_annotations', sequence=cseqinfo['seq']), Markup.escape(cseqinfo['seq']))
        webPage += '<td>' + 'na' + '</td><tr>'
    webPage += '</table>'
    return webPage


@Site_Main_Flask_Obj.route('/forgot_password_submit', methods=['POST', 'GET'])
def forgot_password_submit():
    """
    this page will send the forgoten password to the user via mail. It is called from the reset_password page when submitting the username/email
    input:
    dataid : string
        user name or email

    output:
    """
    debug(3, 'forgot_password_submit')
    usermail = ''
    if request.method == 'GET':
        usermail = request.args['useremail']
    else:
        usermail = request.form['useremail']

    json_user = {'user': usermail}
    debug(3, 'posting to dbbact-server forgot_password')
    httpRes = requests.post(dbbact_server_address + '/users/forgot_password', json=json_user)
    debug(3, 'got result')
    if httpRes.status_code == 200:
        webpage = render_header(title='Password Recovery')
        webpage += render_template('recover_form.html')
    else:
        webpage = render_template('done_fail.html', mes='Failed to reset password', error=httpRes.text)
    return webpage


@Site_Main_Flask_Obj.route('/change_password', methods=['POST', 'GET'])
def change_password():
    """Ask for username, verification code and new password
    and submit the change password request to dbBact
    """
    debug(3, 'change_password')
    webpage = render_header(title='Password Recovery')
    webpage += render_template('recover_form.html')
    return webpage


@Site_Main_Flask_Obj.route('/recover_user_password', methods=['POST', 'GET'])
def recover_user_password():
    """
    this function will update new user password in the db
    input:
    dataid : string
        user email

    output:
    """
    debug(3, 'recover_user_password')
    usermail = ''
    if request.method == 'GET':
        usermail = request.args['user']
        recoverycode = request.args['recoverycode']
        newpassword = request.args['newpassword']
    else:
        usermail = request.form['user']
        recoverycode = request.form['recoverycode']
        newpassword = request.form['newpassword']

    json_user = {}
    json_user['user'] = usermail
    json_user['recoverycode'] = recoverycode
    json_user['newpassword'] = newpassword

    httpRes = requests.post(dbbact_server_address + '/users/recover_password', json=json_user)
    if httpRes.status_code == 200:
        debug(3, 'recover_password for user %s succeeded' % usermail)
        webpage = render_template('update_password_success.html')
    else:
        debug(5, 'recover password for user %s failed' % usermail)
        webpage = render_template('done_fail.html', mes='Failed to reset password', error=httpRes.text)
    return webpage


@Site_Main_Flask_Obj.route('/user_info/<string:username>')
def user_info(username):
    """
    get the information about a user
    input:
    dataid : int
        the user id

    output:
    """
    rdata = {}
    rdata['username'] = username
    if username is None:
        return "Error: Invalid user"

    debug(1, 'get user info for user %s' % username)
    # get the experiment details
    httpRes = requests.post(dbbact_server_address + '/users/get_user_public_information', json=rdata)
    if httpRes.status_code == 200:
        userInfo = httpRes.json()
        username = userInfo.get('name')
        name = userInfo.get('username', 'NA')
        if len(username) == 0:
            username = 'NA'
        desc = userInfo.get('description')
        if desc is None:
            desc = 'NA'
        email = userInfo.get('email', 'NA')
        userid = userInfo.get('id', 0)
        total_annotations = 'NA'
        total_experiments = 'NA'

        # get user annotation
        forUserId = {'foruserid': userid}
        httpRes = requests.get(dbbact_server_address + '/users/get_user_annotations', json=forUserId)
        if httpRes.status_code == 200:
            annotations = httpRes.json().get('userannotations')
            total_annotations = len(annotations)
            exps = set()
            for cannotation in annotations:
                if 'expid' in cannotation:
                    exps.add(cannotation['expid'])
            total_experiments = len(exps)
        else:
            annotations = None
        webPage = render_header(title=username)
        webPage += render_template('userinfo.html', userid=userid, name=name, username=username, desc=desc, email=email, total_annotations=total_annotations, total_experiments=total_experiments)

        if annotations is not None:
            webPage += draw_annotation_details(annotations)
        webPage += render_template('footer.html')
        return webPage
    else:
        message = Markup('Failed to get user information:<br><br><blockquote><code>%s</code></blockquote>'
                         % httpRes.content)
        return(render_header(title='Not found') +
               render_template('error.html', title='Not found',
                               message=message) +
               render_template('footer.html'))


def draw_annotation_details(annotations, seqannotations=None, term_info=None, show_relative_freqs=False, include_word_cloud=True, include_ratio=True, ignore_exp=[], sequences=None):
    '''Draw the wordcloud and details for a list of annotations
    Converts the annotations list to dict, creates the seqannotations and calls draw_group_annotation_details()

    Parameters
    ----------
    annotations: list of dict
        list of annotations (each annotation is a dict of annotation details) - from dbbact /sequences/get_annotations
    seqannotations: list of tuples of (seqid (int), list [annotationids]) or None
        if None, ????
    term_info : dict of dict or None (optional)
        None (default) to skip relative word cloud.
        Otherwise need to have information about all ontology terms to be drawn
        dict of {term: dict} where
            term : ontology term (str)
            dict: pairs of:
                'total_annotations' : int
                'total_sequences' : int
    include_word_cloud: bool, optional
        True (default) to draw the wordcloud, False to just draw the term tables
    ignore_exp: list of int, optional
        the experiment ids to exclude from the analysis
    sequences: list of str or None, optional
        if not None, this is the sequence the annotations are for (used for the download annotations button)

    Returns
    -------
    html part for wordcloud and term tables
    '''
    if sequences is None:
        num_seqs = 0
    else:
        num_seqs = len(sequences)
    if annotations is None:
        num_anno = 0
    else:
        num_anno = len(annotations)
    debug(2, 'draw_annotation_details for %d sequences, %d annotations' % (num_seqs, num_anno))
    annotations_dict = {}
    for cannotation in annotations:
        annotations_dict[str(cannotation['annotationid'])] = cannotation
    if seqannotations is None:
        debug(2, 'no seqannotations provided, querying')
        if sequences is not None:
            res = requests.get(get_dbbact_server_address() + '/sequences/get_fast_annotations', json={'sequences': sequences})
            if res.status_code != 200:
                msg = 'error getting annotations for sequences : %s' % Markup.escape(res.content)
                debug(6, msg)
                return 'Error %s encountered' % msg
            res = res.json()
            orig_seqannotations = res['seqannotations']
            # now make sure we only keep the sequence annotations found in our annotations
            seqannotations = []
            for cseqid, cseqanno in orig_seqannotations:
                ok_anno = [str(x) for x in cseqanno if str(x) in annotations_dict]
                seqannotations.append( (cseqid, ok_anno) )
        else:
            seqannotations = (((0, list(annotations_dict.keys())),))
    wpart = draw_group_annotation_details(annotations_dict, seqannotations=seqannotations, term_info=term_info, include_word_cloud=include_word_cloud, ignore_exp=ignore_exp, sequences=sequences)
    return wpart


def draw_wordcloud_fscore(fscores, recall=None, precision=None, term_count={}):
    '''Draw the wordcloud for the terms in fscores and return an html section with the image embedded

    Parameters
    ----------
    fscores: dict of {term(str): fscore(float)}
    recall: dict of {term(str): recall score(float)}, optional
    precision: dict of {term(str): precision score(float)}, optional
    term_count: dict of {term(str): number of experiments with term(float)}, optional
        used to determine the color intensity

    Returns
    -------
    str: the html part with embedded wordcloud image
    '''
    wordcloud_image = draw_cloud(fscores, recall, precision, term_count)
    wordcloudimage = urllib.parse.quote(wordcloud_image)
    wpart = ''
    if wordcloudimage:
        wpart += render_template('wordcloud.html', wordcloudimage=urllib.parse.quote(wordcloud_image))
    else:
        wpart += '<p></p>'
    return wpart


def draw_annotation_table(annotations, include_ratio=True):
    '''Draw the annotations list table. Note annotations are written according to the list order

    It uses the annottable.html template for the table.
    If the 'website_sequences' is filled in each annotation, it also writes it (XXX/YYY) where XXX is the 'website_sequences'

    Parameters
    ----------
    annotations: list of annotations(dict)
        The annotations to list. Can contain the 'website_sequences' to indicate how many sequence were found for this annotation
        in the query sequence list.

    Returns
    -------
    str: the HTML part for the annotations table
    '''
    wpart = ''
    for dataRow in annotations:
        wpart += '  <tr>'
        # add the experimentid info+link
        expid = dataRow.get('expid', 'not found')
        wpart += "<td><a href=%s>%s</a></td>" % (url_for('.experiment_info', expid=expid), Markup.escape(expid))

        # add user name+link
        userid = dataRow.get('userid', 'not found')
        username = dataRow.get('username', 'not found')
        wpart += "<td><a href=%s>%s</a></td>" % (url_for('.user_info', username=username), Markup.escape(username))

        # add the annotation description
        cdesc = getannotationstrings(dataRow)
        annotationid = dataRow.get('annotationid', -1)
        wpart += "<td><a href=%s>%s</a></td>" % (url_for('.annotation_info', annotationid=annotationid), Markup.escape(cdesc))

        # add the annotation date
        wpart += '<td>%s</td>' % Markup.escape(dataRow['date'])

        # add the annotation region
        wpart += '<td>%s</td>' % Markup.escape(dataRow['primer'])

        # add the sequences
        annotationid = dataRow.get('annotationid', -1)
        num_sequences = dataRow.get('num_sequences', '?')
        if 'website_sequences' in dataRow:
            observed_sequences = len(dataRow['website_sequences'])
            if include_ratio is True:
                sequences_string = '%s / %s' % (observed_sequences, num_sequences)
            else:
                sequences_string = '%s' % (num_sequences)
        else:
            observed_sequences = '?'
            sequences_string = '%s' % num_sequences
        wpart += "<td><a href=%s>%s</a></td>" % (url_for('.annotation_seqs', annotationid=annotationid), Markup.escape(sequences_string))

        # add the annotation review status
        review_status = dataRow.get('review_status', -1)
        if review_status == 0:
            review_status_str = 'pending'
        elif review_status == 1:
            review_status_str = 'approved'
        else:
            review_status_str = 'NA'
        wpart += '<td>%s</td>' % Markup.escape(review_status_str)

        # add the number of flags
        if len(dataRow['flags']) > 0:
            flags = '%s' % len(dataRow['flags'])
        else:
            flags = 'No'
        wpart += '<td>%s</td>' % Markup.escape(flags)

        wpart += '</tr>\n'
    # wpart += '</table>\n'
    # wpart += '</div>\n'
    # return wpart
    ppart = render_template('annottable.html', details=wpart)
    return ppart


def draw_ontology_score_list(scores, section_id, description=None, max_terms=100, term_set=None):
    '''Create table entries for ontology terms sorted by score

    Parameters
    ----------
    scores: dict of {term (str): score (float)}
        the scores to sort the terms by
    section_id: str
        the name of the section (for the tabs - i.e. 'recall' etc.)
    max_terms: int, optional
        maximal number of terms to add ot the list, or None to show all
    term_set: set of str, optional
        if not None, only show terms in this set

    Returns
    -------
    str
        the html code for the section
    list of (str, float)
        the terms and scores sorted by score (after using the max_terms and term_set filters)
    '''
    # wpart = '<div id="%s" class="tab-pane" style="margin-top: 20px; margin-bottom: 20px;">\n' % section_id
    wpart = '<div id="%s" class="tab-pane" role="tabpanel" style="margin-top: 20px; margin-bottom: 20px;">\n' % section_id
    if description is not None:
        wpart += description

    wpart += '<table style="width: 90%;">\n'
    wpart += '<col><col width="100px">\n'
    wpart += '<tr><th>Term</th><th>Score</th></tr>\n'

    data = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    if max_terms is not None:
        data = data[:max_terms]
    
    if term_set is not None:
        data = [x for x in data if x[0] in term_set]

    for cterm, cscore in data:
        if len(cterm) == 0:
            debug(3, 'Empty term in ontology score list')
            continue
        if cterm[0] == '-':
            ctermlink = cterm[1:]
            cterm = 'LOWER IN %s' % Markup.escape(ctermlink)
        else:
            ctermlink = cterm
        wpart += '<tr><td><a href=%s>%s</a></td><td>%f</td></tr>\n' % (url_for('.ontology_info', term=ctermlink), Markup.escape(cterm), cscore)
    wpart += '</table>\n'
    wpart += '</div>\n'
    return wpart, data


@Site_Main_Flask_Obj.route('/annotation_seq_download/<int:annotationid>')
def annotation_seq_download(annotationid):
    '''return a download of the sequences of the annotation as fasta
    '''
    # get the experiment annotations
    res = requests.get(get_dbbact_server_address() + '/annotations/get_full_sequences', json={'annotationid': annotationid})
    annotation = res.json()
    seqs = annotation.get('sequences')
    if seqs is None:
        debug(6, 'No sequences found')
        webPageTemp = render_header(title='Error') + render_template('error_page.html', error_str='No sequences found')
        return(webPageTemp, 400)
    output = ''
    for idx, cseq in enumerate(seqs):
        output += str(Markup.escape('>%d %s\n%s\n' % (idx, cseq.get('taxonomy', ''), cseq['seq'])))
    response = make_response(output)
    response.headers["Content-Disposition"] = "attachment; filename=annotation-%d-sequences.fa" % annotationid
    return response


@Site_Main_Flask_Obj.route('/download_fscores_sequences_form',methods=['POST'])
def download_fscores_sequences_form():
    debug(1,'download_fscores_sequences_form')
    sequences = request.form['sequences']
    sequences = sequences.split(',')
    debug(2,'got %d sequences for download_fscores_sequences_form' % len(sequences))

    res = requests.get(get_dbbact_server_address() + '/sequences/get_fast_annotations',
                       json={'sequences': sequences})
    if res.status_code != 200:
        msg = 'error getting annotations for sequences : %s' % Markup.escape(res.content)
        debug(6, msg)
        return msg, msg

    res = res.json()
    annotations = res['annotations']
    seqannotations = res['seqannotations']
    if len(seqannotations) == 0:
        msg = 'None of the %d sequences were found in dbBact. Are these >100bp long 16S sequences?\nNote dbBact is populated mostly by EMP V4 (515F) amplicon sequences.' % len(seqs)
        debug(3, msg)
        return msg, msg
    term_info = res['term_info']

    # get the experiment annotations
    ignore_exp = []
    dbc = dbbact_calour.dbbact.DBBact(dburl=get_dbbact_server_address(), test_version=False)
    fscores, recall, precision, term_count, reduced_f = dbc.get_enrichment_score(annotations, seqannotations, ignore_exp=ignore_exp, term_info=term_info)

    output = 'term\tf-score\trecall\tprecision\tcount\n'
    for cterm,cfscore in reduced_f.items():
        output += '%s\t%s\t%s\t%s\t%s\n' % (cterm, cfscore, recall.get(cterm, 0), precision.get(cterm, 0), term_count.get(cterm, 0))
    debug(1, output)
    response = make_response(output)
    response.headers["Content-Disposition"] = "attachment; filename=sequence-scores.tsv"
    return response


def _get_color(word, font_size, position, orientation, font_path, random_state, fscore, recall, precision, term_count):
    '''Get the color for a wordcloud term based on the term_count and higher/lower

    If term starts with "-", it is lower in and colored red. otherwise colored blue
    If we have term_count, we use it to color from close to white(low count) to full color (>=10 experiments)

    Parameters
    ----------
    fscores: dict of {term(str): fscore(float)}
        between 0 and 1
    recall: dict of {term(str): recall score(float)}, optional
    precision: dict of {term(str): precision score(float)}, optional
    term_count: dict of {term(str): number of experiments with term(float)}, optional
        used to determine the color intensity


    Returns
    -------
    str: the color in hex "0#RRGGBB"
    '''
    if word in term_count:
        count = min(term_count[word], 10)
    else:
        count = 10

    if word[0] == '-':
        cmap = mpl.cm.get_cmap('Oranges')
        rgba = cmap(float(0.4 + count / 40), bytes=True)
        rgba = cmap(float(0.3 + count / 20), bytes=True)
    else:
        cmap = mpl.cm.get_cmap('Purples')
        rgba = cmap(float(0.4 + count / 40), bytes=True)
        rgba = cmap(float(0.3 + count / 20), bytes=True)

    red = format(rgba[0], '02x')
    green = format(rgba[1], '02x')
    blue = format(rgba[2], '02x')
    return '#%s%s%s' % (red, green, blue)


def draw_cloud(fscores, recall={}, precision={}, term_count={}, local_save_name=None):
    '''
    Draw a wordcloud for a list of terms

    Parameters
    ----------
    fscores: dict of {term (str): fscore(float)}
        the f-score for each term (determines the size)
    recall : dict of {term (str): recall(float)}
        the recall score for each term (determines the green color)
    precision : dict of {term (str): precision(float)}
        the precision score for each term (determines the blue color)
    term_count: dict of {term (str): number of experiments where the term was observed(int)}
        The number of unique experiments where the term annotations are coming from
    local_save_name : str or None (optional)
        if str, save also as pdf to local file local_save_name (for MS figures) in high res

    Returns
    -------
    BytesIO file with the image
    '''
    from wordcloud import WordCloud
    import matplotlib.pyplot as plt

    debug(1, 'draw_cloud for %d words' % len(fscores))
    if len(fscores) == 0:
        debug(3, 'no words for wordcloud')
        return ''

    # normalize the fractions to a scale max=1
    # also sanitize the terms to remove new lines (wordcloud doesn't like them)
    new_scores = {}
    if fscores is not None:
        maxval = max(fscores.values())
        debug(1, 'normalizing fscores. maxval is %f' % maxval)
        for ckey, cval in fscores.items():
            sanitized_ckey = ckey.replace('\n', ' ')
            new_scores[sanitized_ckey] = fscores[ckey] / maxval
    fscores = new_scores

    # wc = WordCloud(background_color="white", relative_scaling=0.5, stopwords=set(),colormap="Blues")
    if local_save_name is not None:
        wc = WordCloud(width=400 * 3, height=200 * 3, background_color="white", relative_scaling=0.5, stopwords=set(), color_func=lambda *x, **y: _get_color(*x, **y, fscore=fscores, recall=recall, precision=precision, term_count=term_count))
    else:
        wc = WordCloud(background_color="white", relative_scaling=0.5, stopwords=set(), color_func=lambda *x, **y: _get_color(*x, **y, fscore=fscores, recall=recall, precision=precision, term_count=term_count))

    if isinstance(fscores, str):
        debug(2, 'generating from words list')
        wordcloud = wc.generate(fscores)
    elif isinstance(fscores, dict):
        debug(2, 'generating from frequency dict')
        wordcloud = wc.generate_from_frequencies(fscores)
    else:
        debug(4, 'unknown type for generate_wordcloud!')
        return None

    fig = plt.figure()
    plt.imshow(wordcloud)
    plt.axis("off")
    fig.tight_layout()
    figfile = BytesIO()
    if local_save_name is not None:
        debug(1, 'saving wordcloud to local file %s' % local_save_name)
        fig.savefig(local_save_name, format='pdf', bbox_inches='tight')
    fig.savefig(figfile, format='png', bbox_inches='tight')
    figfile.seek(0)  # rewind to beginning of file
    figdata_png = base64.b64encode(figfile.getvalue())
    figfile.close()
    plt.close()
    return figdata_png


def get_term_info_for_annotations(annotations):
    '''
    Get the statistics about each term in annotations

    Parameters
    ----------
    annotations: list of annotations

    Returns
    -------
    term_info: dict of XXX
    The statistics about each term
    '''
    terms = get_annotations_terms(annotations)
    res = requests.get(get_dbbact_server_address() + '/ontology/get_term_stats', json={'terms': terms})
    if res.status_code != 200:
        debug(6, 'error encountered in get_term_stats: %s' % res.reason)
        return []
    ans = res.json()
    term_info = ans.get('term_info')
    return term_info


@Site_Main_Flask_Obj.route('/reset_password', methods=['POST', 'GET'])
def reset_password():
    """
    Title: Reset password via mail
    URL: /reset password
    Method: POST
    """
    webpage = render_header(title='Reset Password')
    webpage += render_template('reset_password.html')
    return webpage


@Site_Main_Flask_Obj.route('/about', methods=['POST', 'GET'])
def about():
    """
    Title: About us
    URL: /about
    Method: POST
    """
    webpage = render_header(title='About Us')
    webpage += render_template('about.html')
    return webpage


"""
Auto complete tests
"""
@Site_Main_Flask_Obj.route('/add_annotation', methods=['POST', 'GET'])
def add_data():
    """
    Title: Add a new annotation to dbBact. ALPHA!!!
    URL: /about
    Method: POST
    """

    # res = requests.get(get_dbbact_server_address() + '/ontology/get_all_descriptions')
    res = requests.get(get_dbbact_server_address() + '/ontology/get_all_terms')
    if res.status_code != 200:
        debug(6, 'failed to get list of ontologies')
        list_of_ont = {}
    else:
        import json
        list_of_ont = json.dumps(res.json())

    res = requests.get(get_dbbact_server_address() + '/ontology/get_all_synonyms')
    if res.status_code != 200:
        debug(6, 'failed to get list of synonyms')
        list_of_synonym = {}
    else:
        import json
        list_of_synonym = json.dumps(res.json())

    webpage = render_header()
    webpage += render_template('add_data.html', syn_list=list_of_synonym, ont_list=list_of_ont, display='{{display}}', group='{{group}}', query='{{query}}')
    return webpage


"""
Auto complete tests
"""
@Site_Main_Flask_Obj.route('/add_data2', methods=['POST', 'GET'])
def add_data2():
    """
    Title: About us
    URL: /about
    Method: POST
    """
    
    res = requests.get(get_dbbact_server_address() + '/ontology/get_all_terms')
    # res = requests.get(get_dbbact_server_address() + '/ontology/get_all_descriptions')
    if res.status_code != 200:
           debug(6, 'failed to get list of ontologies')
           parents = []
    else:
           import json
           list_of_ont = json.dumps(res.json())
    
    res = requests.get(get_dbbact_server_address() + '/ontology/get_all_synonyms')
    if res.status_code != 200:
           debug(6, 'failed to get list of synonyms')
           parents = []
    else:
           import json
           list_of_synonym = json.dumps(res.json())
    
    webpage = render_header()
    webpage += render_template('add_data2.html',syn_list=list_of_synonym,ont_list=list_of_ont,display='{{display}}',group='{{group}}',query='{{query}}')
    return webpage


"""
Auto complete tests
"""
@Site_Main_Flask_Obj.route('/auto_complete_test', methods=['POST', 'GET'])
def auto_complete_test():
    """
    Title: About us
    URL: /about
    Method: POST
    """
    
    res = requests.get(get_dbbact_server_address() + '/ontology/get_all_descriptions')
    if res.status_code != 200:
           debug(6, 'failed to get list of ontologies')
           parents = []
    else:
           import json
           list_of_ont = json.dumps(res.json())
    
    res = requests.get(get_dbbact_server_address() + '/ontology/get_all_synonyms')
    if res.status_code != 200:
           debug(6, 'failed to get list of synonyms')
           parents = []
    else:
           import json
           list_of_synonym = json.dumps(res.json())

    webpage = render_template('demo-autocomplete.html',syn_list=list_of_synonym,ont_list=list_of_ont,display='{{display}}',group='{{group}}',query='{{query}}')
    return webpage


def error_message(title, message):
    '''
    '''
    return(render_header(title=title) +
           render_template('error.html', title=title,
                           message=message) +
           render_template('footer.html'))


# Handle old dbbact rest-api requests (which are now supposed to go to api.dbbact.org instead)
@Site_Main_Flask_Obj.route('/REST-API/<path:path>', methods=['POST', 'GET'])
def old_dbbact(path):
    res = {}
    annotation1 = {}
    annotation1['description'] = 'ERROR - please update dbbact-calour client'
    annotation1['annotationtype'] = 'other'
    annotation1['details'] = ['all', 'na']
    annotation2 = {}
    annotation2['description'] = 'use pip install --upgrade --force-reinstall --no-deps git+git://github.com/amnona/dbbact-calour'
    annotation2['annotationtype'] = 'other'
    annotation2['details'] = ['all', 'na']
    res['annotations'] = [annotation1, annotation2]
    res['term_info'] = {'error': {'total_sequences': 1, 'total_annotations': 1}}
    return json.dumps(res)


def draw_group_annotation_details(annotations, seqannotations, term_info, include_word_cloud=True, ignore_exp=[], local_save_name=None, sequences=None):
    '''
    Create wordcloud and table entries for a list of annotations

    Parameters
    ----------
    annotations : dict of {annotationid(str): annotationdetails(dict)}
        dict of annotation details (from REST API). NOTE: key is str of annotationID!
    seqannotations: list of (seqid, [annotation ids])
        list of annotations for each sequence.
    term_info : dict of dict or None, optional
        None (default) to skip relative word cloud.
        Otherwise need to have information about all ontology terms to be drawn
        dict of {term: dict} where
            term : ontology term (str)
            dict: pairs of:
                'total_annotations' : int
                'total_sequences' : int
    include_word_cloud: bool (optional)
        True to plot the wordcloud. False to not plot it
    ignore_exp : list of int (optional)
        list of experiment ids to ignore when calculating the score.
    sequences : list of str or None (optional)
        the sequences for which the annotation details are obtained.
        ??? If None, calculate from seqannotations

    Returns
    -------
    wpart : str
        html code for the annotations table
    '''
    # The output webpage part
    wpart = ''

    # calculate the score for each term
    if seqannotations is None:
        num_seqanno = 0
    else:
        num_seqanno = len(seqannotations)
    if sequences is None:
        num_seqs = 0
    else:
        num_seqs = len(sequences)
    if annotations is None:
        num_anno = 0
    else:
        num_anno = len(annotations)
    debug(2, 'calculating fscore using %d annotations, %d seqannotations, ignore_exp=%s and %d sequences' % (num_anno, num_seqanno, ignore_exp, num_seqs))
    dbc = dbbact_calour.dbbact.DBBact(dburl=get_dbbact_server_address(), test_version=False)
    fscores, recall, precision, term_count, reduced_f = dbc.get_enrichment_score(annotations, seqannotations, ignore_exp=ignore_exp, term_info=term_info)

    # draw the wordcloud for the group terms
    if include_word_cloud is True:
        debug(2, 'drawing term pair word cloud')
        # wpart += draw_wordcloud_fscore(fscores, recall, precision, term_count)
        wpart += draw_wordcloud_fscore(reduced_f, recall, precision, term_count)
        wpart += draw_download_button(sequences=sequences)

    wpart += render_template('tabs.html')

    tterms = set()
    cpart, tscores = draw_ontology_score_list(fscores, section_id='fscores', description='term enrichment score')
    wpart += cpart
    tterms.update(cterm for cterm,score in tscores)
    cpart, tscores = draw_ontology_score_list(recall, section_id='recall', description='Fraction of dbbact annotations with this term covered by the query')
    wpart += cpart
    tterms.update(cterm for cterm,score in tscores)
    cpart, tscores = draw_ontology_score_list(precision, section_id='precision', description='Fraction of annotations for the query sequences containing the term')
    wpart += cpart
    tterms.update(cterm for cterm,score in tscores)
    # add the term count tab, only for the terms listed in the fscore/precision/recall parts
    cpart, tscores = draw_ontology_score_list(term_count, section_id='term_count', description='Number of experiments associating the term to the sequence',term_set=tterms, max_terms=None)
    wpart += cpart

    # draw the annotation table
    # first we need to sort the annotations by total sequences from submitted list in each annotation
    annotation_count = defaultdict(int)
    annotation_seq_list = defaultdict(list)
    for cseq, cseq_annotations in seqannotations:
        for cannotation in cseq_annotations:
            annotation_count[cannotation] += 1
            annotation_seq_list[cannotation].append(cseq)
    normalized_annotation_count = {}
    for x in annotation_count:
        normalized_annotation_count[x] = annotation_count[x] / (annotations[str(x)]['num_sequences'] + 50)
    sorted_annotation_count = sorted(annotation_count, key=normalized_annotation_count.get, reverse=True)
    sorted_annotations = [annotations[str(x)] for x in sorted_annotation_count]
    # set up the website_sequences field so we'll see XXX/YYY in the table
    for cannotation, cseqlist in annotation_seq_list.items():
        annotations[str(cannotation)]['website_sequences'] = cseqlist
    wpart += draw_annotation_table(sorted_annotations)

    # draw the ontology term list
    # wpart += draw_ontology_list(sorted_annotations, term_info, fscores)

    wpart += '    </div>\n'
    wpart += '  </div>\n'
    return wpart


def get_alert_text(filename='dbbact_website_alert.txt'):
    try:
        with open(filename) as fl:
            alert_text = fl.readlines()
            if len(alert_text) > 0:
                alert_text = [x.strip() for x in alert_text]
                alert_text = ';'.join(alert_text)
            else:
                alert_text = ''
    except:
        with open(filename, 'w') as fl:
            debug(5, 'empty alert file created: %s' % filename)
        alert_text = ''
    return alert_text


def render_header(title='dbBact', alert_text=True, **kwargs):
    '''Render the header template

    Parameters
    ----------
    title: str, optional
        title of the page
    alert_text: str or bool, optional
        if True, display an alert text if the "DBBACT_WEBSITE_ALERT_MSG" is set and not empty.
        if False, never show an alert
        if str, show the alert provided in str

    Returns
    -------
    str: the rendered html header
    '''
    if alert_text is True:
        alert_text = get_alert_text()
    elif alert_text is False:
        alert_text = ''

    webPageTemp = render_template('header.html', header_color=get_dbbact_server_color(), title=title, alert_text=alert_text)
    return webPageTemp


@Site_Main_Flask_Obj.route('/term_info/<string:term>')
def draw_term_info(term):
    """
    get the information about a given term (either term name (i.e. 'feces' or term_id (i.e. 'gaz:000001')))

    Parameters
    ----------
    term: str
        the term_id (i.e. 'gaz:000001') to get the info about
    show_ontology_tree: bool, optional
        if True, draw the cytoscape.js ontology tree graph
    """
    rdata = {}
    rdata['terms'] = [term]
    rdata['relation'] = 'both'
    if term is None:
        return "Error: Invalid term"

    debug(1, 'get term info for term %s' % term)
    # get the experiment details
    httpRes = requests.get(dbbact_server_address + '/ontology/get_family_graph', json=rdata)
    if httpRes.status_code == 200:
        termInfo = httpRes.json()['family']
        if len(termInfo['nodes']) == 0:
            return 'no results for term %s' % term
        # get list of all parent terms
        all_terms = []
        for cnode in termInfo['nodes']:
            if 'name' in cnode:
                all_terms.append(cnode['name'])
        # get the info about the number of annotations/experiments/sequences per term
        rdata = {'terms': all_terms}
        httpRes = requests.get(dbbact_server_address + '/ontology/get_term_stats', json=rdata)
        if httpRes.status_code == 200:
            term_counts = httpRes.json()['term_info']
        else:
            term_counts = {}

        # create the node and edge data for the cytoscpae graph
        nodes = ''
        for cnode in termInfo['nodes']:
            cnode = cnode.copy()
            cnumexp = 1
            cnumanno = 1
            if 'name' in cnode:
                if cnode['name'] in term_counts:
                    cnumexp += term_counts[cnode['name']]['total_experiments']
                    cnumanno += term_counts[cnode['name']]['total_annotations']
            else:
                cnode['name'] = 'NA'
            if cnumexp > 50:
                cnumexp = 50
            cnode['num_exp'] = cnumexp
            cnode['num_anno'] = cnumanno
            nodes += '{ data: %s},' % cnode
        edges = ''
        for cedge in termInfo['links']:
            edges += '{ data: %s },' % cedge

        return render_template('term_info_graph.html', term=term, nodes=nodes, edges=edges, node_link_url=url_for('.ontology_info', term=''))
    return "term %s not found" % term


@Site_Main_Flask_Obj.route('/download', methods=['POST', 'GET'])
def download():
    '''return a list of all the weekly database dump files available for download
    '''
    data_dir = os.path.join(current_app.root_path, 'data_dump')
    onlyfiles = [f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f)) and (f.endswith('.psql') or f.endswith('.txt'))]

    webPage = render_header(title='Download')
    # webPage += render_template('userinfo.html', userid=userid, name=name, username=username, desc=desc, email=email, total_annotations=total_annotations, total_experiments=total_experiments)
    webPage += '''
    <h1>dbBact snapshot download</h1>
    Files are PostgreSQL 9.5.10 binary dumps.<br>
    For locally installing the dbBact rest-api server, see documentation at: <a href='https://github.com/amnona/dbbact-server'>https://github.com/amnona/dbbact-server</a><br>
    <div id="annot-table-div" class="tab-pane active" role="tabpanel" aria-labelledby="annot-table-div-header" style="margin-top: 20px; margin-bottom: 20px;">
    <table id="annot-table" class="table responsive" style="width: 100%;">
    <thead>
      <tr>
        <th>File Name</th>
        <th>Size</th>
      </tr>
    </thead>
    <tbody>'''

    onlyfiles.sort(reverse=True)
    for cfile in onlyfiles:
        csize = os.stat(os.path.join(data_dir, cfile)).st_size
        if csize < 1000:
            csize_str = '%d Bytes' % csize
        elif csize < 1000000:
            csize_str = '%d KB' % (csize / 1000)
        elif csize < 1000000000:
            csize_str = '%d MB' % (csize / 1000000)
        else:
            csize_str = '%d GB' % (csize / 1000000000)

        webPage += "<tr><td><a href='%s'>" % url_for('.get_file', filename=cfile) + cfile + "</a></td><td>%s</tr>" % csize_str
    webPage += '''</tbody>
    </table>
    </div>'''

    webPage += render_template('footer.html')
    return webPage

    return str(onlyfiles)


@Site_Main_Flask_Obj.route('/faq_prepare_sequences', methods=['POST', 'GET'])
def faq_prepare_sequences():
    webpage = render_header(title='Preparing sequences for dbBact query')
    webpage += render_template('prepare_sequences.html')
    return webpage


@Site_Main_Flask_Obj.route('/get_file/<string:filename>', methods=['POST', 'GET'])
def get_file(filename):
    '''return a file from the data_dump directory.
    '''
    return send_from_directory('data_dump', filename, as_attachment=True)


def trim_primers_from_sequence(cseq, primers={'AGAGTTTGATC[AC]TGG[CT]TCAG': 'v1', 'CCTACGGG[ACGT][CGT]GC[AT][CG]CAG': 'v3', 'GTGCCAGC[AC]GCCGCGGTAA': 'v4'}, max_start=25, min_primer_len=10):
    '''Try to trim primers from a given sequence and return the trimmed sequence if successful

    assumes the sequence is not in dbBact.

    Looks for V1/V3/V4 forward primers

    Parameters
    ----------
    cseq: str
        sequence to trim ('ACGT')
    primers: dict of {primer(str): region_name(str)}, optional
        the primers to test for
    max_start: int, optional
        maximal start position for the primer (i.e. do not return if primer starts after position max_start)
    min_primer_len: int, optional
        trim primers to keep only min_primer_len last chars

    Returns
    -------
    trimmed: str or None
        the trimmed sequence (str) if region identified
        None if no region identified
    msg: str
        the information how the trimmed sequence was obtained
    '''
    # trim the primers if needed
    if min_primer_len is not None:
        new_primers = {}
        for k, v in primers.items():
            pos = len(k)
            numchars = 0
            newp = ''
            while True:
                if numchars >= min_primer_len:
                    break
                pos = pos - 1
                if pos < 0:
                    break
                if k[pos] != ']':
                    newp = k[pos] + newp
                    numchars += 1
                    continue
                while k[pos] != '[':
                    newp = k[pos] + newp
                    pos = pos - 1
                newp = k[pos] + newp
                numchars += 1
            new_primers[newp] = v
        primers = new_primers

    cseq = cseq.upper()

    # test for primers in the sequence
    for cprimer in primers.keys():
        ccseq = cseq[:max_start + len(cprimer)]
        match = re.search(cprimer, ccseq)
        if match is not None:
            trimmed = cseq[match.end():]
            return trimmed, 'found primer %s (region %s). trimmed first %d nucleotides.' % (cprimer, primers[cprimer], match.end())

    # didn't find primer, so let's check if it seems to start after a few nucleotides with a known region sequence
    for ctrim in range(10):
        trimmed = cseq[ctrim + 1:]
        if test_if_sequence_exists(trimmed):
            return trimmed, 'Left-trimmed first %d nucleotides' % (ctrim + 1)

    return None, 'no region identified'


def test_if_sequence_exists(seq, use_sequence_translator=True):
    '''test if a sequence exists in dbBact (sequencesTable)

    Parameters:
    -----------
    seq: str
        the sequence to search for ('ACGT')

    Returns:
    --------
    True if exists, False if doesn't exist
    '''
    rdata = {'sequence': seq, 'use_sequence_translator': use_sequence_translator}
    httpResTax = requests.get(dbbact_server_address + '/sequences/getid', json=rdata)
    if httpResTax.status_code == requests.codes.ok:
        match_ids = httpResTax.json().get('seqId')
        if match_ids is None:
            return False
        if len(match_ids) > 0:
            return True
    return False


def get_close_sequences(sequence, max_mismatches = 2):
    '''Get sequences that are close to the given sequence
    
    Parameters
    ----------
    sequence: str
        the sequence to search for ('ACGT')
    max_mismatches: int, optional
        the maximal number of mismatches to allow
        
    Returns
    -------
    list of {'sequence': str, 'seq_id':int, 'mismatches': int} (one entery per similar sequence)
    '''
    debug(2, 'get_close_sequences for sequence %s' % sequence)
    rdata = {'sequence': sequence, 'max_mismatches': max_mismatches}
    httpResTax = requests.get(dbbact_server_address + '/sequences/get_close_sequences', json=rdata)
    if httpResTax.status_code == requests.codes.ok:
        res = httpResTax.json()
        debug(2, 'Found %d close sequences' % len(res['similar_seqs']))
        debug(1, str(res['similar_seqs']))
        return res['similar_seqs']
    debug(7, 'error encountered when trying to get close sequences for sequence %s' % sequence)
    return []


def get_term_seq_scores(term, annotations=None, num_to_show=10):
    '''Calculate the F-score for all sequences in dbBact that are associated with the term

    Parameters
    ----------
    term: str
        the dbBact term
    annotations: list of dict or None, optional
        if not None, the dbBact annotations containing the term.
        if None, the function will query dbBact for the term annotations
    num_to_show: int, optional
        the max number of top positive/negative associations to plot

    Returns
    -------
    part of webpage with top positive/negative associated sequences
    '''
    # get all the child terms for the term of interest (so we can test each annotation positive/negative also for them)
    term = term.lower()
    res = requests.get(get_dbbact_server_address() + '/ontology/get_term_children', json={'term': term, 'only_annotated': 'true'})
    children = set(res.json()['terms'].values())
    debug(2, 'found %d children terms (including main term) with annotations for term %s' % (len(children), term))

    if annotations is None:
        res = requests.get(get_dbbact_server_address() + '/ontology/get_annotations', params={'term': term, 'get_children': 'true'})
        if res.status_code != 200:
            msg = 'error getting annotations for ontology term %s: %s' % (Markup.escape(term), res.content)
            debug(6, msg)
            return msg, msg
        annotations = res.json()['annotations']

    if len(annotations) == 0:
        debug(1, 'ontology term %s not found' % Markup.escape(term))
        return 'term not found', 'term not found'

    annotation_ids = [x['annotationid'] for x in annotations]

    # get the sequences for each annotation
    res = requests.get(get_dbbact_server_address() + '/annotations/get_list_sequences', json={'annotation_ids': annotation_ids})
    if res.status_code != 200:
        msg = 'error getting annotation sequences for ontology term %s: %s' % (Markup.escape(term), res.content)
        debug(6, msg)
        return msg, msg
    annotation_seqs = res.json()['annotation_seqs']

    # divide the annotations into positive correlation (COMMON/HIGH/DOMINANT) and negative correlation (LOW) with the term
    annotation_seqs_low = {}
    annotation_seqs_high = {}
    annotation_seqs_common = {}
    for cannotation in annotations:
        term_context = None
        cannotation_id = str(cannotation['annotationid'])
        details = cannotation['details']
        for cdetail in details:
            cterm = cdetail[1]
            # does this annotation detail contain the term or one of it's children?
            if cterm not in children:
                continue
            ctype = cdetail[0]
            if ctype == 'low':
                term_context = 'low'
                break
            elif ctype == 'high':
                term_context = 'high'
                break
            elif ctype == 'all':
                term_context = 'all'
            else:
                term_context = 'other'

        if term_context is None:
            debug(6, 'term %s does not appear in any context in annotation %s' % (term, cannotation_id))
            continue
        if term_context == 'low':
            annotation_seqs_low[cannotation_id] = annotation_seqs[cannotation_id]
        if term_context == 'high':
            annotation_seqs_high[cannotation_id] = annotation_seqs[cannotation_id]
        if term_context in ['high', 'all']:
            annotation_seqs_common[cannotation_id] = annotation_seqs[cannotation_id]

    seq_annotations_low = defaultdict(list)
    seq_annotations_high = defaultdict(list)
    seq_annotations_common = defaultdict(list)
    seqs = set()
    for cannotation_id, cseqlist in annotation_seqs_low.items():
        for cseq in cseqlist:
            seq_annotations_low[cseq].append(cannotation_id)
        seqs.update(cseqlist)
    for cannotation_id, cseqlist in annotation_seqs_high.items():
        for cseq in cseqlist:
            seq_annotations_high[cseq].append(cannotation_id)
        seqs.update(cseqlist)
    for cannotation_id, cseqlist in annotation_seqs_common.items():
        for cseq in cseqlist:
            seq_annotations_common[cseq].append(cannotation_id)
        seqs.update(cseqlist)

    # get the total number of annotations per sequence
    seqlist = list(seqs)
    res = requests.get(get_dbbact_server_address() + '/sequences/get_info', json={'seqids': seqlist})
    if res.status_code != 200:
        msg = 'failed getting number of annotations per sequence'
        debug(6, msg)
        return msg
    seq_info = res.json()['sequences']

    tot_seq_annotations = {}
    for idx, cdat in enumerate(seq_info):
        tot_seq_annotations[seqlist[idx]] = cdat['total_annotations']

    res = requests.get(get_dbbact_server_address() + '/ontology/get_term_stats', json={'terms': [term, '-' + term]})
    if res.status_code != 200:
        msg = 'failed getting term %s stats' % term
        debug(6, msg)
        return msg
    res = res.json()['term_info']
    if term not in res:
        msg = 'get_term_stats failed for term %s' % term
        debug(6, msg)
        return msg
    num_term_annotations_high = res[term]['total_annotations']
    num_term_annotations_low = res['-' + term]['total_annotations']

    precision_low = {}
    precision_high = {}
    precision_common = {}
    for cseq in seqs:
        precision_low[cseq] = len(seq_annotations_low[cseq]) / (tot_seq_annotations[cseq] + 1)
        precision_high[cseq] = len(seq_annotations_high[cseq]) / (tot_seq_annotations[cseq] + 1)
        precision_common[cseq] = len(seq_annotations_common[cseq]) / (tot_seq_annotations[cseq] + 1)

    recall_low = {}
    recall_high = {}
    recall_common = {}
    for cseq in seqs:
        recall_low[cseq] = len(seq_annotations_low[cseq]) / (num_term_annotations_low + 1)
        recall_high[cseq] = len(seq_annotations_high[cseq]) / (num_term_annotations_high + 1)
        recall_common[cseq] = len(seq_annotations_common[cseq]) / (num_term_annotations_high + 1)

    fscore_low = defaultdict(float)
    fscore_high = defaultdict(float)
    fscore_common = defaultdict(float)
    for cseq in seqs:
        if recall_low[cseq] + precision_low[cseq] > 0:
            # fscore_low[cseq] = 2 * (recall_low[cseq] * precision_low[cseq]) / (recall_low[cseq] + precision_low[cseq])
            fscore_low[cseq] = recall_low[cseq]
        if recall_high[cseq] + precision_high[cseq] > 0:
            # fscore_high[cseq] = 2 * (recall_high[cseq] * precision_high[cseq]) / (recall_high[cseq] + precision_high[cseq])
            fscore_high[cseq] = recall_high[cseq]
        if recall_common[cseq] + precision_common[cseq] > 0:
            fscore_common[cseq] = 2 * (recall_common[cseq] * precision_common[cseq]) / (recall_common[cseq] + precision_common[cseq])

    sfscore_high = {}
    sfscore_low = {}
    sfscore_common = {}
    for idx, cseq in enumerate(seqlist):
        sfscore_high[idx] = fscore_high[cseq]
        sfscore_low[idx] = fscore_low[cseq]
        sfscore_common[idx] = fscore_common[cseq]

    sorted_f_high = sorted(sfscore_high.items(), key=lambda kv: kv[1], reverse=True)
    sorted_f_low = sorted(sfscore_low.items(), key=lambda kv: kv[1], reverse=True)
    sorted_f_common = sorted(sfscore_common.items(), key=lambda kv: kv[1], reverse=True)

    seq_info_high = []
    for cpos in range(min([num_to_show, len(sorted_f_high)])):
        cidx, cscore = sorted_f_high[cpos]
        seq_info[cidx]['fscore_high'] = cscore
        seq_info_high.append(seq_info[cidx])

    seq_info_low = []
    for cpos in range(min([num_to_show, len(sorted_f_low)])):
        cidx, cscore = sorted_f_low[cpos]
        seq_info[cidx]['fscore_low'] = cscore
        seq_info_low.append(seq_info[cidx])

    seq_info_common = []
    for cpos in range(min([num_to_show, len(sorted_f_common)])):
        cidx, cscore = sorted_f_common[cpos]
        seq_info[cidx]['fscore_common'] = cscore
        seq_info_common.append(seq_info[cidx])

    return draw_highlow_sequences_info(seqinfo_high=seq_info_high, seqinfo_low=seq_info_low, seqinfo_common=seq_info_common)


def draw_highlow_sequences_info(seqinfo_high, seqinfo_low, seqinfo_common):
    '''write the table entries for each sequence (sequence, total counts etc.)
    '''
    # webPage = render_template('header.html')
    datalow = ''
    for cseqinfo in seqinfo_low:
        if cseqinfo['fscore_low'] == 0:
            continue
        cseqinfo['seq'] = cseqinfo['seq'].upper()
        datalow += "<tr>"
        datalow += '<td>' + cseqinfo['taxonomy'] + '</td>'
        datalow += '<td><a href=%s>%s</a></td>' % (url_for('.sequence_annotations', sequence=cseqinfo['seq']), Markup.escape(cseqinfo['seq']))
        datalow += '<td>' + '%f' % cseqinfo['fscore_low'] + '</td></tr>'

    datahigh = ''
    for cseqinfo in seqinfo_high:
        if cseqinfo['fscore_high'] == 0:
            continue
        cseqinfo['seq'] = cseqinfo['seq'].upper()
        datahigh += "<tr>"
        datahigh += '<td>' + cseqinfo['taxonomy'] + '</td>'
        datahigh += '<td><a href=%s>%s</a></td>' % (url_for('.sequence_annotations', sequence=cseqinfo['seq']), Markup.escape(cseqinfo['seq']))
        datahigh += '<td>' + '%f' % cseqinfo['fscore_high'] + '</td></tr>'

    datacommon = ''
    for cseqinfo in seqinfo_common:
        if cseqinfo['fscore_common'] == 0:
            continue
        cseqinfo['seq'] = cseqinfo['seq'].upper()
        datacommon += "<tr>"
        datacommon += '<td>' + cseqinfo['taxonomy'] + '</td>'
        datacommon += '<td><a href=%s>%s</a></td>' % (url_for('.sequence_annotations', sequence=cseqinfo['seq']), Markup.escape(cseqinfo['seq']))
        datacommon += '<td>' + '%f' % cseqinfo['fscore_common'] + '</td></tr>'

    webPage = render_template('seqlist-highlow.html', datalow=datalow, datahigh=datahigh, datacommon=datacommon)
    return webPage


@Site_Main_Flask_Obj.route('/interactive', methods=['POST', 'GET'])
def interactive():
    '''Interactive analysis of experiment data
    Used to generate the following:
    * interactive heatmap
    * f-score wordcloud
    * term enrichment analysis
    
    When submitted, goes to interactive_experiment_get_data'''
    webPage = render_header() + render_template('interactive.html') + render_template('footer.html')
    return webPage


def _load_experiment_table():
    '''load the table and metadata files using the session object

    we use the following session variables:
    table_format - the format of the table file (biom or qiime2)
    table_tmp_file_name - the name of the table file
    metadata_tmp_file_name - the name of the metadata file
    
    Returns
    -------
    err: str
        the error enounctered or empty if ok
    table : calour.AmpliconExperiment
    '''
    import calour as ca
    import tempfile
    import os

    table_format = session['table_format']
    table_tmp_file_name = session['table_tmp_file_name']
    metadata_tmp_file_name = session['metadata_tmp_file_name']
    repseqs_tmp_file_name = session.get('repseqs_tmp_file_name', None)

    if table_format == 'biom':
        table = ca.read_amplicon(table_tmp_file_name, metadata_tmp_file_name, min_reads=1000, normalize=10000)
    else:
        table = ca.read_qiime2(table_tmp_file_name, metadata_tmp_file_name, rep_seq_file=repseqs_tmp_file_name, min_reads=1000, normalize=10000)
        if len(table.feature_metadata.index.values[0]) == 32:
            return 'The table seems to contain hashed ids. Please also supply the rep-seqs file', table

    seq_len = len(table.feature_metadata.index.values[0])

    if seq_len == 32:
        return 'The table seems to contain hashed ids. Please also supply the rep-seqs file', table

    if seq_len < 100:
        return 'the table contains too short sequences (length=%d). Minimal length should be 100bp' % seq_len, table

    debug(2,'loaded table %s. sequence length is %d' % (table, seq_len))
    return '', table


@Site_Main_Flask_Obj.route('/interactive_experiment_get_data', methods=['POST', 'GET'])
def interactive_experiment_get_data():
    '''Store the uploaded files and redirect to analysis selection (to prevent resubmitting the form).
    Called when submitting the interactive() form.
    We store the following files locally as a temp file and store the tmp file name in the session object:
    table-file
    metadata-file
    experiment_name

    Then redirect to interactive_experiment_details
    '''
    import tempfile
    import os

    # get the form data
    form = request.form
    # get the name parameter from the form
    experiment_name = form['experiment_name']

    # get the two files (table and metadata)
    table_file = request.files['table-file']
    metadata_file = request.files['metadata-file']
    repseqs_file = request.files['repseqs-file']
    if table_file.filename.endswith('.biom'):
        table_format = 'biom'
    elif table_file.filename.endswith('.qza'):
        table_format = 'qza'
    else:
        debug(2, 'bad file name (%s). Dont know the file type' % table_file.filename)
        return 'table file format not recognized. Please select a .biom or .qza file'

    # save the uploaded files
    with tempfile.NamedTemporaryFile(suffix='.'+table_format, dir='./dbbact_website/tmp_files', delete=False) as table_file_name:
        table_file.save(table_file_name.name)
        table_tmp_file_name = table_file_name.name
 
    with tempfile.NamedTemporaryFile(suffix='.txt',dir='./dbbact_website/tmp_files',delete=False) as metadata_tmp_file:
        metadata_file.save(metadata_tmp_file.name)
        metadata_tmp_file_name = metadata_tmp_file.name

    repseqs_tmp_file_name = None
    if repseqs_file:
        with tempfile.NamedTemporaryFile(suffix='.qza',dir='./dbbact_website/tmp_files',delete=False) as repseqs_tmp_file:
            repseqs_file.save(repseqs_tmp_file.name)
            repseqs_tmp_file_name = repseqs_tmp_file.name

    # # clear the previous tmp files of the user session
    # try:
    #     if 'table_tmp_file_name' in session:
    #         os.remove(session['table_tmp_file_name'])
    #     if 'metadata_tmp_file_name' in session:
    #         os.remove(session['metadata_tmp_file_name'])
    #     if 'repseqs_tmp_file_name' in session:
    #         if session['repseqs_tmp_file_name']:
    #             os.remove(session['repseqs_tmp_file_name'])
    # except:
    #     pass

    # and set the new tmp file names in the session
    session['table_tmp_file_name'] = table_tmp_file_name
    session['metadata_tmp_file_name'] = metadata_tmp_file_name
    session['table_format'] = table_format
    session['experiment_name'] = experiment_name
    session['repseqs_tmp_file_name'] = repseqs_tmp_file_name

    return redirect(url_for('.interactive_experiment_details'))

@Site_Main_Flask_Obj.route('/interactive_experiment_get_example', methods=['POST', 'GET'])
def interactive_experiment_get_example():
    '''Prepare the analysis example dataset (when user presses the example link from the interactive webpage)
    Set the temp file names and redirect to interactive_experiment_details
    '''
    # and set the new tmp file names to point to the example files
    base_dir = './dbbact_website/static/analysis_example_data'
    session['table_tmp_file_name'] = os.path.join(base_dir, 'table.biom')
    session['metadata_tmp_file_name'] = os.path.join(base_dir, 'metadata.txt')
    session['table_format'] = 'biom'
    session['experiment_name'] = 'Chronic fatigue syndrome example'
    session['repseqs_tmp_file_name'] = None

    # and redirect
    return redirect(url_for('.interactive_experiment_details'))


@Site_Main_Flask_Obj.route('/interactive_experiment_details', methods=['POST', 'GET'])
def interactive_experiment_details():
    '''After loading the table and metadata, show the summary page for the experiment
    From this page the user can select the analysis to perform
    '''
    # load the table and metadata into a Calour AmpliconExperiment object
    import calour as ca

    table_format = session['table_format']
    table_tmp_file_name = session['table_tmp_file_name']
    metadata_tmp_file_name = session['metadata_tmp_file_name']
    experiment_name = session.get('experiment_name', 'NA')

    err, table = _load_experiment_table()
    if err:
        return render_header() + '<h1>Error encountered</h1><br>' + err

    # create a dict of field names and values in the table sample_metadata
    # keeping only fields with more than 1 value and less than all values
    field_data = {}
    for cfield in table.sample_metadata.columns:
        cvals = list(table.sample_metadata[cfield].unique())
        if len(cvals) <= 1:
            continue
        if len(cvals) == len(table.sample_metadata):
            continue
        cvals = [str(x) for x in cvals]
        cvals = sorted(cvals)
        field_data[cfield] = cvals

    return render_header()+render_template('interactive_exp_details.html', num_samples = len(table.sample_metadata), num_features = len(table.feature_metadata), metadata_fields=field_data, field_data=field_data, experiment_name=experiment_name) + render_template('footer.html')


@Site_Main_Flask_Obj.route('/wordcloud_analysis', methods=['POST', 'GET'])
def wordcloud_analysis():
    '''Draw a webpage with the experiment wordcloud
    Uses the experiment stored in session

    Parameters
    ----------
    method : str
        The method to use for weighting the ASVs. options:
        'none': usee all the ASVs with at least one read in the experiment
        'prevalence': use all the ASVs with at least 30% prevalence in the experiment
        'linear': weigh the per-ASV f-scores by the meam ASV prevalence in the experiment
        'log': weigh the per-ASV f-scores by the log of the mean ASV prevalence in the experiment
        'rank': weigh the per-ASV f-scores by the rank of the mean ASV prevalence in the experiment
    '''
    import matplotlib.pyplot as plt

    method = request.args.get('method', 'prevalence')

    err, table = _load_experiment_table()
    if err:
        return render_header() + '<h1>Error encountered</h1><br>' + err
    if method == 'prevalence':
        debug(3, 'filtering table by prevalence')
        table = table.filter_prevalence(0.3)
        freq_weighted = False
    else:
        return 'method %s not supported' % method

    dbc = dbbact_calour.dbbact.DBBact(dburl=get_dbbact_server_address(), test_version=False)
    debug(2,'drawing wordcloud')
    f = dbc.draw_wordcloud(table, freq_weighted=freq_weighted)
    figfile = BytesIO()
    f.savefig(figfile, format='png', bbox_inches='tight')
    figfile.seek(0)  # rewind to beginning of file
    figdata_png = base64.b64encode(figfile.getvalue())
    figfile.close()
    figdat_md5 = hashlib.md5(figdata_png).hexdigest()
    figdata_tmp_file_name = os.path.join('./dbbact_website/tmp_files', figdat_md5)
    f.savefig(figdata_tmp_file_name, format='png', bbox_inches='tight', dpi=300)
    if 'files_md5' not in session:
        session['files_md5'] = {}
    session['files_md5'][figdat_md5] = 'wordcloud.png'
    plt.close(f.figure)
    debug(2,'got the figure data')
    # table.filter_prevalence(0.1).export_html(output_file='./dbbact_website/templates/tmp_heatmap.html')
    # return render_template('tmp_heatmap.html')

    # create an html page with a download button that downloads the e.feature_metadata pandas dataframe
    return render_header() + render_template('wordcloud.html', wordcloudimage=urllib.parse.quote(figdata_png), figure_link='/download_file/%s' % figdat_md5) + render_template('footer.html')


@Site_Main_Flask_Obj.route('/interactive_metadata_submit', methods=['POST', 'GET'])
def interactive_metadata_submit():
    form = request.form
    # get the metadata field to use
    metadata_field = form['metadata_field']
    # get the metadata values to use for group1 and group2
    metadata_group1 = form['group1_values'].split(',')
    metadata_group2 = form['group2_values'].split(',')

    session['metadata_field'] = metadata_field
    session['metadata_group1'] = metadata_group1
    session['metadata_group2'] = metadata_group2

    debug(2, 'metadata_field=%s, metadata_group1=%s, metadata_group2=%s' % (metadata_field, metadata_group1, metadata_group2))
    # we redirect to prevent form resubmission on refresh
    return redirect(url_for('.interactive_metadata_submit2'))


@Site_Main_Flask_Obj.route('/interactive_metadata_submit2', methods=['POST', 'GET'])
def interactive_metadata_submit2():
    '''The differential abundance+term enrichment analysis
    assumes the metadata_field, metadata_group1 and metadata_group2 are already in the session (from interactive_metadata_submit)
    '''
    import matplotlib.pyplot as plt

    err, table = _load_experiment_table()
    if err:
        return render_header() + '<h1>Error encountered</h1><br>' + err


    metadata_field = session['metadata_field']
    metadata_group1 = session['metadata_group1']
    metadata_group2 = session['metadata_group2']
    # table_file_name = session['table_tmp_file_name']
    # metadata_file_name = session['metadata_tmp_file_name']
    # table_format = session['table_format']

    # import calour as ca
    # if table_format == 'biom':
    #     table = ca.read_amplicon(table_file_name, metadata_file_name, min_reads=1000, normalize=10000)
    # else:
    #     table = ca.read_qiime2(table_file_name, metadata_file_name, min_reads=1000, normalize=10000)

    # convert all the sample_metadata to strings (to overcome the form sending everything as strings)
    for cfield in table.sample_metadata.columns:
        table.sample_metadata[cfield] = table.sample_metadata[cfield].astype(str)

    dd=table.diff_abundance(metadata_field, metadata_group1, metadata_group2, random_seed=2023)
    debug(2,'got the diff abundance table %s' % dd)
    if dd is None:
        debug(2, 'no significant features found')
        return('No significant features found')

    if len(dd.feature_metadata) == 0:
        debug(2, 'no significant features found')
        return('No significant features found')

    debug(2, 'found %d differentially abundant features' % len(dd.feature_metadata))
    # save the differential abundance results to a tmp file and store it in the session
    if session.get('dd_tmp_file_name', None) is not None:
        try:
            os.remove(session['dd_tmp_file_name']+'.biom')
        except:
            pass
        try:
            os.remove(session['dd_tmp_file_name']+'_features.txt')
        except:
            pass
        try:
            os.remove(session['dd_tmp_file_name']+'_samples.txt')
        except:
            pass

    dd_tmp_file_name = tempfile.mktemp(suffix='_diff', dir='./dbbact_website/tmp_files')
    dd.save(dd_tmp_file_name)
    session['dd_tmp_file_name'] = dd_tmp_file_name
    debug(2, 'saved diff abundance results to tmp file %s' % dd_tmp_file_name)
    # save the diff abundance features table to a tmp file and store it in the session
    dd_features = dd.feature_metadata
    dd_md5 = hashlib.md5(dd_features.to_csv().encode()).hexdigest()
    dd_features_tmp_file_name = os.path.join('./dbbact_website/tmp_files', dd_md5)
    dd_features.to_csv(dd_features_tmp_file_name, sep='\t', index=False)
    if 'files_md5' not in session:
        session['files_md5'] = {}
    session['files_md5'][dd_md5] = 'diff_abundace_features.tsv'


    f,e = dd.plot_diff_abundance_enrichment()
    debug(2,'got the term enrichment plot')
    figfile = BytesIO()
    f.figure.savefig(figfile, format='png', bbox_inches='tight')
    figfile.seek(0)  # rewind to beginning of file
    figdata_png = base64.b64encode(figfile.getvalue())
    figfile.close()
    # save also as a tmp file for later download
    figdat_md5 = hashlib.md5(figdata_png).hexdigest()
    figdata_tmp_file_name = os.path.join('./dbbact_website/tmp_files', figdat_md5)
    f.figure.savefig(figdata_tmp_file_name, format='png', bbox_inches='tight', dpi=300)
    if 'files_md5' not in session:
        session['files_md5'] = {}
    session['files_md5'][figdat_md5] = 'enriched_terms.png'
    plt.close(f.figure)
    debug(2,'got the figure data')

    # create an html page with a download button that downloads the e.feature_metadata pandas dataframe
    group1_samples = len(table.sample_metadata[table.sample_metadata[metadata_field].isin(metadata_group1)])
    group2_samples = len(table.sample_metadata[table.sample_metadata[metadata_field].isin(metadata_group2)])
    group1_size = np.sum(dd.feature_metadata['_calour_stat']>0)
    group2_size = np.sum(dd.feature_metadata['_calour_stat']<0)
    # put the e.feature_metadata in a list where each entry is a list containing the index value, the _calour_effect_size and the _calour_qval
    # this is used to create the table in the html page

    MAX_TERMS_TO_SHOW = 20    
    showneg = np.sum(e.feature_metadata['odif']<0)
    if showneg > MAX_TERMS_TO_SHOW:
        showneg = MAX_TERMS_TO_SHOW
    showpos = np.sum(e.feature_metadata['odif']>0)
    if showpos > MAX_TERMS_TO_SHOW:
        showpos = MAX_TERMS_TO_SHOW

    # concatenate the frst showneg rows with the last showpos rows
    ee = e.feature_metadata[:showneg]
    ee = ee.append(e.feature_metadata[-showpos:])
    # reverse the order
    ee = ee.iloc[::-1]

    term_results = [[x[1]['term'], x[1]['odif'], x[1]['pvals']] for x in ee.iterrows()]

    # calculate the md5 checksum of the e.feature_metadata dataframe
    # this is used to check if the user downloaded the same file as the one we created
    # (to prevent malicious users from creating a fake file with the same name)
    md5sum = hashlib.md5(e.feature_metadata.to_csv(sep='\t', index=False).encode('utf-8')).hexdigest()
    if 'files_md5' not in session:
        session['files_md5'] = {}
    session['files_md5'][md5sum] = 'enriched_terms.tsv'
    # save the dataframe to a file named by the md5sum in the tmp_files folder
    enriched_tsv = os.path.join('./dbbact_website/tmp_files', md5sum)
    e.feature_metadata.to_csv(enriched_tsv, sep='\t', index=False)

    return render_header() + render_template('enriched-terms.html', toptermimage=urllib.parse.quote(figdata_png), 
                                             field=metadata_field, group1=metadata_group1, group2=metadata_group2,group1_samples=group1_samples,group2_samples=group2_samples,group1_size=group1_size,group2_size=group2_size,terms=term_results,
                                             enriched_terms_link='/download_file/%s' % md5sum, figure_link='/download_file/%s' % figdat_md5, diff_abundance_link='/download_file/%s' % dd_md5) + render_template('footer.html')


@Site_Main_Flask_Obj.route('/download_file/<string:md5sum>')
def download_file(md5sum):
    '''Download the enriched terms table

    Parameters
    ----------
    md5sum : str
        The md5sum of the file to download
    '''
    debug(4,'downloading file for hash %s' % md5sum)
    # if md5sum not in session.get('files_md5',{}).keys():
    #     debug(4, 'hash %s not in session' % md5sum)
    #     print(session.items())
    #     return 'Error: file hash %s not found' % md5sum
    # # get the file name from the session hash
    # visible_file_name = session['files_md5'][md5sum]
    # return the file
    # return send_file('tmp_files/%s' % md5sum, as_attachment=True, download_name=visible_file_name)
    return send_file('tmp_files/%s' % md5sum, as_attachment=True)


@Site_Main_Flask_Obj.route('/single_term_analysis/<string:term>', methods=['POST', 'GET'])
def single_term_analysis(term):
    '''Show a single term venn and f-score scatter

    NOTE: the function uses the session variable 'dd_tmp_file_name' to get the differential abundance results

    Parameters
    ----------
    term : str
        The term to show
    '''
    import matplotlib.pyplot as plt
    import calour as ca
    import calour_utils as cu

    metadata_field = session['metadata_field']
    metadata_group1 = session['metadata_group1']
    metadata_group2 = session['metadata_group2']

    # prepare the search term used for the dbBact query
    # we remove the * at the end and add a - at the beginning if it is 'LOWER IN'
    # and remove any spaces at the beginning or end
    search_term = term.lower()
    if search_term.startswith('lower in'):
        search_term = '-' + search_term[9:]
    if search_term.endswith('*'):
        search_term = search_term[:-1]
    search_term = search_term.strip()

    dd_tmp_file_name = session['dd_tmp_file_name']
    dd = ca.read_amplicon(dd_tmp_file_name+'.biom', dd_tmp_file_name+'_sample.txt', feature_metadata_file=dd_tmp_file_name+'_feature.txt', normalize=None, min_reads=0)
    debug(2, 'loaded diff abundance results from tmp file %s' % dd_tmp_file_name)

    dbc = dbbact_calour.dbbact.DBBact(dburl=get_dbbact_server_address(), test_version=False)
    group1_size = np.sum(dd.feature_metadata['_calour_stat']>0)
    group2_size = np.sum(dd.feature_metadata['_calour_stat']<0)
    f, overlaps, venn_pval = dbc.plot_term_venn_all(search_term,dd,max_size=np.max([group1_size,group2_size])*1.5, term_names=term)
    figfile = BytesIO()
    f.savefig(figfile, format='png', bbox_inches='tight')
    figfile.seek(0)  # rewind to beginning of file
    fig_venn_data_png = base64.b64encode(figfile.getvalue())
    figfile.close()
    # save for download button
    fig_venn_data_md5 = hashlib.md5(fig_venn_data_png).hexdigest()
    fig_venn_tmp_file_name = os.path.join('./dbbact_website/tmp_files', fig_venn_data_md5)
    f.figure.savefig(fig_venn_tmp_file_name, format='png', bbox_inches='tight', dpi=300)
    if 'files_md5' not in session:
        session['files_md5'] = {}
    session['files_md5'][fig_venn_data_md5] = 'venn_%s.png' % search_term
    plt.close(f.figure)
    debug(2,'got the venn figure data')

    # create the per-sample f-score plot
    # load the original biom table
    err, table = _load_experiment_table()
    if err:
        return render_header() + '<h1>Error encountered</h1><br>' + err

    # calculate the f-scores for the term in each sample. it is added as a new sample_metadata field
    newexp = dbc.get_term_sample_fscores(table, search_term, ignore_exp=True, transform=None)
    # plot the per-sample f-scores
    labels, vals, f = cu.plot_violin_category(newexp, session['metadata_field'], '_dbbact_fscore_'+search_term, show_labels=False)
    figfile = BytesIO()
    f.savefig(figfile, format='png', bbox_inches='tight')
    figfile.seek(0)  # rewind to beginning of file
    fig_fscores_data_png = base64.b64encode(figfile.getvalue())
    figfile.close()
    # save for download button
    fig_fscores_data_md5 = hashlib.md5(fig_fscores_data_png).hexdigest()
    fig_fscores_tmp_file_name = os.path.join('./dbbact_website/tmp_files', fig_fscores_data_md5)
    f.figure.savefig(fig_fscores_tmp_file_name, format='png', bbox_inches='tight', dpi=300)
    if 'files_md5' not in session:
        session['files_md5'] = {}
    session['files_md5'][fig_fscores_data_md5] = 'fscores_%s.png' % search_term
    plt.close(f.figure)
    sample_fscore_pval = scipy.stats.kruskal(*vals)[1]

    return render_header() + render_template('enriched-term-info.html', venn_image=urllib.parse.quote(fig_venn_data_png), sample_fscore_image=urllib.parse.quote(fig_fscores_data_png),
                                             term=term, field=metadata_field, group1=metadata_group1, group1_samples=group1_size, group2=metadata_group2, group2_samples=group2_size,
                                             venn_pval=venn_pval, sample_fscore_pval=sample_fscore_pval,
                                             venn_figure_link='/download_file/%s' % fig_venn_data_md5, fscores_figure_link='/download_file/%s' % fig_fscores_data_md5) + render_template('footer.html')


@Site_Main_Flask_Obj.route('/interactive_heatmap_submit', methods=['POST', 'GET'])
def interactive_heatmap_submit():
    import matplotlib.pyplot as plt

    form = request.form
    metadata_field = form['metadata_field_heatmap']
    if metadata_field.lower() == 'none':
        metadata_field = None
    debug(2,'generating heatmap for metadata field %s' % metadata_field)

    debug(2, 'loading table')
    err, table = _load_experiment_table()
    if err:
        return render_header() + '<h1>Error encountered</h1><br>' + err

    debug(2, 'clustering data')
    table = table.cluster_features(10)
    if metadata_field is not None:
        table = table.sort_samples(metadata_field)

    debug(2, 'drawing heatmap')
    table.export_html(output_file='./dbbact_website/templates/tmp_heatmap.html',sample_field=metadata_field)
    return render_template('tmp_heatmap.html')


@Site_Main_Flask_Obj.route('/get_sequences_stats', methods=['POST', 'GET'])
def get_sequences_stats():
    # TODO: Add taxonomy information (both list of exact matches and the RDP predicted taxonomy)
    alldat = request.get_json()
    if 'sequences' not in alldat:
        return 'Error: no sequences field provided in json'
    seqs = alldat.get('sequences')
    debug(2, 'getting annotations for %d sequences' % len(seqs))
    res = requests.get(get_dbbact_server_address() + '/sequences/get_fast_annotations',
                       json={'sequences': seqs})
    if res.status_code != 200:
        msg = 'error getting annotations for sequences : %s' % Markup.escape(res.content)
        debug(6, msg)
        return msg, msg
    res = res.json()
    annotations = res['annotations']
    seqannotations = res['seqannotations']
    if len(seqannotations) == 0:
        msg = 'None of the %d sequences were found in dbBact. Are these >100bp long 16S sequences?\nNote dbBact is populated mostly by EMP V4 (515F) amplicon sequences.' % len(seqs)
        debug(3, msg)
        return msg, msg
    term_info = res['term_info']

    debug(2,'got %d annotations' % len(annotations))

    # get the fscores for the sequences
    dbc = dbbact_calour.dbbact.DBBact(dburl=get_dbbact_server_address(), test_version=False)
    fscores, recall, precision, term_count, reduced_f = dbc.get_enrichment_score(annotations, seqannotations, term_info=term_info)

    res = requests.get(get_dbbact_server_address() + '/experiments/get_experiments_list')
    if res.status_code != 200:
        msg = 'error getting experiments list'
        debug(6, msg)
        return msg, msg
    res = res.json()
    exps = {}
    for cexp in res['explist']:
        cname = 'NA'
        for cdet in cexp[1]:
            if cdet[0] == 'name':
                cname = cdet[1]
                break
        exps[cexp[0]]=cname

    desc = []
    for cid,canno in annotations.items():
        ctext = getannotationstrings(canno, use_descriptions=False)
        cdesc = canno.get('description', None)
        cexp = exps.get(canno['expid'], 'NA')
        clink = 'https://dbbact.org/annotation_info/%s' % cid
        cdat = {'text': ctext, 'experiment': cexp, 'description': cdesc, 'link': clink}
        desc.append(cdat)

    # get the rdp taxonomy for the sequences
    taxonomies = []
    for cseq in seqs:
        rdata = {}
        rdata['sequence'] = cseq
        taxStr = "na"
        httpResTax = requests.get(dbbact_server_address + '/sequences/get_taxonomy_str', json=rdata)
        if httpResTax.status_code == requests.codes.ok:
            taxStr = httpResTax.json().get('taxonomy')
        else:
            debug(6, 'Error getting taxonomy for sequence %s: %s' % (cseq, httpResTax.content))
        taxonomies.append(taxStr)

    # Get the species and taxonomies based on 100% matching to the whole sequence database (i.e. silva)
    species = []
    for cseq in seqs:
        seq_species = []
        rdata = {}
        rdata['sequence'] = cseq
        httpResTax = requests.get(dbbact_server_address + '/sequences/get_whole_seq_taxonomy', json=rdata)
        if httpResTax.status_code == requests.codes.ok:
            seq_species = httpResTax.json().get('species')
        else:
            debug(6, 'Error getting whole sequence taxonomy for sequence %s: %s' % (cseq, httpResTax.content))
        species.append(seq_species)

    res = {'fscores': fscores, 'annotations': desc, 'rdp_taxonomy': taxonomies, 'matched_taxonomy': species}
    return res

@Site_Main_Flask_Obj.route('/sequences_fscores', methods=['POST', 'GET'])
def sequences_fscores():
    '''Return the term f-scores associated with a set of sequences
    
    json parameters:
        sequences: list of str
            the sequences (ACGT) to get the f-scores for
            NOTE: the sequences should be >100bp long 16S sequences from one of the supported dbBact regions
        score_method: str, optional
            'fscore' (default), 'recall', 'precision'
            the values to get the per-terms scores for
            precision: the fraction of per-sequence annotations that contain the term
            recall: the fraction of all annotations containing the term that are associated with the sequence
            fscore combines the recall and precision
        normalize: bool, optional
            True (default) to normalize the f-scores to a scale of 0-1 (max=1)
            False to return the raw f-scores
    Returns
    -------
    'fscores': dict of str:float
        the f-scores for each term (mean over all the sequences)
    '''
    alldat = request.get_json()
    if 'sequences' not in alldat:
        debug(6, 'Error: no sequences field provided in json')
        return 'Error: no sequences field provided in json', 400
    seqs = alldat.get('sequences')
    score_method = alldat.get('score_method', 'fscore')
    normalize = alldat.get('normalize', True)

    debug(2, 'getting annotations for %d sequences' % len(seqs))
    res = requests.get(get_dbbact_server_address() + '/sequences/get_fast_annotations',
                       json={'sequences': seqs})
    if res.status_code != 200:
        msg = 'error getting annotations for sequences : %s' % Markup.escape(res.content)
        debug(6, msg)
        return msg, 400
    res = res.json()

    annotations = res['annotations']
    term_info = res['term_info']
    seqannotations = res['seqannotations']
    if len(seqannotations) == 0:
        msg = 'None of the %d sequences were found in dbBact. Are these >100bp long 16S sequences?\nNote dbBact is populated mostly by EMP V4 (515F) amplicon sequences.' % len(seqs)
        debug(3, msg)
        return msg, 400

    debug(2,'got %d annotations' % len(annotations))

    # get the fscores for the sequences
    dbc = dbbact_calour.dbbact.DBBact(dburl=get_dbbact_server_address(), test_version=False)
    fscores, recall, precision, term_count, reduced_f = dbc.get_enrichment_score(annotations, seqannotations, term_info=term_info)

    if len(fscores) == 0:
        msg = 'no terms found for wordcloud'
        debug(4, msg)
        return msg, 400

    if score_method == 'fscore':
        score = fscores
    elif score_method == 'recall':
        score = recall
    elif score_method == 'precision':
        score = precision
    else:
        msg = 'score_method %s not supported. options are: fscore, recall, precision' % score_method
        debug(4, msg)
        return msg, 400

    # delete entries with emtpy key in scores
    new_scores = {}
    for ckey, cval in score.items():
        if ckey == '':
            continue
        new_scores[ckey] = cval
    score = new_scores

    # normalize the fractions to a scale max=1 if needed
    if normalize:
        new_scores = {}
        if score is not None:
            maxval = max(score.values())
            debug(1, 'normalizing score. maxval is %f' % maxval)
            for ckey, cval in score.items():
                new_scores[ckey] = score[ckey] / maxval
        score = new_scores

    return score
