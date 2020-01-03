#!/usr/bin/python

import json, requests
import vcf as pyvcf
import argparse
import sys

import requests
import json
import time
import sys

class EnsemblRestClient(object):
    def __init__(self, reqs_per_sec=15):
        #self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
        self.repeat = 0

    def perform_rest_action(self, server, endpoint, hdrs=None, parameters=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if parameters is None:
            parameters = {}
        data = None
        x=0
        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = requests.get(server + endpoint, headers=hdrs, params=parameters)
            request.raise_for_status()
            response = request.text
            if response:
                data = json.loads(response)
            self.req_count += 1

        except requests.exceptions.HTTPError as error:
            # check if we are being rate limited by the server
            if int(error.response.status_code) == 429:
                if 'Retry-After' in error.response.headers:
                    retry = error.response.headers['Retry-After']
                    time.sleep(float(retry))
                    data=self.perform_rest_action(server, endpoint, hdrs, parameters)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.response.status_code} Reason: {1.response.reason}\n'.format(server+endpoint, error))

        except requests.exceptions.ConnectionError as error:
            time.sleep(1)
            data=self.perform_rest_action(server, endpoint, hdrs, parameters)

        if data is None:
            self.repeat += 1
            if self.repeat <= 5:
                time.sleep(1)
                data=self.perform_rest_action(server, endpoint, hdrs, parameters)
            else:
                sys.stderr.write("Too many tries to connect to the Ensembl database")
        return data


#SERVER='http://grch37.rest.ensembl.org'
#ENDPOINT="/overlap/region/human/"+str(CHROM)+":"+str(POS)+"-"+str(POS)
#HEADERS={"Content-Type" : "application/json"}
#PARAMS={"feature": "transcript"}

#genes_data=EnsemblRestClient().perform_rest_action(SERVER, ENDPOINT,
#HEADERS, PARAMS)



parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('vcf', help='VCF file')
parser.add_argument('-o', '--offset', default=0, type=int, help='Offset [default: 0]')
parser.add_argument('-g', '--genombuild', required=True, type=str, help='Genome build version (default: %(default)s)')
parser.add_argument('-f', '--flank', default=200, type=int, help='Flank [default: 200]')
parser.add_argument('-m', '--mask',action='store_true')
args = parser.parse_args()
vcf = args.vcf
offset = args.offset
flank = args.flank
mask = args.mask

if args.genomebuild  == '38':
    server='http://rest.ensembl.org'
else:
    server='http://grch'+args.genomebuild+'.rest.ensembl.org'

def mask_seq( chr, start, end, strand, seq ):
    ext = "/overlap/region/human/"+str(chr)+":"+str(start)+"-"+str(end)+":"+str(strand)+"?feature=variation"
    headers={ "Content-Type" : "application/json"}

    data=EnsemblRestClient().perform_rest_action(server, ext, headers)
    masked_seq = seq

    for snp in data:
        snp_pos= snp['start']-start
        if str(strand) == '-1':
            snp_pos = len(masked_seq)-snp_pos-1
        masked_seq = masked_seq[:snp_pos]+'n'+masked_seq[snp_pos+1:];
    return( masked_seq )


def get_seq(chr, start, end, strand):
    ext = "/info/assembly/homo_sapiens/"+str(chr)
    headers={ "Content-Type" : "application/json"}

    data=EnsemblRestClient().perform_rest_action(server, ext, headers)
    chrlength = data['length']

    ext = "/sequence/region/human/"+str(chr)+":"+str(start)+"-"+str(end)+":"+str(strand)+""
    if start < 1:
        sys.stderr.write("\tFailed to get seq, because POS is too close to START of the chr\n")
        seq = False
    elif end > chrlength:
        sys.stderr.write("\tFailed to get seq, because ENDPOS is too close to END of the chr\n")
        seq = False
    else:
        data=EnsemblRestClient().perform_rest_action(server, ext, headers)

        seq = data['seq']
        if mask:
            masked_seq = mask_seq(chr, start, end, strand, seq )
        else:
            masked_seq = seq
    return( seq, masked_seq )

# def reverse_complement( seq ):
#     reverse_seq = seq[::-1]
#     complement = {'A':'T', 'C':'G', 'G':'C','T':'A'}
#     bases = list(reverse_seq)
#     bases = [complement.get(base,'n') for base in bases]
#     return( ''.join(bases) )

def alt_convert( record ):
    orientation = None
    remoteOrientation = None
    if record.INFO['SVTYPE'] == 'INS':
        orientation = False
        remoteOrientation = True
    if record.INFO['SVTYPE'] == 'DEL':
        orientation = False
        remoteOrientation = True
    elif record.INFO['SVTYPE'] == 'DUP':
        orientation = True
        remoteOrientation = False
    elif 'INV3' in record.INFO:
    	orientation = False
    	remoteOrientation = False
    elif 'INV5' in record.INFO:
        orientation = True
        remoteOrientation = True
    elif 'STRANDS' in record.INFO:
        if record.INFO['STRANDS'][0] == '++':
            orientation = False
            remoteOrientation = False
        elif record.INFO['STRANDS'][0] == '--':
            orientation = True
            remoteOrientation = True
        elif record.INFO['STRANDS'][0] == '-+':
            orientation = True
            remoteOrientation = False
        elif record.INFO['STRANDS'][0] == '+-':
            orientation = False
            remoteOrientation = True

    if 'CT' in record.INFO:
        if record.INFO['CT'] == '3to5':
            orientation = False
            remoteOrientation = True
        elif record.INFO['CT'] == '3to3':
            orientation = False
            remoteOrientation = False
        elif record.INFO['CT'] == '5to3':
            orientation = True
            remoteOrientation = False
        elif record.INFO['CT'] == '5to5':
            orientation = True
            remoteOrientation = True
    if orientation is None or remoteOrientation is None:
        sys.stderr.write( "\tFailed to convert orientation"+"\n" )
        record = False
    else:
        record.ALT = [ pyvcf.model._Breakend( record.CHROM, record.INFO['END'], orientation, remoteOrientation, record.REF, True ) ]
    return( record )


vcf_reader = pyvcf.Reader(open(vcf, 'r'))
for record in vcf_reader:
    if record.CHROM.startswith("#"):
        continue
    seq1 = ''
    seq2 = ''
    if isinstance( record.ALT[0], pyvcf.model._Substitution ):
        seq1a, masked_seq1a = '', ''
        seq2a, masked_seq2a = '', ''
        if offset != 0:
            seq1a, masked_seq1a = get_seq(record.CHROM, record.POS-offset, record.POS-1, 1)
            seq2a, masked_seq2a = get_seq(record.CHROM, record.POS+1, record.POS+offset, 1)
        seq1b, masked_seq1b = get_seq(record.CHROM, record.POS-offset-flank-1, record.POS-offset-1, 1)
        seq2b, masked_seq2b = get_seq(record.CHROM, record.POS+offset+len(record.REF), record.POS+len(record.REF)+offset+flank, 1)
    print(">"+record.CHROM+"_"+str(record.POS))
    print( masked_seq1b+"["+masked_seq1a+str(record.ALT[0]).lower()+masked_seq2a+"]"+masked_seq2b )

    #sys.stderr.write(">"+record.CHROM+"_"+str(record.POS)+"\n")
    #sys.stderr.write( seq1b+"["+seq1a+str(record.ALT[0]).lower()+seq2a+"]"+seq2b+"\n" )


# for record in vcf_reader:
#     sys.stderr.write( str(record) )
#     print( type(record.ALT[0]))
#     if isinstance(record.ALT[0], pyvcf.model._SV) or isinstance( record.ALT[0], pyvcf.model._Substitution ):
#         record = alt_convert( record )
#     if not record:
#         continue
#     seq1 = ''
#     seq2 = ''
#     if record.ALT[0].orientation:
#         seq1 = get_seq(record.CHROM, record.POS+offset, record.POS+offset+flank, -1)
#     else:
#         seq1 = get_seq(record.CHROM, record.POS-flank-offset, record.POS-offset, 1)
#     if not seq1:
#         continue
#     if record.ALT[0].remoteOrientation:
#         seq2 = get_seq(record.ALT[0].chr, record.ALT[0].pos+offset, record.ALT[0].pos+offset+flank, 1)
#     else:
#         seq2 = get_seq(record.ALT[0].chr, record.ALT[0].pos-offset-flank, record.ALT[0].pos-offset, -1)
#     if not seq2:
#         continue
#     sys.stderr.write("\tSucces\n")
#     print(">"+record.ID)
#     print( seq1+"[]"+seq2 )
